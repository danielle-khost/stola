"""Graph-based CNEE orthology classification."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed

ORTHOLOGY_COLUMNS = [
    "query_cnee_id",
    "partner_ref_cnees",
    "orthology_category",
]


@dataclass(frozen=True)
class Projection:
    ref_id: str
    query_chrom: str
    query_start: int
    query_end: int
    query_strand: str
    projection_name: str
    sequence_identity: float | None
    fraction_gap: float


@dataclass(frozen=True)
class ProjectionEvidence:
    projections: list[Projection]
    refs_with_candidates: set[str]


@dataclass
class EdgeMetrics:
    projection_metrics: list[tuple[float | None, float]] = field(default_factory=list)

    def add_projection(self, projection: Projection) -> None:
        # Keep each projection's own (identity, gap) pair. Solidity is judged
        # per single projection, never blended across projections landing on
        # the same ref/query edge.
        self.projection_metrics.append(
            (projection.sequence_identity, projection.fraction_gap)
        )

    @property
    def projection_count(self) -> int:
        return len(self.projection_metrics)


@dataclass
class QueryCluster:
    node_id: str
    chrom: str
    start: int
    end: int
    strands: set[str] = field(default_factory=set)
    refs: set[str] = field(default_factory=set)
    edge_metrics: dict[str, EdgeMetrics] = field(default_factory=dict)
    projection_count: int = 0
    query_cnee_id: str = ""
    orthology_category: str = ""

    @property
    def strand(self) -> str:
        # Return a single cluster strand, or "." when merged projections disagree.
        return next(iter(self.strands)) if len(self.strands) == 1 else "."


@dataclass(frozen=True)
class OrthologyResult:
    orthology_rows: int
    query_bed_rows: int


def classify_orthology(
    stats_path: str | Path,
    cnee_bed: str | Path,
    out_orthology: str | Path,
    out_query_bed: str | Path,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    output_mode: str = "strict",
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs."""
    # Build graph evidence, filter solid edges, and write orthology calls.
    if min_sequence_identity < 0 or min_sequence_identity > 1:
        raise ValueError("min_sequence_identity must be between 0 and 1")
    if max_fraction_gap < 0 or max_fraction_gap > 1:
        raise ValueError("max_fraction_gap must be between 0 and 1")
    if output_mode not in {"all", "strict"}:
        raise ValueError("output_mode must be 'all' or 'strict'")
    try:
        import networkx as nx
    except ImportError as exc:
        raise RuntimeError(
            "networkx is required for orthology classification. "
            "Activate the biopython conda environment or install networkx."
        ) from exc

    ref_ids = [interval.cnee_uid for interval in parse_bed(cnee_bed)]
    evidence = _read_projection_evidence(stats_path)
    all_clusters = _cluster_query_projections(evidence.projections)
    all_graph = _build_orthology_graph(nx, ref_ids, all_clusters)
    all_solid_graph = _build_solid_graph(
        nx,
        all_graph,
        min_sequence_identity=min_sequence_identity,
        max_fraction_gap=max_fraction_gap,
    )
    if output_mode == "strict":
        projections = [
            projection
            for projection in evidence.projections
            if _projection_is_solid(
                projection,
                min_sequence_identity=min_sequence_identity,
                max_fraction_gap=max_fraction_gap,
            )
        ]
        clusters = _cluster_query_projections(projections)
        graph = _build_orthology_graph(nx, ref_ids, clusters)
    else:
        clusters = all_clusters
        refs_with_solid = _refs_with_query_edges(all_solid_graph, clusters)
        graph = _build_ref_filtered_graph(nx, all_graph, refs_with_solid)

    components = [set(component) for component in nx.connected_components(graph)]
    _assign_orthology_categories(components, clusters)
    orthology_count = _write_orthology_tsv(
        graph,
        all_graph,
        clusters,
        ref_ids,
        evidence.refs_with_candidates,
        out_orthology,
    )
    bed_count = _write_query_bed(clusters, out_query_bed)
    return OrthologyResult(orthology_rows=orthology_count, query_bed_rows=bed_count)


def _read_projection_evidence(stats_path: str | Path) -> ProjectionEvidence:
    # Read candidate evidence and aligned query coordinates usable as graph edges.
    projections: list[Projection] = []
    refs_with_candidates: set[str] = set()
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_uid",
            "aligned_bases",
            "query_chrom",
            "query_strand",
            "query_aligned_start",
            "query_aligned_end",
            "query_projection_name",
            "gap_fraction_within_chain_span",
            "seq_identity",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            refs_with_candidates.add(row["cnee_uid"])
            if int(row["aligned_bases"]) <= 0:
                continue
            if row["query_aligned_start"] == "." or row["query_aligned_end"] == ".":
                continue
            projections.append(
                Projection(
                    ref_id=row["cnee_uid"],
                    query_chrom=row["query_chrom"],
                    query_start=int(row["query_aligned_start"]),
                    query_end=int(row["query_aligned_end"]),
                    query_strand=row["query_strand"],
                    projection_name=row["query_projection_name"],
                    sequence_identity=_parse_optional_float(row["seq_identity"]),
                    fraction_gap=float(row["gap_fraction_within_chain_span"]),
                )
            )
    return ProjectionEvidence(
        projections=projections,
        refs_with_candidates=refs_with_candidates,
    )


def _cluster_query_projections(projections: list[Projection]) -> list[QueryCluster]:
    # Merge overlapping query-side projections into provisional query CNEE loci.
    by_chrom: dict[str, list[Projection]] = defaultdict(list)
    for projection in projections:
        by_chrom[projection.query_chrom].append(projection)

    clusters: list[QueryCluster] = []
    for chrom in sorted(by_chrom):
        current: QueryCluster | None = None
        for projection in sorted(
            by_chrom[chrom],
            key=lambda item: (item.query_start, item.query_end, item.ref_id),
        ):
            if current is None or projection.query_start >= current.end:
                if current is not None:
                    clusters.append(current)
                current = QueryCluster(
                    node_id=f"query:{len(clusters) + 1:06d}",
                    chrom=chrom,
                    start=projection.query_start,
                    end=projection.query_end,
                )
            current.end = max(current.end, projection.query_end)
            current.strands.add(projection.query_strand)
            current.refs.add(projection.ref_id)
            _add_cluster_edge_metrics(current, projection)
            current.projection_count += 1
        if current is not None:
            clusters.append(current)
    return clusters


def _add_cluster_edge_metrics(cluster: QueryCluster, projection: Projection) -> None:
    # Record this projection's own metrics on the ref/query edge.
    metrics = cluster.edge_metrics.setdefault(projection.ref_id, EdgeMetrics())
    metrics.add_projection(projection)


def _build_orthology_graph(nx, ref_ids: list[str], clusters: list[QueryCluster]):
    # Build a ref/query graph with sequence identity and gap fraction edge scores.
    graph = nx.Graph()
    for ref_id in ref_ids:
        graph.add_node(_ref_node(ref_id), kind="ref", ref_id=ref_id)
    for cluster in clusters:
        graph.add_node(cluster.node_id, kind="query")
        for ref_id in cluster.refs:
            metrics = cluster.edge_metrics[ref_id]
            graph.add_edge(
                _ref_node(ref_id),
                cluster.node_id,
                projection_metrics=metrics.projection_metrics,
                projection_count=metrics.projection_count,
            )
    return graph


def _build_solid_graph(
    nx,
    graph,
    min_sequence_identity: float,
    max_fraction_gap: float,
):
    # Copy all nodes but keep only edges that pass the current solid-match thresholds.
    solid_graph = nx.Graph()
    solid_graph.add_nodes_from(graph.nodes(data=True))
    for left, right, data in graph.edges(data=True):
        if _edge_is_solid(data, min_sequence_identity, max_fraction_gap):
            solid_graph.add_edge(left, right, **data)
    return solid_graph


def _build_ref_filtered_graph(nx, graph, ref_ids_to_keep: set[str]):
    # Keep all query edges for refs that have at least one solid projection.
    filtered_graph = nx.Graph()
    filtered_graph.add_nodes_from(graph.nodes(data=True))
    for left, right, data in graph.edges(data=True):
        ref_node = left if str(left).startswith("ref:") else right
        if _node_ref_id(ref_node) in ref_ids_to_keep:
            filtered_graph.add_edge(left, right, **data)
    return filtered_graph


def _refs_with_query_edges(graph, clusters: list[QueryCluster]) -> set[str]:
    # Collect reference IDs connected to any currently emitted query cluster.
    return {
        _node_ref_id(node)
        for cluster in clusters
        for node in graph.neighbors(cluster.node_id)
        if str(node).startswith("ref:")
    }


def _refs_with_any_query_edges(graph) -> set[str]:
    # Collect reference IDs that have any aligned query edge in a graph.
    refs: set[str] = set()
    for node, data in graph.nodes(data=True):
        if data.get("kind") != "ref":
            continue
        if any(str(neighbor).startswith("query:") for neighbor in graph.neighbors(node)):
            refs.add(_node_ref_id(node))
    return refs


def _solid(
    sequence_identity: float | None,
    fraction_gap: float,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # Apply the solid-match criteria to one (identity, gap) pair.
    if fraction_gap >= max_fraction_gap:
        return False
    if sequence_identity is None:
        return True
    return sequence_identity > min_sequence_identity


def _edge_is_solid(
    edge_data: dict,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # An edge is solid only if at least one single projection meets both
    # criteria together; metrics are never blended across projections.
    return any(
        _solid(sequence_identity, fraction_gap, min_sequence_identity, max_fraction_gap)
        for sequence_identity, fraction_gap in edge_data["projection_metrics"]
    )


def _projection_is_solid(
    projection: Projection,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # Apply the solid-match criteria to one candidate projection.
    return _solid(
        projection.sequence_identity,
        projection.fraction_gap,
        min_sequence_identity,
        max_fraction_gap,
    )


def _assign_orthology_categories(
    components: list[set[str]], clusters: list[QueryCluster]
) -> None:
    # Assign query IDs and broad orthology categories from graph components.
    cluster_by_node = {cluster.node_id: cluster for cluster in clusters}
    component_counter = 0
    for component in sorted(
        components,
        key=lambda nodes: sorted(nodes)[0],
    ):
        ref_nodes = sorted(node for node in component if str(node).startswith("ref:"))
        query_nodes = sorted(node for node in component if str(node).startswith("query:"))
        if not ref_nodes or not query_nodes:
            continue
        component_counter += 1
        ref_ids = [_node_ref_id(node) for node in ref_nodes]
        if len(ref_nodes) == 1 and len(query_nodes) == 1:
            category = "one2one"
            query_ids = [f"q_{_clean_id(ref_ids[0])}"]
        elif len(ref_nodes) == 1 and len(query_nodes) > 1:
            category = "one2many"
            query_ids = [
                f"q_{_clean_id(ref_ids[0])}_{index}"
                for index in range(1, len(query_nodes) + 1)
            ]
        elif len(ref_nodes) > 1 and len(query_nodes) == 1:
            category = "many2one"
            query_ids = [f"q_many2one_{component_counter:06d}"]
        else:
            category = "complex/many2many"
            query_ids = [
                f"q_complex_{component_counter:06d}_{index}"
                for index in range(1, len(query_nodes) + 1)
            ]

        for query_node, query_id in zip(
            _sort_query_nodes_by_coordinates(query_nodes, cluster_by_node), query_ids
        ):
            cluster = cluster_by_node[query_node]
            cluster.query_cnee_id = query_id
            cluster.orthology_category = category


def _write_orthology_tsv(
    solid_graph,
    all_graph,
    clusters: list[QueryCluster],
    ref_ids: list[str],
    refs_with_candidates: set[str],
    out_path: str | Path,
) -> int:
    # Write compact query/ref/category rows, including lost/missing references.
    count = 0

    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ORTHOLOGY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for cluster in sorted(clusters, key=lambda item: (item.chrom, item.start, item.end)):
            if not cluster.query_cnee_id:
                continue
            partner_refs = sorted(
                _node_ref_id(node)
                for node in solid_graph.neighbors(cluster.node_id)
                if str(node).startswith("ref:")
            )
            writer.writerow(
                {
                    "query_cnee_id": cluster.query_cnee_id,
                    "partner_ref_cnees": ",".join(partner_refs),
                    "orthology_category": cluster.orthology_category,
                }
            )
            count += 1

        refs_with_solid_query = _refs_with_query_edges(solid_graph, clusters)
        refs_with_aligned_query = _refs_with_any_query_edges(all_graph)
        for ref_id in ref_ids:
            if ref_id in refs_with_solid_query:
                continue
            if ref_id not in refs_with_candidates:
                category = "missing"
            elif ref_id not in refs_with_aligned_query:
                category = "lost"
            else:
                category = "uncertain_loss"
            writer.writerow(
                {
                    "query_cnee_id": "NA",
                    "partner_ref_cnees": ref_id,
                    "orthology_category": category,
                }
            )
            count += 1
    return count


def _write_query_bed(clusters: list[QueryCluster], out_path: str | Path) -> int:
    # Write graph-assigned query CNEE loci as BED6 rows.
    count = 0
    with Path(out_path).open("w", encoding="utf-8") as handle:
        for cluster in sorted(clusters, key=lambda item: (item.chrom, item.start, item.end)):
            if not cluster.query_cnee_id:
                continue
            handle.write(
                "\t".join(
                    [
                        cluster.chrom,
                        str(cluster.start),
                        str(cluster.end),
                        cluster.query_cnee_id,
                        "0",
                        cluster.strand,
                    ]
                )
                + "\n"
            )
            count += 1
    return count


def _sort_query_nodes_by_coordinates(
    query_nodes: list[str], cluster_by_node: dict[str, QueryCluster]
) -> list[str]:
    # Order query graph nodes by genomic coordinates for stable ID assignment.
    return sorted(
        query_nodes,
        key=lambda node: (
            cluster_by_node[node].chrom,
            cluster_by_node[node].start,
            cluster_by_node[node].end,
        ),
    )


def _ref_node(ref_id: str) -> str:
    # Build the graph node name used for a reference CNEE.
    return f"ref:{ref_id}"


def _node_ref_id(node: str) -> str:
    # Recover the original reference CNEE ID from a graph node name.
    return node.split(":", 1)[1]


def _clean_id(value: str) -> str:
    # Make a reference ID safe to embed in generated query CNEE IDs.
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _parse_optional_float(value: str) -> float | None:
    # Convert TSV numeric fields to floats while preserving missing values.
    if value == "." or value == "":
        return None
    return float(value)
