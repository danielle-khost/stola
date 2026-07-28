"""Graph-based CNEE orthology classification with synteny-verified pruning."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed
from .chain import iter_chain_blocks
from .index import load_index
from .models import BedInterval
from .sequence_compare import SequenceComparator, aligned_segments_for_interval
from .summary import is_solid_match
from .synteny import (
    BlockIndex,
    QueryBedRecord,
    QueryNeighborIndex,
    RefNeighborIndex,
    _aligned_bases,
)

ORTHOLOGY_COLUMNS = [
    "query_cnee_id",
    "partner_ref_cnees",
    "orthology_category",
]


@dataclass(frozen=True)
class Projection:
    ref_id: str
    ref_start: int
    ref_end: int
    chain_id: str
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
    index_path: str | Path,
    out_orthology: str | Path,
    out_query_bed: str | Path,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    synteny_check: bool = True,
    sequence_comparator: SequenceComparator | None = None,
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs.

    Builds a ref<->query graph from every aligned projection, then (unless
    synteny_check is False) verifies and prunes projections in ambiguous
    categories -- a non-solid single ref/query edge ("putative one2one(weak)"),
    or any edge inside a one2many/many2one/complex component -- using synteny
    evidence from the reference/query neighbor annotations plus the aligned
    chain itself. Query loci are re-clustered from the surviving projections
    before final categories are assigned, so pruning can split a cluster that
    was only held together by a since-removed projection.
    """
    if min_sequence_identity < 0 or min_sequence_identity > 1:
        raise ValueError("min_sequence_identity must be between 0 and 1")
    if max_fraction_gap < 0 or max_fraction_gap > 1:
        raise ValueError("max_fraction_gap must be between 0 and 1")
    try:
        import networkx as nx
    except ImportError as exc:
        raise RuntimeError(
            "networkx is required for orthology classification. "
            "Activate the biopython conda environment or install networkx."
        ) from exc

    ref_intervals = parse_bed(cnee_bed)
    ref_by_uid = {interval.cnee_uid: interval for interval in ref_intervals}
    ref_ids = [interval.cnee_uid for interval in ref_intervals]
    evidence = _read_projection_evidence(stats_path)
    refs_with_any_projection = {projection.ref_id for projection in evidence.projections}

    initial_clusters = _cluster_query_projections(evidence.projections)
    initial_graph = _build_orthology_graph(nx, ref_ids, initial_clusters)

    if synteny_check:
        surviving = _verify_and_prune_projections(
            nx,
            evidence.projections,
            initial_graph,
            initial_clusters,
            ref_by_uid,
            cnee_bed,
            index_path,
            sequence_comparator,
            min_sequence_identity,
            max_fraction_gap,
        )
        final_clusters = _cluster_query_projections(surviving)
        final_graph = _build_orthology_graph(nx, ref_ids, final_clusters)
    else:
        final_clusters = initial_clusters
        final_graph = initial_graph

    components = [set(component) for component in nx.connected_components(final_graph)]
    _assign_orthology_categories(
        components, final_clusters, final_graph, min_sequence_identity, max_fraction_gap
    )
    orthology_count = _write_orthology_tsv(
        final_clusters,
        ref_ids,
        evidence.refs_with_candidates,
        refs_with_any_projection,
        out_orthology,
    )
    bed_count = _write_query_bed(final_clusters, out_query_bed)
    return OrthologyResult(orthology_rows=orthology_count, query_bed_rows=bed_count)


def _read_projection_evidence(stats_path: str | Path) -> ProjectionEvidence:
    # Read candidate evidence and aligned query coordinates usable as graph edges.
    projections: list[Projection] = []
    refs_with_candidates: set[str] = set()
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_uid",
            "ref_start",
            "ref_end",
            "chain_id",
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
                    ref_start=int(row["ref_start"]),
                    ref_end=int(row["ref_end"]),
                    chain_id=row["chain_id"],
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
    # Merge overlapping query-side projections into provisional query CNEE
    # loci. Grouping by (chrom, strand) first ensures overlapping projections
    # on opposite strands are never merged into the same query locus.
    by_group: dict[tuple[str, str], list[Projection]] = defaultdict(list)
    for projection in projections:
        by_group[(projection.query_chrom, projection.query_strand)].append(projection)

    clusters: list[QueryCluster] = []
    for group_key in sorted(by_group):
        chrom = group_key[0]
        current: QueryCluster | None = None
        for projection in sorted(
            by_group[group_key],
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


def _edge_is_solid(
    edge_data: dict,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # An edge is solid only if at least one single projection meets both
    # criteria together; metrics are never blended across projections.
    return any(
        is_solid_match(sequence_identity, fraction_gap, min_sequence_identity, max_fraction_gap)
        for sequence_identity, fraction_gap in edge_data["projection_metrics"]
    )


# --------------------------------------------------------------------------
# Synteny verification: decide, per projection, whether it survives into the
# final graph. Only applied to "ambiguous" categories -- a non-solid lone
# ref/query edge, or any edge inside a one2many/many2one/complex component --
# a clean solid one2one edge is never second-guessed.
# --------------------------------------------------------------------------


def _components_needing_verification(
    nx, initial_graph, min_sequence_identity: float, max_fraction_gap: float
) -> dict[str, bool]:
    # Map every ref/query node to whether its component's edges need
    # verification: exempt only a lone (1 ref, 1 query) component whose sole
    # edge is already solid.
    needs_verification: dict[str, bool] = {}
    for component in nx.connected_components(initial_graph):
        ref_nodes = [node for node in component if str(node).startswith("ref:")]
        query_nodes = [node for node in component if str(node).startswith("query:")]
        if not ref_nodes or not query_nodes:
            continue
        if len(ref_nodes) == 1 and len(query_nodes) == 1:
            edge_data = initial_graph.edges[ref_nodes[0], query_nodes[0]]
            verify = not _edge_is_solid(edge_data, min_sequence_identity, max_fraction_gap)
        else:
            verify = True
        for node in component:
            needs_verification[node] = verify
    return needs_verification


def _query_neighbor_index_from_clusters(clusters: list[QueryCluster]) -> QueryNeighborIndex:
    by_chrom: dict[str, list[QueryBedRecord]] = defaultdict(list)
    for cluster in clusters:
        by_chrom[cluster.chrom].append(
            QueryBedRecord(
                chrom=cluster.chrom, start=cluster.start, end=cluster.end, query_id=cluster.node_id
            )
        )
    return QueryNeighborIndex.from_records(by_chrom)


def _verify_and_prune_projections(
    nx,
    projections: list[Projection],
    initial_graph,
    initial_clusters: list[QueryCluster],
    ref_by_uid: dict[str, BedInterval],
    cnee_bed: str | Path,
    index_path: str | Path,
    sequence_comparator: SequenceComparator | None,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> list[Projection]:
    needs_verification = _components_needing_verification(
        nx, initial_graph, min_sequence_identity, max_fraction_gap
    )
    if not any(needs_verification.values()):
        return list(projections)

    ref_neighbor_index = RefNeighborIndex.from_bed(cnee_bed)
    query_neighbor_index = _query_neighbor_index_from_clusters(initial_clusters)
    cluster_refs_by_query_node = {
        cluster.node_id: frozenset(cluster.refs) for cluster in initial_clusters
    }

    payload, records = load_index(index_path)
    chain_records_by_id = {record.chain_id: record for record in records}
    chain_path = payload["indexed_chain_path"]
    block_cache: dict[str, BlockIndex] = {}

    surviving: list[Projection] = []
    with Path(chain_path).open("rb") as chain_handle:
        for projection in projections:
            if not needs_verification.get(_ref_node(projection.ref_id), False):
                surviving.append(projection)
                continue
            if _verify_projection(
                projection,
                ref_by_uid,
                ref_neighbor_index,
                query_neighbor_index,
                cluster_refs_by_query_node,
                chain_records_by_id,
                chain_handle,
                block_cache,
                sequence_comparator,
                min_sequence_identity,
                max_fraction_gap,
            ):
                surviving.append(projection)
    return surviving


def _verify_projection(
    projection: Projection,
    ref_by_uid: dict[str, BedInterval],
    ref_neighbor_index: RefNeighborIndex,
    query_neighbor_index: QueryNeighborIndex,
    cluster_refs_by_query_node: dict[str, frozenset[str]],
    chain_records_by_id: dict,
    chain_handle,
    block_cache: dict[str, BlockIndex],
    sequence_comparator: SequenceComparator | None,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    """Score a projection against the 4 synteny-verification criteria.

    Required passes generalizes the two explicitly-specified cases (2-of-4
    normally, 1-of-3 when exactly one flank is missing) as roughly half of
    whichever criteria are actually applicable: `applicable // 2`, floored at 1.
    """
    applicable = 0
    passed = 0

    # Criterion 1: this specific projection is itself a solid match.
    applicable += 1
    if is_solid_match(
        projection.sequence_identity, projection.fraction_gap, min_sequence_identity, max_fraction_gap
    ):
        passed += 1

    # Reference and query coordinates are both normalized to plus-strand
    # genome coordinates, so the anchor that corresponds to the reference's
    # LEFT side flips to the query interval's END when strand is "-".
    forward = projection.query_strand == "+"
    if forward:
        anchors = {"left": (projection.query_start, "before"), "right": (projection.query_end, "after")}
    else:
        anchors = {"left": (projection.query_end, "after"), "right": (projection.query_start, "before")}

    ref_neighbors = {
        "left": ref_neighbor_index.neighbor(projection.ref_id, "left"),
        "right": ref_neighbor_index.neighbor(projection.ref_id, "right"),
    }

    # Criteria 2/3: left/right flank synteny match -- ref-neighbor and
    # query-neighbor (looked up against the *initial*, pre-pruning cluster
    # geometry, so this can't be circular) are partners in the initial graph.
    for side in ("left", "right"):
        ref_neighbor_id, _ = ref_neighbors[side]
        if ref_neighbor_id is None:
            continue
        anchor, direction = anchors[side]
        query_neighbor_node, _ = query_neighbor_index.neighbor(
            projection.query_chrom, anchor, direction
        )
        if query_neighbor_node is None:
            continue
        applicable += 1
        if ref_neighbor_id in cluster_refs_by_query_node.get(query_neighbor_node, frozenset()):
            passed += 1

    # Criterion 4: alignment of the connecting region is solid, required on
    # EVERY side that has a ref-neighbor at all (AND, not OR) -- extended to
    # cover the gap between the candidate and the neighbor plus the neighbor
    # itself, not just the neighbor's own annotation.
    sides_with_ref_neighbor = [side for side in ("left", "right") if ref_neighbors[side][0] is not None]
    if sides_with_ref_neighbor:
        applicable += 1
        record = chain_records_by_id[projection.chain_id]
        if projection.chain_id not in block_cache:
            chain_handle.seek(record.start_byte)
            chain_text = chain_handle.read(record.byte_length).decode("utf-8")
            block_cache[projection.chain_id] = BlockIndex.from_blocks(iter_chain_blocks(chain_text))
        blocks = block_cache[projection.chain_id]
        ref_chrom = ref_by_uid[projection.ref_id].chrom
        all_sides_solid = True
        for side in sides_with_ref_neighbor:
            ref_neighbor_id, _ = ref_neighbors[side]
            neighbor_interval = ref_by_uid[ref_neighbor_id]
            if side == "left":
                window_start, window_end = neighbor_interval.start, projection.ref_start
            else:
                window_start, window_end = projection.ref_end, neighbor_interval.end
            if not _window_is_solid(
                ref_chrom,
                window_start,
                window_end,
                record,
                blocks,
                sequence_comparator,
                min_sequence_identity,
                max_fraction_gap,
            ):
                all_sides_solid = False
                break
        if all_sides_solid:
            passed += 1

    # Relaxed to ~half of whatever's applicable (2-of-4 normally, 1-of-3 when
    # one flank is missing) rather than a strict majority -- the stricter bar
    # was pruning too many genuine one2many/many2one calls. Floored at 1 so a
    # fully isolated candidate (no neighbor on either side, applicable == 1)
    # still has to clear the lone "solid projection" criterion rather than
    # passing unconditionally.
    required = max(1, applicable // 2)
    return passed >= required


def _window_is_solid(
    ref_chrom: str,
    window_start: int,
    window_end: int,
    record,
    blocks: BlockIndex,
    sequence_comparator: SequenceComparator | None,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # "If there is no alignment somehow, treat that as failed" -- a window
    # with zero aligned bases (or a degenerate/empty window) never passes.
    if window_end <= window_start:
        return False
    overlapping_blocks = list(blocks.overlapping(window_start, window_end))
    aligned = _aligned_bases(overlapping_blocks, window_start, window_end)
    if aligned == 0:
        return False
    window_len = window_end - window_start
    fraction_gap = (window_len - aligned) / window_len
    sequence_identity = None
    if sequence_comparator is not None and overlapping_blocks:
        metrics = sequence_comparator.compare(
            ref_chrom,
            window_start,
            window_end,
            record.query_chrom,
            record.query_strand,
            aligned_segments_for_interval(overlapping_blocks, window_start, window_end, record.query_strand),
        )
        compared = int(metrics["seq_compared_bases"])
        if compared > 0:
            sequence_identity = metrics["seq_matches"] / compared
    return is_solid_match(sequence_identity, fraction_gap, min_sequence_identity, max_fraction_gap)


# --------------------------------------------------------------------------
# Category assignment and output.
# --------------------------------------------------------------------------


def _assign_orthology_categories(
    components: list[set[str]],
    clusters: list[QueryCluster],
    graph,
    min_sequence_identity: float,
    max_fraction_gap: float,
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
            edge_data = graph.edges[ref_nodes[0], query_nodes[0]]
            if _edge_is_solid(edge_data, min_sequence_identity, max_fraction_gap):
                category = "one2one"
            else:
                category = "one2one(weak)"
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
    clusters: list[QueryCluster],
    ref_ids: list[str],
    refs_with_candidates: set[str],
    refs_with_any_projection: set[str],
    out_path: str | Path,
) -> int:
    # Write compact query/ref/category rows, including lost/missing/
    # uncertain_loss references. Since clusters here are already the final
    # (post-verification) result, a cluster's own .refs set IS the correct
    # partner list -- no separate filtered-graph lookup needed.
    count = 0
    refs_with_final_query_edge: set[str] = set()

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
            partner_refs = sorted(cluster.refs)
            refs_with_final_query_edge.update(partner_refs)
            writer.writerow(
                {
                    "query_cnee_id": cluster.query_cnee_id,
                    "partner_ref_cnees": ",".join(partner_refs),
                    "orthology_category": cluster.orthology_category,
                }
            )
            count += 1

        for ref_id in ref_ids:
            if ref_id in refs_with_final_query_edge:
                continue
            if ref_id not in refs_with_candidates:
                category = "missing"
            elif ref_id not in refs_with_any_projection:
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
