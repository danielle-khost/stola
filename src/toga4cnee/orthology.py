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


@dataclass
class QueryCluster:
    node_id: str
    chrom: str
    start: int
    end: int
    strands: set[str] = field(default_factory=set)
    refs: set[str] = field(default_factory=set)
    projection_count: int = 0
    query_cnee_id: str = ""
    orthology_category: str = ""

    @property
    def strand(self) -> str:
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
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise RuntimeError(
            "networkx is required for orthology classification. "
            "Activate the biopython conda environment or install networkx."
        ) from exc

    ref_ids = [interval.cnee_id for interval in parse_bed(cnee_bed)]
    projections = _read_aligned_projections(stats_path)
    clusters = _cluster_query_projections(projections)

    graph = nx.Graph()
    for ref_id in ref_ids:
        graph.add_node(_ref_node(ref_id), kind="ref", ref_id=ref_id)
    for cluster in clusters:
        graph.add_node(cluster.node_id, kind="query")
        for ref_id in cluster.refs:
            graph.add_edge(_ref_node(ref_id), cluster.node_id, weight=1)

    components = [set(component) for component in nx.connected_components(graph)]
    _assign_orthology_categories(components, clusters)
    orthology_count = _write_orthology_tsv(graph, clusters, ref_ids, out_orthology)
    bed_count = _write_query_bed(clusters, out_query_bed)
    return OrthologyResult(orthology_rows=orthology_count, query_bed_rows=bed_count)


def _read_aligned_projections(stats_path: str | Path) -> list[Projection]:
    projections: list[Projection] = []
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_id",
            "aligned_bases",
            "query_chrom",
            "query_strand",
            "query_aligned_start",
            "query_aligned_end",
            "query_projection_name",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            if int(row["aligned_bases"]) <= 0:
                continue
            if row["query_aligned_start"] == "." or row["query_aligned_end"] == ".":
                continue
            projections.append(
                Projection(
                    ref_id=row["cnee_id"],
                    query_chrom=row["query_chrom"],
                    query_start=int(row["query_aligned_start"]),
                    query_end=int(row["query_aligned_end"]),
                    query_strand=row["query_strand"],
                    projection_name=row["query_projection_name"],
                )
            )
    return projections


def _cluster_query_projections(projections: list[Projection]) -> list[QueryCluster]:
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
            current.projection_count += 1
        if current is not None:
            clusters.append(current)
    return clusters


def _assign_orthology_categories(
    components: list[set[str]], clusters: list[QueryCluster]
) -> None:
    cluster_by_node = {cluster.node_id: cluster for cluster in clusters}
    component_counter = 0
    for component in sorted(
        components,
        key=lambda nodes: sorted(nodes)[0],
    ):
        ref_nodes = sorted(node for node in component if str(node).startswith("ref:"))
        query_nodes = sorted(node for node in component if str(node).startswith("query:"))
        if not query_nodes:
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
    graph,
    clusters: list[QueryCluster],
    ref_ids: list[str],
    out_path: str | Path,
) -> int:
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
            partner_refs = sorted(
                _node_ref_id(node)
                for node in graph.neighbors(cluster.node_id)
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

        ref_with_query = {
            _node_ref_id(node)
            for cluster in clusters
            for node in graph.neighbors(cluster.node_id)
            if str(node).startswith("ref:")
        }
        for ref_id in ref_ids:
            if ref_id in ref_with_query:
                continue
            writer.writerow(
                {
                    "query_cnee_id": "NA",
                    "partner_ref_cnees": ref_id,
                    "orthology_category": "missing",
                }
            )
            count += 1
    return count


def _write_query_bed(clusters: list[QueryCluster], out_path: str | Path) -> int:
    count = 0
    with Path(out_path).open("w", encoding="utf-8") as handle:
        for cluster in sorted(clusters, key=lambda item: (item.chrom, item.start, item.end)):
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
    return sorted(
        query_nodes,
        key=lambda node: (
            cluster_by_node[node].chrom,
            cluster_by_node[node].start,
            cluster_by_node[node].end,
        ),
    )


def _ref_node(ref_id: str) -> str:
    return f"ref:{ref_id}"


def _node_ref_id(node: str) -> str:
    return node.split(":", 1)[1]


def _clean_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
