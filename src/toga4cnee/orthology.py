"""Graph-based CNEE orthology classification with synteny-verified pruning."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed
from .sequence_compare import SequenceComparator
from .summary import is_solid_match
from .synteny import QueryBedRecord, QueryNeighborIndex, RefNeighborIndex

ORTHOLOGY_COLUMNS = [
    "query_cnee_id",
    "partner_ref_cnees",
    "orthology_category",
]

CULLED_PROJECTION_COLUMNS = [
    "projection_name",
    "ref_cnee_id",
    "initial_orthology_category",
    "query_chrom",
    "query_start",
    "query_end",
    "sequence_identity",
    "gap_fraction",
]

# Verification criteria v2 (see synteny_verification.md): retain a projection
# only if it clears BOTH an align_score and a neighbor-ortholog-fraction bar.
DEFAULT_MIN_ALIGN_SCORE = 0.5
DEFAULT_MIN_NEIGHBOR_FRACTION = 0.5
DEFAULT_VERIFICATION_NUM_NEIGHBORS = 5


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


@dataclass(frozen=True)
class CulledProjection:
    projection: Projection
    initial_orthology_category: str


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
    culled_projection_rows: int


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
    min_align_score: float = DEFAULT_MIN_ALIGN_SCORE,
    min_neighbor_fraction: float = DEFAULT_MIN_NEIGHBOR_FRACTION,
    verification_num_neighbors: int = DEFAULT_VERIFICATION_NUM_NEIGHBORS,
    out_culled_projections: str | Path | None = None,
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs.

    Builds a ref<->query graph from every aligned projection, then (unless
    synteny_check is False) verifies and prunes projections in ambiguous
    categories -- a non-solid single ref/query edge (putative one2one(weak)),
    or any edge inside a one2many/many2one/complex component. "lost" refs
    (zero aligned query bases) are never part of this graph at all, so they
    are never subject to verification -- that category needs different
    criteria entirely and is deliberately left alone for now.

    A projection survives verification only if BOTH:
      - align_score (its own sequence identity relative to the best-scoring
        projection anywhere in its connected component) exceeds
        min_align_score, AND
      - neighbor_fraction (the fraction of its nearest verification_num_neighbors
        *trusted* reference-CNEE neighbors per flank -- i.e. refs that are
        themselves a solid one2one and so can't be pruned -- whose own query
        ortholog falls among this projection's nearest verification_num_neighbors
        trusted query neighbors) exceeds min_neighbor_fraction. A projection
        with no trusted neighbor on either flank fails this criterion
        outright (no synteny evidence to confirm it).

    Query loci are re-clustered from the surviving projections before final
    categories are assigned, so pruning can split a cluster that was only
    held together by a since-removed projection. A ref left with no
    surviving query edge (but that did have some projection) is written to
    out_orthology as "uncertain_loss"; if out_culled_projections is given,
    every individual pruned projection is also logged there for tracking.

    Any surviving single ref/query edge is reported as plain "one2one" --
    once verification has run, a weak edge that survived is just as trusted
    as one that was always solid, so the one2one(weak) label is dropped from
    final output entirely (it still shows up as an *initial*, pre-pruning
    category in out_culled_projections).
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
    ref_ids = [interval.cnee_uid for interval in ref_intervals]
    evidence = _read_projection_evidence(stats_path)
    refs_with_any_projection = {projection.ref_id for projection in evidence.projections}

    initial_clusters = _cluster_query_projections(evidence.projections)
    initial_graph = _build_orthology_graph(nx, ref_ids, initial_clusters)

    culled: list[CulledProjection] = []
    if synteny_check:
        surviving, culled = _verify_and_prune_projections(
            nx,
            evidence.projections,
            initial_graph,
            initial_clusters,
            cnee_bed,
            min_sequence_identity,
            max_fraction_gap,
            min_align_score,
            min_neighbor_fraction,
            verification_num_neighbors,
        )
        final_clusters = _cluster_query_projections(surviving)
        final_graph = _build_orthology_graph(nx, ref_ids, final_clusters)
    else:
        final_clusters = initial_clusters
        final_graph = initial_graph

    components = [set(component) for component in nx.connected_components(final_graph)]
    _assign_orthology_categories(components, final_clusters)
    orthology_count = _write_orthology_tsv(
        final_clusters,
        ref_ids,
        evidence.refs_with_candidates,
        refs_with_any_projection,
        out_orthology,
    )
    bed_count = _write_query_bed(final_clusters, out_query_bed)
    if out_culled_projections is not None:
        _write_culled_projections_tsv(culled, out_culled_projections)
    return OrthologyResult(
        orthology_rows=orthology_count,
        query_bed_rows=bed_count,
        culled_projection_rows=len(culled),
    )


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
# Synteny verification (v2): decide, per projection, whether it survives into
# the final graph. Only applied to "ambiguous" categories -- a non-solid lone
# ref/query edge, or any edge inside a one2many/many2one/complex component --
# a clean solid one2one edge is never second-guessed, and "lost" refs (no
# aligned query bases at all, so no Projection object exists for them) never
# enter this graph in the first place.
# --------------------------------------------------------------------------


def _query_neighbor_index_from_clusters(clusters: list[QueryCluster]) -> QueryNeighborIndex:
    by_chrom: dict[str, list[QueryBedRecord]] = defaultdict(list)
    for cluster in clusters:
        by_chrom[cluster.chrom].append(
            QueryBedRecord(
                chrom=cluster.chrom, start=cluster.start, end=cluster.end, query_id=cluster.node_id
            )
        )
    return QueryNeighborIndex.from_records(by_chrom)


def _initial_verification_state(
    nx, initial_graph, min_sequence_identity: float, max_fraction_gap: float
) -> tuple[dict[str, str], set[str], dict[str, int]]:
    """Categorize every ref from the raw, pre-pruning graph shape.

    Returns (initial_category_by_ref, trusted_ref_ids, component_id_by_ref):
    - initial_category_by_ref: the would-be category for every ref that needs
      verification (one2one(weak)/one2many/many2one/complex/many2many).
      Absence from this dict means the ref doesn't need verification at all
      (a solid lone one2one edge, or a ref with no projections).
    - trusted_ref_ids: refs NOT needing verification, i.e. already a solid
      lone one2one edge. Per synteny_verification.md v2, only these are
      eligible as synteny anchors for scoring *other* candidates -- a ref
      that might itself get pruned can't vouch for one still being verified.
    - component_id_by_ref: every ref's connected-component id in the initial
      graph, so a projection's siblings (for align_score) can be gathered.
    """
    initial_category_by_ref: dict[str, str] = {}
    trusted_ref_ids: set[str] = set()
    component_id_by_ref: dict[str, int] = {}
    for component_id, component in enumerate(nx.connected_components(initial_graph)):
        ref_nodes = [node for node in component if str(node).startswith("ref:")]
        query_nodes = [node for node in component if str(node).startswith("query:")]
        if not ref_nodes or not query_nodes:
            continue
        for node in ref_nodes:
            component_id_by_ref[_node_ref_id(node)] = component_id
        if len(ref_nodes) == 1 and len(query_nodes) == 1:
            edge_data = initial_graph.edges[ref_nodes[0], query_nodes[0]]
            if _edge_is_solid(edge_data, min_sequence_identity, max_fraction_gap):
                trusted_ref_ids.add(_node_ref_id(ref_nodes[0]))
            else:
                initial_category_by_ref[_node_ref_id(ref_nodes[0])] = "one2one(weak)"
        elif len(ref_nodes) == 1:
            initial_category_by_ref[_node_ref_id(ref_nodes[0])] = "one2many"
        elif len(query_nodes) == 1:
            for node in ref_nodes:
                initial_category_by_ref[_node_ref_id(node)] = "many2one"
        else:
            for node in ref_nodes:
                initial_category_by_ref[_node_ref_id(node)] = "complex/many2many"
    return initial_category_by_ref, trusted_ref_ids, component_id_by_ref


def _trusted_ref_neighbor_chain(
    ref_neighbor_index: RefNeighborIndex,
    trusted_ref_ids: set[str],
    ref_id: str,
    side: str,
    count: int,
) -> list[str]:
    # Walk outward along the reference genome from ref_id, skipping any CNEE
    # that isn't itself a trusted anchor, collecting up to `count` ids.
    chrom = ref_neighbor_index.chrom_by_uid.get(ref_id)
    if chrom is None:
        return []
    intervals = ref_neighbor_index.intervals_by_chrom[chrom]
    position = ref_neighbor_index.position_by_uid[ref_id]
    step = -1 if side == "left" else 1
    found: list[str] = []
    index = position + step
    while 0 <= index < len(intervals) and len(found) < count:
        candidate_id = intervals[index].cnee_uid
        if candidate_id in trusted_ref_ids:
            found.append(candidate_id)
        index += step
    return found


def _trusted_query_neighbor_chain(
    trusted_query_index: QueryNeighborIndex,
    trusted_query_cluster_by_node_id: dict[str, QueryCluster],
    chrom: str,
    anchor: int,
    direction: str,
    count: int,
) -> set[str]:
    # Walk outward along the query genome from anchor, through only the
    # trusted-only index, re-anchoring on each found cluster's far edge.
    found: set[str] = set()
    current = anchor
    for _ in range(count):
        node_id, _distance = trusted_query_index.neighbor(chrom, current, direction)
        if node_id is None:
            break
        found.add(node_id)
        cluster = trusted_query_cluster_by_node_id[node_id]
        current = cluster.start if direction == "before" else cluster.end
    return found


def _score_projection(
    projection: Projection,
    trusted_ref_ids: set[str],
    ref_neighbor_index: RefNeighborIndex,
    trusted_query_node_by_ref: dict[str, str],
    trusted_query_index: QueryNeighborIndex,
    trusted_query_cluster_by_node_id: dict[str, QueryCluster],
    max_identity_by_component: dict[int, float],
    component_id_by_ref: dict[str, int],
    num_neighbors: int,
) -> tuple[float, float | None]:
    """Score a single projection's align_score and neighbor_fraction.

    align_score: this projection's own sequence identity relative to the
    best-scoring projection anywhere in its connected component (including
    itself, so the component's best projection always lands at exactly 1.0;
    a lone one2one(weak) projection is trivially its own component max, so it
    always scores 1.0 too).

    neighbor_fraction: of the nearest `num_neighbors` trusted reference-CNEE
    neighbors on EACH flank (so up to 2 * num_neighbors checked total), the
    fraction whose own trusted query ortholog falls among this projection's
    own nearest `num_neighbors` trusted query neighbors. None if neither
    flank has any trusted neighbor at all (no synteny evidence either way).
    """
    component_id = component_id_by_ref[projection.ref_id]
    max_identity = max_identity_by_component.get(component_id, 0.0)
    if max_identity <= 0:
        align_score = 1.0
    else:
        align_score = (projection.sequence_identity or 0.0) / max_identity

    # Reference and query coordinates are both normalized to plus-strand
    # genome coordinates, so the anchor that corresponds to the reference's
    # LEFT side flips to the query interval's END when strand is "-".
    forward = projection.query_strand == "+"
    if forward:
        anchors = {"left": (projection.query_start, "before"), "right": (projection.query_end, "after")}
    else:
        anchors = {"left": (projection.query_end, "after"), "right": (projection.query_start, "before")}

    total_checked = 0
    total_matched = 0
    for side in ("left", "right"):
        ref_neighbors = _trusted_ref_neighbor_chain(
            ref_neighbor_index, trusted_ref_ids, projection.ref_id, side, num_neighbors
        )
        if not ref_neighbors:
            continue
        anchor, direction = anchors[side]
        query_neighbor_ids = _trusted_query_neighbor_chain(
            trusted_query_index,
            trusted_query_cluster_by_node_id,
            projection.query_chrom,
            anchor,
            direction,
            num_neighbors,
        )
        for neighbor_ref_id in ref_neighbors:
            total_checked += 1
            neighbor_query_node = trusted_query_node_by_ref.get(neighbor_ref_id)
            if neighbor_query_node is not None and neighbor_query_node in query_neighbor_ids:
                total_matched += 1

    neighbor_fraction = (total_matched / total_checked) if total_checked else None
    return align_score, neighbor_fraction


def _verify_and_prune_projections(
    nx,
    projections: list[Projection],
    initial_graph,
    initial_clusters: list[QueryCluster],
    cnee_bed: str | Path,
    min_sequence_identity: float,
    max_fraction_gap: float,
    min_align_score: float,
    min_neighbor_fraction: float,
    num_neighbors: int,
) -> tuple[list[Projection], list[CulledProjection]]:
    initial_category_by_ref, trusted_ref_ids, component_id_by_ref = _initial_verification_state(
        nx, initial_graph, min_sequence_identity, max_fraction_gap
    )
    if not initial_category_by_ref:
        return list(projections), []

    # Every trusted ref is, by definition, a solid lone one2one -- exactly
    # one query cluster, so this lookup is never ambiguous for trusted refs.
    query_node_by_ref: dict[str, str] = {}
    for cluster in initial_clusters:
        for ref_id in cluster.refs:
            query_node_by_ref.setdefault(ref_id, cluster.node_id)
    trusted_query_node_by_ref = {
        ref_id: query_node_by_ref[ref_id] for ref_id in trusted_ref_ids if ref_id in query_node_by_ref
    }
    trusted_cluster_ids = set(trusted_query_node_by_ref.values())
    trusted_query_clusters = [c for c in initial_clusters if c.node_id in trusted_cluster_ids]
    trusted_query_index = _query_neighbor_index_from_clusters(trusted_query_clusters)
    trusted_query_cluster_by_node_id = {c.node_id: c for c in trusted_query_clusters}

    ref_neighbor_index = RefNeighborIndex.from_bed(cnee_bed)

    projections_by_component: dict[int, list[Projection]] = defaultdict(list)
    for projection in projections:
        component_id = component_id_by_ref.get(projection.ref_id)
        if component_id is not None:
            projections_by_component[component_id].append(projection)
    max_identity_by_component = {
        component_id: max((p.sequence_identity or 0.0) for p in comp_projections)
        for component_id, comp_projections in projections_by_component.items()
    }

    surviving: list[Projection] = []
    culled: list[CulledProjection] = []
    for projection in projections:
        initial_category = initial_category_by_ref.get(projection.ref_id)
        if initial_category is None:
            surviving.append(projection)
            continue
        align_score, neighbor_fraction = _score_projection(
            projection,
            trusted_ref_ids,
            ref_neighbor_index,
            trusted_query_node_by_ref,
            trusted_query_index,
            trusted_query_cluster_by_node_id,
            max_identity_by_component,
            component_id_by_ref,
            num_neighbors,
        )
        passes = (
            align_score > min_align_score
            and neighbor_fraction is not None
            and neighbor_fraction > min_neighbor_fraction
        )
        if passes:
            surviving.append(projection)
        else:
            culled.append(
                CulledProjection(projection=projection, initial_orthology_category=initial_category)
            )
    return surviving, culled


# --------------------------------------------------------------------------
# Category assignment and output.
# --------------------------------------------------------------------------


def _assign_orthology_categories(
    components: list[set[str]],
    clusters: list[QueryCluster],
) -> None:
    # Assign query IDs and broad orthology categories from graph components.
    # A surviving lone ref/query edge is always "one2one" -- if it started
    # weak, it only got here by passing synteny verification, so there's
    # nothing left to distinguish it from an edge that was always solid.
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


def _write_culled_projections_tsv(culled: list[CulledProjection], out_path: str | Path) -> int:
    # Log every individual projection synteny verification pruned, for
    # tracking -- these refs end up "uncertain_loss" in out_orthology (unless
    # a sibling projection of the same ref survived), but this is the only
    # place the specific culled projection(s) themselves are recorded.
    count = 0
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CULLED_PROJECTION_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for item in sorted(
            culled,
            key=lambda c: (c.projection.ref_id, c.projection.query_chrom, c.projection.query_start),
        ):
            projection = item.projection
            writer.writerow(
                {
                    "projection_name": projection.projection_name,
                    "ref_cnee_id": projection.ref_id,
                    "initial_orthology_category": item.initial_orthology_category,
                    "query_chrom": projection.query_chrom,
                    "query_start": projection.query_start,
                    "query_end": projection.query_end,
                    "sequence_identity": (
                        "." if projection.sequence_identity is None else projection.sequence_identity
                    ),
                    "gap_fraction": projection.fraction_gap,
                }
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
