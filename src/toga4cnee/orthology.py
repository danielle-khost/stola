"""Graph-based CNEE orthology classification with synteny-verified pruning."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .assembly_gaps import DEFAULT_MIN_GAP_LENGTH, NRunIndex, scan_query_assembly_gaps
from .bed import parse_bed
from .chain import iter_chain_blocks
from .index import load_index
from .sequence_compare import SequenceComparator
from .summary import is_solid_match
from .synteny import BlockIndex, QueryBedRecord, QueryNeighborIndex, RefNeighborIndex, _flanking_blocks

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

REF_CNEE_COLUMNS = [
    "ref_cnee_id",
    "orthology_category",
    "num_projections_final",
    "num_projections_pruned",
    "pre_pruning_category",
]

# Verification criteria v2 (see synteny_verification.md): retain a projection
# only if it clears BOTH an align_score and a neighbor-ortholog-fraction bar.
DEFAULT_MIN_ALIGN_SCORE = 0.5
DEFAULT_MIN_NEIGHBOR_FRACTION = 0.5
DEFAULT_VERIFICATION_NUM_NEIGHBORS = 5

# Fallback ref category for a "lost" ref (chain overlap, zero aligned bases)
# when the query assembly's gap-marker convention can't be confidently
# inferred (or no query FASTA was supplied at all) -- see
# orthology_classification.md and assembly_gaps.py.
UNCERTAIN_LOSS_TYPE_CATEGORY = "missing/loss"

# Prefix prepended to every graph-assigned query CNEE id (e.g. "q_ce1234").
DEFAULT_QUERY_ID_PREFIX = "q_"


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
class LostGapCandidate:
    # A ref's best-ranked candidate chain among its zero-aligned-bases rows --
    # the chain whose gap gets examined for the probable_loss/assembly_gap_loss
    # distinction. Not yet known to be "lost" overall (this ref might still
    # have a *different* candidate chain with real aligned bases).
    ref_id: str
    ref_chrom: str
    ref_start: int
    ref_end: int
    chain_id: str


@dataclass(frozen=True)
class ProjectionEvidence:
    projections: list[Projection]
    refs_with_candidates: set[str]
    zero_aligned_best_row_by_ref: dict[str, LostGapCandidate]


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
    ref_cnee_rows: int


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
    query_fasta: str | Path | None = None,
    min_gap_length: int = DEFAULT_MIN_GAP_LENGTH,
    out_ref_cnees: str | Path | None = None,
    query_id_prefix: str = DEFAULT_QUERY_ID_PREFIX,
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs.

    Builds a ref<->query graph from every aligned projection, then (unless
    synteny_check is False) verifies and prunes projections in ambiguous
    categories -- a non-solid single ref/query edge (putative one2one(weak)),
    or any edge inside a one2many/many2one/complex component. "lost" refs
    (zero aligned query bases) are never part of this graph at all, so they
    are never subject to verification.

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

    Independently of verification, when query_fasta is given, every "lost"
    ref (zero aligned query bases anywhere) is further split into
    probable_loss or assembly_gap_loss by checking whether its chain gap's
    query-side span coincides with an inferred contig-boundary N-run in the
    query assembly (see assembly_gaps.py). When query_fasta is omitted, or
    the query assembly's gap-length convention can't be confidently
    inferred, or a given ref's own chain gap is itself ambiguous, this falls
    back to the collapsed "missing/loss" category instead of guessing. This
    sub-categorization, plus per-ref pruning counts, is written to
    out_ref_cnees if given (see output_files.md) -- out_orthology's own
    "lost"/"unaligned" labels are unaffected either way.

    Every graph-assigned query CNEE id is prefixed with query_id_prefix
    (default "q_", e.g. "q_ce1234").
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
    initial_category_by_ref: dict[str, str] = {}
    if synteny_check:
        surviving, culled, initial_category_by_ref = _verify_and_prune_projections(
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
    _assign_orthology_categories(components, final_clusters, query_id_prefix)
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

    ref_cnee_count = 0
    if out_ref_cnees is not None:
        lost_gap_candidates = {
            ref_id: candidate
            for ref_id, candidate in evidence.zero_aligned_best_row_by_ref.items()
            if ref_id not in refs_with_any_projection
        }
        n_run_index: NRunIndex | None = None
        if query_fasta is not None and lost_gap_candidates:
            gap_convention, candidate_index = scan_query_assembly_gaps(query_fasta, min_gap_length)
            if gap_convention.confidence != "none":
                n_run_index = candidate_index
        lost_category_by_ref = _classify_lost_refs(lost_gap_candidates, index_path, n_run_index)
        ref_cnee_count = _write_ref_cnee_tsv(
            ref_ids,
            final_clusters,
            initial_category_by_ref,
            culled,
            evidence.refs_with_candidates,
            refs_with_any_projection,
            lost_category_by_ref,
            out_ref_cnees,
        )

    return OrthologyResult(
        orthology_rows=orthology_count,
        query_bed_rows=bed_count,
        culled_projection_rows=len(culled),
        ref_cnee_rows=ref_cnee_count,
    )


def _read_projection_evidence(stats_path: str | Path) -> ProjectionEvidence:
    # Read candidate evidence and aligned query coordinates usable as graph edges.
    projections: list[Projection] = []
    refs_with_candidates: set[str] = set()
    zero_aligned_best_row_by_ref: dict[str, LostGapCandidate] = {}
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_uid",
            "ref_chrom",
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
            "rank_by_chain_score",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            refs_with_candidates.add(row["cnee_uid"])
            if int(row["aligned_bases"]) <= 0:
                # Track the best-ranked zero-aligned candidate chain per ref --
                # only meaningful for refs that turn out to have NO aligned
                # projection at all anywhere (filtered against
                # refs_with_any_projection by the caller).
                if int(row["rank_by_chain_score"]) == 1:
                    zero_aligned_best_row_by_ref[row["cnee_uid"]] = LostGapCandidate(
                        ref_id=row["cnee_uid"],
                        ref_chrom=row["ref_chrom"],
                        ref_start=int(row["ref_start"]),
                        ref_end=int(row["ref_end"]),
                        chain_id=row["chain_id"],
                    )
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
        zero_aligned_best_row_by_ref=zero_aligned_best_row_by_ref,
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
) -> tuple[list[Projection], list[CulledProjection], dict[str, str]]:
    initial_category_by_ref, trusted_ref_ids, component_id_by_ref = _initial_verification_state(
        nx, initial_graph, min_sequence_identity, max_fraction_gap
    )
    if not initial_category_by_ref:
        return list(projections), [], initial_category_by_ref

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
    return surviving, culled, initial_category_by_ref


# --------------------------------------------------------------------------
# "Lost" ref sub-categorization: a ref whose CNEE sits entirely inside a
# chain gap (zero aligned query bases) is split into probable_loss (the gap
# isn't a query assembly artifact) vs assembly_gap_loss (it is), using the
# query genome's inferred contig-boundary N-runs. Independent of synteny
# verification -- this always runs (when a query FASTA is available),
# regardless of synteny_check.
# --------------------------------------------------------------------------


def _classify_lost_refs(
    lost_gap_candidates: dict[str, LostGapCandidate],
    index_path: str | Path,
    n_run_index: NRunIndex | None,
) -> dict[str, str]:
    """Categorize each "lost" ref as probable_loss / assembly_gap_loss.

    Falls back to UNCERTAIN_LOSS_TYPE_CATEGORY ("missing/loss") whenever the
    distinction can't be made confidently: n_run_index is None (no query
    FASTA was supplied, or its gap-length convention couldn't be inferred
    confidently), or this particular ref's own chain gap is itself ambiguous
    (no aligned block on one side at all, i.e. it sits at the very edge of
    its chain's aligned span).
    """
    if not lost_gap_candidates:
        return {}

    category_by_ref: dict[str, str] = {}
    if n_run_index is None:
        return {ref_id: UNCERTAIN_LOSS_TYPE_CATEGORY for ref_id in lost_gap_candidates}

    payload, records = load_index(index_path)
    record_by_chain_id = {record.chain_id: record for record in records}
    chain_path = payload["indexed_chain_path"]
    block_cache: dict[str, BlockIndex] = {}

    with Path(chain_path).open("rb") as chain_handle:
        for ref_id, candidate in lost_gap_candidates.items():
            record = record_by_chain_id[candidate.chain_id]
            if candidate.chain_id not in block_cache:
                chain_handle.seek(record.start_byte)
                chain_text = chain_handle.read(record.byte_length).decode("utf-8")
                block_cache[candidate.chain_id] = BlockIndex.from_blocks(iter_chain_blocks(chain_text))
            blocks = block_cache[candidate.chain_id]

            preceding, following = _flanking_blocks(blocks, candidate.ref_start, candidate.ref_end)
            if preceding is None or following is None:
                category_by_ref[ref_id] = UNCERTAIN_LOSS_TYPE_CATEGORY
                continue

            if record.query_strand == "+":
                query_gap_start, query_gap_end = preceding.query_end, following.query_start
            else:
                query_gap_start, query_gap_end = following.query_end, preceding.query_start

            if query_gap_end > query_gap_start and n_run_index.overlaps(
                record.query_chrom, query_gap_start, query_gap_end
            ):
                category_by_ref[ref_id] = "assembly_gap_loss"
            else:
                category_by_ref[ref_id] = "probable_loss"
    return category_by_ref


# --------------------------------------------------------------------------
# Category assignment and output.
# --------------------------------------------------------------------------


def _assign_orthology_categories(
    components: list[set[str]],
    clusters: list[QueryCluster],
    query_id_prefix: str,
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
            query_ids = [f"{query_id_prefix}{_clean_id(ref_ids[0])}"]
        elif len(ref_nodes) == 1 and len(query_nodes) > 1:
            category = "one2many"
            query_ids = [
                f"{query_id_prefix}{_clean_id(ref_ids[0])}_{index}"
                for index in range(1, len(query_nodes) + 1)
            ]
        elif len(ref_nodes) > 1 and len(query_nodes) == 1:
            category = "many2one"
            query_ids = [f"{query_id_prefix}many2one_{component_counter:06d}"]
        else:
            category = "complex/many2many"
            query_ids = [
                f"{query_id_prefix}complex_{component_counter:06d}_{index}"
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
    # Write compact query/ref/category rows, including lost/unaligned/
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
                category = "unaligned"
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


def _write_ref_cnee_tsv(
    ref_ids: list[str],
    final_clusters: list[QueryCluster],
    initial_category_by_ref: dict[str, str],
    culled: list[CulledProjection],
    refs_with_candidates: set[str],
    refs_with_any_projection: set[str],
    lost_category_by_ref: dict[str, str],
    out_path: str | Path,
) -> int:
    """Write one row per reference CNEE (see output_files.md).

    orthology_category is the ref-based final category: one2one/one2many/
    many2one/complex/many2many for a ref with a surviving query edge,
    otherwise unaligned/probable_loss/assembly_gap_loss/"missing/loss"/
    uncertain_loss. pre_pruning_category is the same value for anything
    verification never touches (unaligned, the loss subtypes); for anything
    that went through verification it's the raw pre-pruning graph shape,
    which can differ from the final category (e.g. an initial
    complex/many2many that pruned down to a lone one2one edge).
    """
    final_category_by_ref: dict[str, str] = {}
    final_count_by_ref: dict[str, int] = defaultdict(int)
    for cluster in final_clusters:
        if not cluster.query_cnee_id:
            continue
        for ref_id in cluster.refs:
            final_category_by_ref[ref_id] = cluster.orthology_category
            final_count_by_ref[ref_id] += 1

    pruned_count_by_ref: dict[str, int] = defaultdict(int)
    for item in culled:
        pruned_count_by_ref[item.projection.ref_id] += 1

    count = 0
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=REF_CNEE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for ref_id in ref_ids:
            if ref_id in final_category_by_ref:
                category = final_category_by_ref[ref_id]
                pre_pruning_category = initial_category_by_ref.get(ref_id, category)
            elif ref_id not in refs_with_candidates:
                category = "unaligned"
                pre_pruning_category = "unaligned"
            elif ref_id not in refs_with_any_projection:
                category = lost_category_by_ref.get(ref_id, UNCERTAIN_LOSS_TYPE_CATEGORY)
                pre_pruning_category = category
            else:
                category = "uncertain_loss"
                pre_pruning_category = initial_category_by_ref.get(ref_id, "uncertain_loss")
            writer.writerow(
                {
                    "ref_cnee_id": ref_id,
                    "orthology_category": category,
                    "num_projections_final": final_count_by_ref.get(ref_id, 0),
                    "num_projections_pruned": pruned_count_by_ref.get(ref_id, 0),
                    "pre_pruning_category": pre_pruning_category,
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
