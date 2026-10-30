"""Graph-based CNEE orthology classification with relative-score filtering."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

import networkx as nx

from .assembly_gaps import DEFAULT_MIN_GAP_LENGTH, NRunIndex, scan_query_assembly_gaps
from .bed import parse_bed
from .chain import BlockIndex, _flanking_blocks, iter_chain_blocks, read_chain_text
from .index import ChainIndex, as_chain_index
from .models import ChainRecord
from .sequence_compare import SequenceComparator
from .spurious_filter import (
    ChainEnvelopeCache,
    DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    DEFAULT_FLANK_WINDOW,
    DEFAULT_RELATIVE_SCORE_THRESHOLD,
    DEFAULT_SOFT_DEPTH_LIMIT,
    ProjectionVerdict,
    build_query_envelopes,
    decide_locus,
    score_locus,
)

ORTHOLOGY_COLUMNS = [
    "query_cnee_id",
    "partner_ref_cnees",
    "orthology_category",
    "sequence_identity",
    "gap_fraction",
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
    "rejection_reason",
    "locus_depth",
    "local_score",
    "best_competitor_local_score",
    "relative_score",
    "flank_any_coverage",
]

REF_CNEE_COLUMNS = [
    "ref_cnee_id",
    "orthology_category",
    "num_projections_final",
    "num_projections_pruned",
    "pre_pruning_category",
]

# Fallback ref category for a "lost" ref (chain overlap, zero aligned bases)
# when the query assembly's gap-marker convention can't be confidently
# inferred (or no query FASTA was supplied at all) -- see assembly_gaps.py.
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
    verdict: ProjectionVerdict | None = None


@dataclass
class EdgeMetrics:
    projection_metrics: list[tuple[float | None, float]] = field(default_factory=list)

    def add_projection(self, projection: Projection) -> None:
        # Keep each projection's own (identity, gap) pair rather than
        # blending projections that land on the same ref/query edge.
        self.projection_metrics.append(
            (projection.sequence_identity, projection.fraction_gap)
        )

    @property
    def projection_count(self) -> int:
        return len(self.projection_metrics)

    def representative(self) -> tuple[float | None, float]:
        # The (identity, gap) pair to report for this edge in output: the
        # highest-identity projection among the (possibly several) surviving
        # projections that landed on it.
        return max(
            self.projection_metrics,
            key=lambda pair: pair[0] if pair[0] is not None else -1.0,
        )


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
    chain_index: ChainIndex | str | Path,
    out_orthology: str | Path,
    out_query_bed: str | Path,
    sequence_comparator: SequenceComparator | None = None,
    filter_multi_projections: bool = True,
    soft_depth_limit: int = DEFAULT_SOFT_DEPTH_LIMIT,
    relative_score_threshold: float = DEFAULT_RELATIVE_SCORE_THRESHOLD,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    flank_any_coverage_threshold: float = DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    out_culled_projections: str | Path | None = None,
    query_fasta: str | Path | None = None,
    min_gap_length: int = DEFAULT_MIN_GAP_LENGTH,
    out_ref_cnees: str | Path | None = None,
    query_id_prefix: str = DEFAULT_QUERY_ID_PREFIX,
    progress: Callable[[str], None] | None = None,
) -> OrthologyResult:
    """Classify query CNEE loci from overlap stats and write TSV/BED outputs.

    Builds a ref<->query graph from every aligned projection, then (unless
    filter_multi_projections is False) filters projections belonging to any
    ref CNEE whose initial graph shape is ambiguous -- one2many, many2one,
    or complex/many2many (a lone ref/query edge is never touched: it is
    always reported as "one2one", with its own sequence_identity/
    gap_fraction exposed as output columns instead of being gated on them).
    "lost" refs (zero aligned query bases) are never part of this graph at
    all, so they are never subject to filtering either.

    Filtering (see spurious_filter.py) works per query locus -- the merged
    query span of the eligible projections landing there -- scoring every
    chain aligned at that locus once:
      - soft depth filter: of the chains whose query span overlaps the
        locus, keep the soft_depth_limit best by chain score that actually
        have aligned bases over it; a projection on any other chain is
        rejected (repeat-driven pileup).
      - exactly one kept chain: accept iff the projection's flanks
        (flank_window each side) have aligned-block coverage from ANY chain
        above flank_any_coverage_threshold.
      - otherwise: local score = sequence identity over the locus x aligned
        fraction over the locus plus flank_window each side; accept iff
        this projection's chain scores above relative_score_threshold
        relative to the locus's best chain. Several chains can pass at one
        locus, so many2one/one2many calls can survive filtering.

    filter_multi_projections=True requires sequence_comparator (real ref/
    query genome sequence access is needed to score candidates).

    Query loci are re-clustered from the surviving projections before final
    categories are assigned, so filtering can split a cluster that was only
    held together by a since-removed projection. A ref left with no
    surviving query edge (but that did have some projection) is written to
    out_orthology as "uncertain_loss"; if out_culled_projections is given,
    every individual rejected projection is also logged there for tracking,
    along with why it was rejected.

    Independently of filtering, when query_fasta is given, every "lost" ref
    (zero aligned query bases anywhere) is further split into probable_loss
    or assembly_gap_loss by checking whether its chain gap's query-side span
    coincides with an inferred contig-boundary N-run in the query assembly
    (see assembly_gaps.py). When query_fasta is omitted, or the query
    assembly's gap-length convention can't be confidently inferred, or a
    given ref's own chain gap is itself ambiguous, this falls back to the
    collapsed "missing/loss" category instead of guessing. This sub-
    categorization, plus per-ref pruning counts, is written to out_ref_cnees
    if given (see the README's Output section) -- out_orthology's own
    "lost"/"unaligned" labels are unaffected either way.

    Every graph-assigned query CNEE id is prefixed with query_id_prefix
    (default "q_", e.g. "q_ce1234"). A trailing "_" is appended automatically
    if query_id_prefix doesn't already end with one, so "aphWoo" and
    "aphWoo_" both produce ids like "aphWoo_ce1234".

    progress, if given, is called with a one-line message at each sub-step,
    and every 10% of the way through filtering.
    """
    report = progress or (lambda message: None)
    if filter_multi_projections and soft_depth_limit < 1:
        raise ValueError("soft_depth_limit must be at least 1")
    if query_id_prefix and not query_id_prefix.endswith("_"):
        query_id_prefix = f"{query_id_prefix}_"

    ref_intervals = parse_bed(cnee_bed)
    ref_ids = [interval.cnee_uid for interval in ref_intervals]
    evidence = _read_projection_evidence(stats_path)
    refs_with_any_projection = {projection.ref_id for projection in evidence.projections}
    report(
        f"Read {len(evidence.projections)} aligned projections of "
        f"{len(refs_with_any_projection)} reference CNEEs"
    )

    initial_clusters = _cluster_query_projections(evidence.projections)
    initial_graph = _build_orthology_graph(ref_ids, initial_clusters)

    # Chain index is needed both by filtering (chain-envelope search + block
    # parsing for scoring) and by the lost-ref gap-boundary check below;
    # loaded at most once (if given as a path) and shared between them.
    need_chain_index = filter_multi_projections or out_ref_cnees is not None
    chain_records: list[ChainRecord] = []
    indexed_chain_path: Path | None = None
    if need_chain_index:
        chain_index = as_chain_index(chain_index)
        chain_records = chain_index.records
        indexed_chain_path = chain_index.indexed_chain_path

    culled: list[CulledProjection] = []
    initial_category_by_ref: dict[str, str] = {}
    if filter_multi_projections:
        initial_category_by_ref = _initial_category_by_ref(initial_graph)
        if initial_category_by_ref:
            if sequence_comparator is None:
                raise ValueError(
                    "sequence_comparator is required when filter_multi_projections is True "
                    "(real ref/query sequence access is needed to score candidate projections "
                    "against their competitors)"
                )
            surviving, culled = _filter_multi_projections(
                evidence.projections,
                initial_category_by_ref,
                chain_records,
                indexed_chain_path,
                sequence_comparator,
                soft_depth_limit,
                relative_score_threshold,
                flank_window,
                flank_any_coverage_threshold,
                report,
            )
        else:
            surviving = list(evidence.projections)
        report("Assigning final orthology categories ...")
        final_clusters = _cluster_query_projections(surviving)
        final_graph = _build_orthology_graph(ref_ids, final_clusters)
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
        if lost_gap_candidates:
            report(
                f"Classifying {len(lost_gap_candidates)} lost reference CNEEs as probable "
                "or assembly-gap losses ..."
            )
        n_run_index: NRunIndex | None = None
        if query_fasta is not None and lost_gap_candidates:
            gap_convention, candidate_index = scan_query_assembly_gaps(query_fasta, min_gap_length)
            if gap_convention.confidence != "none":
                n_run_index = candidate_index
        lost_category_by_ref = _classify_lost_refs(
            lost_gap_candidates, chain_records, indexed_chain_path, n_run_index
        )
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


def _iter_overlapping_projection_groups(projections: list[Projection]) -> Iterator[list[Projection]]:
    # Merge overlapping query-side projections into query loci. Grouping by
    # (chrom, strand) first ensures overlapping projections on opposite
    # strands are never merged into the same query locus.
    by_group: dict[tuple[str, str], list[Projection]] = defaultdict(list)
    for projection in projections:
        by_group[(projection.query_chrom, projection.query_strand)].append(projection)

    for group_key in sorted(by_group):
        current: list[Projection] = []
        current_end = 0
        for projection in sorted(
            by_group[group_key],
            key=lambda item: (item.query_start, item.query_end, item.ref_id),
        ):
            if current and projection.query_start >= current_end:
                yield current
                current = []
            if not current:
                current_end = projection.query_end
            current.append(projection)
            current_end = max(current_end, projection.query_end)
        if current:
            yield current


def _cluster_query_projections(projections: list[Projection]) -> list[QueryCluster]:
    # Provisional query CNEE loci, one per group of overlapping projections.
    clusters: list[QueryCluster] = []
    for group in _iter_overlapping_projection_groups(projections):
        cluster = QueryCluster(
            node_id=f"query:{len(clusters) + 1:06d}",
            chrom=group[0].query_chrom,
            start=group[0].query_start,
            end=max(projection.query_end for projection in group),
        )
        for projection in group:
            cluster.strands.add(projection.query_strand)
            cluster.refs.add(projection.ref_id)
            _add_cluster_edge_metrics(cluster, projection)
            cluster.projection_count += 1
        clusters.append(cluster)
    return clusters


def _add_cluster_edge_metrics(cluster: QueryCluster, projection: Projection) -> None:
    # Record this projection's own metrics on the ref/query edge.
    metrics = cluster.edge_metrics.setdefault(projection.ref_id, EdgeMetrics())
    metrics.add_projection(projection)


def _build_orthology_graph(ref_ids: list[str], clusters: list[QueryCluster]):
    # Build a ref/query graph, recording each edge's projection metrics.
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


# --------------------------------------------------------------------------
# Relative-score filtering: applied only to projections belonging to a ref
# CNEE whose initial graph shape is ambiguous (one2many/many2one/complex/
# many2many). A lone ref/query edge is never filtered -- it is always
# reported as "one2one" (see spurious_filter.py for the per-locus decision
# itself).
# --------------------------------------------------------------------------


def _initial_category_by_ref(initial_graph) -> dict[str, str]:
    """Every ref's pre-filtering graph-shape category.

    A ref absent from this dict is either a lone one2one edge (never
    filtered, always kept) or has no projections at all.
    """
    initial_category_by_ref: dict[str, str] = {}
    for component in nx.connected_components(initial_graph):
        ref_nodes = [node for node in component if str(node).startswith("ref:")]
        query_nodes = [node for node in component if str(node).startswith("query:")]
        if not ref_nodes or not query_nodes:
            continue
        if len(ref_nodes) == 1 and len(query_nodes) == 1:
            continue
        if len(ref_nodes) == 1:
            category = "one2many"
        elif len(query_nodes) == 1:
            category = "many2one"
        else:
            category = "complex/many2many"
        for node in ref_nodes:
            initial_category_by_ref[_node_ref_id(node)] = category
    return initial_category_by_ref


def _filter_multi_projections(
    projections: list[Projection],
    initial_category_by_ref: dict[str, str],
    chain_records: list[ChainRecord],
    indexed_chain_path: str | Path,
    sequence_comparator: SequenceComparator,
    soft_depth_limit: int,
    relative_score_threshold: float,
    flank_window: int,
    flank_any_coverage_threshold: float,
    report: Callable[[str], None],
) -> tuple[list[Projection], list[CulledProjection]]:
    records_by_id = {record.chain_id: record for record in chain_records}
    envelopes = build_query_envelopes(chain_records)
    cache = ChainEnvelopeCache(indexed_chain_path, records_by_id)
    try:
        surviving: list[Projection] = []
        eligible: list[Projection] = []
        for projection in projections:
            if projection.ref_id in initial_category_by_ref:
                eligible.append(projection)
            else:
                surviving.append(projection)

        culled: list[CulledProjection] = []
        loci = list(_iter_overlapping_projection_groups(eligible))
        report(
            f"Filtering {len(eligible)} projections of {len(initial_category_by_ref)} "
            f"multi-projection reference CNEEs, at {len(loci)} query loci ..."
        )
        # Report at every 10% of loci scored. With under 10 loci, milestones
        # landing on the same locus collapse into one message (the highest).
        milestones = {(len(loci) * tenth + 9) // 10: tenth * 10 for tenth in range(1, 11)}
        # Every eligible projection sharing a query locus is judged against
        # the same set of chains, scored once at that locus.
        for locus_number, locus_projections in enumerate(loci, start=1):
            scores = score_locus(
                locus_projections,
                records_by_id,
                envelopes,
                cache,
                sequence_comparator.ref,
                sequence_comparator.query,
                flank_window,
                soft_depth_limit,
            )
            verdicts = decide_locus(
                locus_projections, scores, relative_score_threshold, flank_any_coverage_threshold
            )
            for projection, verdict in zip(locus_projections, verdicts):
                if verdict.accepted:
                    surviving.append(projection)
                else:
                    culled.append(
                        CulledProjection(
                            projection=projection,
                            initial_orthology_category=initial_category_by_ref[projection.ref_id],
                            verdict=verdict,
                        )
                    )
            if locus_number in milestones:
                report(
                    f"Scored {locus_number}/{len(loci)} query loci ({milestones[locus_number]}%)"
                )
        return surviving, culled
    finally:
        cache.close()


# --------------------------------------------------------------------------
# "Lost" ref sub-categorization: a ref whose CNEE sits entirely inside a
# chain gap (zero aligned query bases) is split into probable_loss (the gap
# isn't a query assembly artifact) vs assembly_gap_loss (it is), using the
# query genome's inferred contig-boundary N-runs. Independent of filtering --
# this always runs (when a query FASTA is available), regardless of
# filter_multi_projections.
# --------------------------------------------------------------------------


def _classify_lost_refs(
    lost_gap_candidates: dict[str, LostGapCandidate],
    chain_records: list[ChainRecord],
    indexed_chain_path: str | Path,
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

    record_by_chain_id = {record.chain_id: record for record in chain_records}
    block_cache: dict[str, BlockIndex] = {}

    with Path(indexed_chain_path).open("rb") as chain_handle:
        for ref_id, candidate in lost_gap_candidates.items():
            record = record_by_chain_id[candidate.chain_id]
            if candidate.chain_id not in block_cache:
                block_cache[candidate.chain_id] = BlockIndex.from_blocks(
                    iter_chain_blocks(read_chain_text(chain_handle, record))
                )
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
    # A surviving lone ref/query edge is always "one2one" -- there is
    # nothing left to distinguish it from an edge that always had a single
    # candidate, since no accept/reject decision is ever applied to it.
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


def query_id_prefix_from_one2one(query_cnee_id: str, ref_cnee_id: str) -> str | None:
    """Recover a run's query id prefix (with its trailing "_") from a one2one
    pair, whose query id is the prefix followed by the cleaned ref id (see
    _assign_orthology_categories). None if the id doesn't have that shape."""
    suffix = _clean_id(ref_cnee_id)
    if not suffix or not query_cnee_id.endswith(suffix):
        return None
    return query_cnee_id[: -len(suffix)]


def _write_orthology_tsv(
    clusters: list[QueryCluster],
    ref_ids: list[str],
    refs_with_candidates: set[str],
    refs_with_any_projection: set[str],
    out_path: str | Path,
) -> int:
    # Write compact query/ref/category rows, including lost/unaligned/
    # uncertain_loss references. Since clusters here are already the final
    # (post-filtering) result, a cluster's own .refs set IS the correct
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
            identities: list[str] = []
            gaps: list[str] = []
            for ref_id in partner_refs:
                identity, gap = cluster.edge_metrics[ref_id].representative()
                identities.append("." if identity is None else str(identity))
                gaps.append(str(gap))
            writer.writerow(
                {
                    "query_cnee_id": cluster.query_cnee_id,
                    "partner_ref_cnees": ",".join(partner_refs),
                    "orthology_category": cluster.orthology_category,
                    "sequence_identity": ",".join(identities),
                    "gap_fraction": ",".join(gaps),
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
                    "sequence_identity": ".",
                    "gap_fraction": ".",
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
    # Log every individual projection filtering rejected, for tracking --
    # these refs end up "uncertain_loss" in out_orthology (unless a sibling
    # projection of the same ref survived), but this is the only place the
    # specific rejected projection(s) themselves, and why, are recorded.
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
            verdict = item.verdict
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
                    "rejection_reason": verdict.outcome if verdict else ".",
                    "locus_depth": verdict.locus_depth if verdict else ".",
                    "local_score": _optional_verdict_value(verdict, "local_score"),
                    "best_competitor_local_score": _optional_verdict_value(
                        verdict, "best_competitor_local_score"
                    ),
                    "relative_score": _optional_verdict_value(verdict, "relative_score"),
                    "flank_any_coverage": _optional_verdict_value(verdict, "flank_any_coverage"),
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
    """Write one row per reference CNEE (see the README's Output section).

    orthology_category is the ref-based final category: one2one/one2many/
    many2one/complex/many2many for a ref with a surviving query edge,
    otherwise unaligned/probable_loss/assembly_gap_loss/"missing/loss"/
    uncertain_loss. pre_pruning_category is the same value for anything
    filtering never touches (unaligned, the loss subtypes, and lone one2one
    edges); for anything that went through filtering it's the raw
    pre-filtering graph shape, which can differ from the final category
    (e.g. an initial complex/many2many that filtered down to a lone one2one
    edge).
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


def _optional_verdict_value(verdict: ProjectionVerdict | None, field_name: str) -> float | str:
    # A verdict metric for TSV output, or "." when that branch never computed it.
    value = None if verdict is None else getattr(verdict, field_name)
    return "." if value is None else value
