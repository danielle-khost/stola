"""Depth/score spurious-hit filtering for multi-projection ref CNEEs.

Replaces the old align_score/neighbor-synteny verification
(see synteny.py's RefNeighborIndex/QueryNeighborIndex, still used
independently by the exploratory synteny reports) for anything ambiguous --
a ref CNEE with more than one candidate query locus, or more than one ref
CNEE sharing a query locus.

Per candidate projection, at its own query-side locus (query_start..query_end,
already the aligned span from overlap_stats.tsv):

1. locus_depth = 1 + (# other chains, genome-wide, whose own query envelope
   overlaps this locus's core window extended by flank_window on each side).
   locus_depth >= depth_filter_max: reject outright (repeat-driven pileup).
2. locus_depth == 1 (no competing chain at all): accept iff the surrounding
   flanks (flank_window each side, core excluded) have aligned-block coverage
   from ANY chain (not just this one) above flank_any_coverage_threshold --
   a real ortholog's neighborhood is usually alignable even when this exact
   hit is uncontested; reject otherwise.
3. Otherwise (real competition, not a pileup): score this projection's own
   chain and every competitor's own chain over the core window (arithmetic
   match/mismatch/affine-gap score -- see score_window_for_chain), normalize
   by the locus's best score, and accept immediately if this projection is
   the (tied-)best and every competitor is well behind (R < cnee_r_threshold).
   Reject immediately if this projection itself is well behind
   (R < cnee_r_threshold). Otherwise the race is close: fall back to the
   same score computed over the two flanking windows (raw, not normalized;
   "missing" coverage counts as a gap run), accepting only if this
   projection's flank score beats every competitor's.

The "genome-wide chain search" and vectorized block parsing here are what
made a full-dataset validation run (140,579 candidate rows) complete in
~90s rather than the ~80 minutes a naive per-line chain.py-style parse would
take -- see parse_query_blocks_vectorized's docstring.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .models import ChainBlock, ChainRecord
from .sequence_compare import AlignedSegment

DEFAULT_DEPTH_FILTER_MAX = 20
DEFAULT_FLANK_WINDOW = 10_000
DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD = 0.5
DEFAULT_CNEE_R_THRESHOLD = 0.7
DEFAULT_MATCH_SCORE = 1.0
DEFAULT_MISMATCH_SCORE = -1.0
DEFAULT_OPEN_GAP_SCORE = -5.0
DEFAULT_EXTEND_GAP_SCORE = -1.0

# Outcome labels, shared vocabulary with the exploratory scripts this was
# validated against (exploratory_analysis/depth_and_score_filtering.py).
DISCARDED_REPEAT_DRIVEN_DEPTH = "discarded_repeat_driven_depth"
ACCEPTED_ISOLATED_WITH_FLANK_SUPPORT = "accepted_isolated_with_flank_support"
REJECTED_ISOLATED_NO_FLANK_SUPPORT = "rejected_isolated_no_flank_support"
ACCEPTED_CLEAR_CNEE_WINNER = "accepted_clear_cnee_winner"
REJECTED_CNEE_CLEAR_LOSER = "rejected_cnee_clear_loser"
ACCEPTED_FLANK_TIEBREAK_WINNER = "accepted_flank_tiebreak_winner"
REJECTED_FLANK_TIEBREAK_LOSER = "rejected_flank_tiebreak_loser"


@dataclass(frozen=True)
class QueryBlocks:
    """A chain's aligned blocks, sorted by query position."""

    query_starts: np.ndarray
    query_ends: np.ndarray
    ref_starts: np.ndarray
    ref_ends: np.ndarray

    def overlapping(self, start: int, end: int) -> list[ChainBlock]:
        if end <= start or self.query_starts.size == 0:
            return []
        lo = int(np.searchsorted(self.query_ends, start, side="right"))
        hi = int(np.searchsorted(self.query_starts, end, side="left"))
        if hi <= lo:
            return []
        return [
            ChainBlock(
                ref_start=int(self.ref_starts[i]), ref_end=int(self.ref_ends[i]),
                query_start=int(self.query_starts[i]), query_end=int(self.query_ends[i]),
            )
            for i in range(lo, hi)
        ]


def parse_query_blocks_vectorized(header_fields: list[str], body_lines: list[str]) -> QueryBlocks:
    """Vectorized replacement for chain.py's iter_chain_blocks.

    iter_chain_blocks walks the chain body one line at a time in pure
    Python, building a ChainBlock dataclass instance per line (with
    defensive min()/max() calls that are actually always no-ops, since a
    block's start is provably < its end by construction). Profiling a real
    genome-wide run showed that loop consuming ~80% of total wall time --
    almost entirely interpreter overhead, not the actual sequence-comparison
    work done elsewhere in this module.

    This computes the exact same (ref_start, ref_end, query_start, query_end)
    per block, but as whole-chain numpy arrays: every body line's (size, dt,
    dq) is tokenized in ONE big str.split() (the ragged last line, which
    omits dt/dq, is parsed separately) and cast to an int64 array in one
    shot, then both coordinate axes are obtained via np.cumsum instead of a
    running Python counter. Blocks are monotonic along each axis (a single
    chain's blocks never reorder or overlap), so getting them into ascending
    query order needs at most one np.flip -- never an O(n log n) sort.
    """
    ref_strand = header_fields[4]
    ref_size = int(header_fields[3])
    ref_raw_start = int(header_fields[5])
    query_strand = header_fields[9]
    query_size = int(header_fields[8])
    query_raw_start = int(header_fields[10])

    n = len(body_lines)
    if n == 0:
        raise ValueError("Cannot parse a chain with no body lines")

    if n > 1:
        tokens = " ".join(body_lines[:-1]).split()
        regular = np.array(tokens, dtype=np.int64).reshape(-1, 3)
        sizes, dts, dqs = regular[:, 0], regular[:, 1], regular[:, 2]
    else:
        sizes = np.empty(0, dtype=np.int64)
        dts = np.empty(0, dtype=np.int64)
        dqs = np.empty(0, dtype=np.int64)

    last_size = int(body_lines[-1].split()[0])
    sizes = np.concatenate([sizes, [last_size]])
    dts = np.concatenate([dts, [0]])
    dqs = np.concatenate([dqs, [0]])

    prefix_ref = np.concatenate(([0], np.cumsum(sizes + dts)[:-1]))
    prefix_query = np.concatenate(([0], np.cumsum(sizes + dqs)[:-1]))

    if ref_strand == "+":
        ref_start = ref_raw_start + prefix_ref
        ref_end = ref_start + sizes
    elif ref_strand == "-":
        ref_end = (ref_size - ref_raw_start) - prefix_ref
        ref_start = ref_end - sizes
    else:
        raise ValueError(f"Invalid chain target strand: {ref_strand}")

    if query_strand == "+":
        query_start = query_raw_start + prefix_query
        query_end = query_start + sizes
    elif query_strand == "-":
        query_end = (query_size - query_raw_start) - prefix_query
        query_start = query_end - sizes
    else:
        raise ValueError(f"Invalid chain query strand: {query_strand}")

    if query_strand == "-":
        ref_start, ref_end, query_start, query_end = (
            arr[::-1] for arr in (ref_start, ref_end, query_start, query_end)
        )
    return QueryBlocks(query_start, query_end, ref_start, ref_end)


class ChainEnvelopeCache:
    """Lazily parse+cache each chain's aligned blocks, by chain_id."""

    def __init__(self, indexed_chain_path: str | Path, records_by_id: dict[str, ChainRecord]) -> None:
        self.handle = Path(indexed_chain_path).open("rb")
        self.records_by_id = records_by_id
        self.blocks_by_id: dict[str, QueryBlocks] = {}

    def get(self, chain_id: str) -> tuple[ChainRecord, QueryBlocks]:
        record = self.records_by_id[chain_id]
        if chain_id not in self.blocks_by_id:
            self.handle.seek(record.start_byte)
            chain_text = self.handle.read(record.byte_length).decode("utf-8")
            lines = [line for line in chain_text.splitlines() if line]
            header_fields = lines[0].split()
            self.blocks_by_id[chain_id] = parse_query_blocks_vectorized(header_fields, lines[1:])
        return record, self.blocks_by_id[chain_id]

    def close(self) -> None:
        self.handle.close()


def build_query_envelopes(records: list[ChainRecord]) -> dict[str, dict[str, np.ndarray]]:
    """Per-query-chromosome chain envelope arrays, for cheap overlap search."""
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for record in records:
        by_chrom.setdefault(record.query_chrom, []).append(
            (record.query_start, record.query_end, record.chain_id)
        )
    envelopes: dict[str, dict[str, np.ndarray]] = {}
    for chrom, rows in by_chrom.items():
        starts, ends, chain_ids = zip(*rows)
        envelopes[chrom] = {
            "starts": np.array(starts, dtype=np.int64),
            "ends": np.array(ends, dtype=np.int64),
            "chain_ids": np.array(chain_ids),
        }
    return envelopes


def find_competitor_ids(
    query_chrom: str, own_chain_id: str, window_start: int, window_end: int,
    envelopes: dict[str, dict[str, np.ndarray]],
) -> list[str]:
    envelope = envelopes.get(query_chrom)
    if envelope is None:
        return []
    mask = (
        (envelope["starts"] < window_end)
        & (envelope["ends"] > window_start)
        & (envelope["chain_ids"] != own_chain_id)
    )
    return [str(cid) for cid in envelope["chain_ids"][mask]]


def aligned_segments_for_query_interval(
    blocks: list[ChainBlock], query_start: int, query_end: int, query_strand: str
) -> list[AlignedSegment]:
    """Query-anchored mirror of sequence_compare.py's aligned_segments_for_interval.

    Given blocks known to overlap [query_start, query_end) on the query
    side, return the paired ref/query sub-intervals clipped to that window.
    Correct because a chain block is an ungapped 1:1 run: moving across it
    in the query-increasing direction corresponds to a fixed
    (query_strand-dependent) direction in the reference.
    """
    segments: list[AlignedSegment] = []
    for block in blocks:
        overlap_start = max(query_start, block.query_start)
        overlap_end = min(query_end, block.query_end)
        if overlap_end <= overlap_start:
            continue
        offset_start = overlap_start - block.query_start
        offset_end = overlap_end - block.query_start
        if query_strand == "+":
            ref_start = block.ref_start + offset_start
            ref_end = block.ref_start + offset_end
        else:
            ref_start = block.ref_end - offset_end
            ref_end = block.ref_end - offset_start
        segments.append(
            AlignedSegment(
                ref_start=min(ref_start, ref_end), ref_end=max(ref_start, ref_end),
                query_start=overlap_start, query_end=overlap_end,
            )
        )
    return segments


def gap_run_lengths(covered_spans: list[tuple[int, int]], window_start: int, window_end: int) -> list[int]:
    """Lengths of maximal uncovered sub-intervals of [window_start, window_end)."""
    ordered = sorted(covered_spans)
    runs: list[int] = []
    cursor = window_start
    for start, end in ordered:
        if start > cursor:
            runs.append(start - cursor)
        cursor = max(cursor, end)
    if window_end > cursor:
        runs.append(window_end - cursor)
    return runs


def affine_gap_penalty(run_lengths: list[int], open_gap_score: float, extend_gap_score: float) -> float:
    return sum(
        open_gap_score + (length - 1) * extend_gap_score for length in run_lengths if length > 0
    )


def score_window_for_chain(
    ref_accessor, query_accessor, ref_chrom: str, query_chrom: str, query_strand: str,
    window_start: int, window_end: int, blocks: list[ChainBlock],
    match_score: float, mismatch_score: float, open_gap_score: float, extend_gap_score: float,
) -> float:
    """Arithmetic match/mismatch/affine-gap score for one chain over one query-anchored window.

    Every window base not covered by one of this chain's own blocks --
    whether an internal chain gap, or the chain simply not reaching that far
    -- becomes part of a gap run. Floored at 0 (a heavily gapped/mismatched
    chain reads as "no real signal here" rather than an arbitrarily negative
    number), mirroring the local-alignment convention of resetting a
    score that would go negative.
    """
    if window_end <= window_start:
        return 0.0

    segments = aligned_segments_for_query_interval(blocks, window_start, window_end, query_strand)
    if not segments:
        return max(0.0, affine_gap_penalty([window_end - window_start], open_gap_score, extend_gap_score))

    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required for sequence comparison. Activate the biopython conda environment."
        ) from exc

    total_matches = 0
    total_mismatches = 0
    covered_spans: list[tuple[int, int]] = []
    for segment in segments:
        ref_part = ref_accessor.fetch(ref_chrom, segment.ref_start, segment.ref_end)
        query_part = query_accessor.fetch(query_chrom, segment.query_start, segment.query_end)
        if query_strand == "-":
            query_part = str(Seq(query_part).reverse_complement()).upper()
        compared = min(len(ref_part), len(query_part))
        if compared == 0:
            continue
        ref_arr = np.frombuffer(ref_part[:compared].encode("ascii"), dtype=np.uint8)
        query_arr = np.frombuffer(query_part[:compared].encode("ascii"), dtype=np.uint8)
        matches = int(np.count_nonzero(ref_arr == query_arr))
        total_matches += matches
        total_mismatches += compared - matches
        covered_spans.append((segment.query_start, segment.query_start + compared))

    gap_runs = gap_run_lengths(covered_spans, window_start, window_end)
    raw_score = (
        total_matches * match_score
        + total_mismatches * mismatch_score
        + affine_gap_penalty(gap_runs, open_gap_score, extend_gap_score)
    )
    return max(0.0, raw_score)


def any_chain_flank_coverage(
    query_chrom: str, own_chain_id: str,
    upstream_start: int, upstream_end: int, downstream_start: int, downstream_end: int,
    envelopes: dict[str, dict[str, np.ndarray]], cache: ChainEnvelopeCache,
) -> float:
    """Fraction of the two flank sub-windows covered by ANY chain's aligned blocks."""
    total_window = max(0, upstream_end - upstream_start) + max(0, downstream_end - downstream_start)
    if total_window <= 0:
        return float("nan")

    envelope = envelopes.get(query_chrom)
    covered = 0
    for start, end in ((upstream_start, upstream_end), (downstream_start, downstream_end)):
        if end <= start:
            continue
        chain_ids = {own_chain_id}
        if envelope is not None:
            mask = (envelope["starts"] < end) & (envelope["ends"] > start)
            chain_ids.update(str(cid) for cid in envelope["chain_ids"][mask])
        track = np.zeros(end - start, dtype=bool)
        for chain_id in chain_ids:
            _, blocks = cache.get(chain_id)
            for block in blocks.overlapping(start, end):
                b_start = max(start, block.query_start) - start
                b_end = min(end, block.query_end) - start
                if b_end > b_start:
                    track[b_start:b_end] = True
        covered += int(track.sum())
    return covered / total_window


@dataclass(frozen=True)
class ProjectionVerdict:
    accepted: bool
    branch: str  # "depth_filtered" | "isolated" | "competing"
    outcome: str
    locus_depth: int
    n_competing_chains: int
    own_cnee_score: float | None = None
    best_competitor_cnee_score: float | None = None
    r_own_cnee: float | None = None
    r_best_competitor_cnee: float | None = None
    own_flank_score: float | None = None
    best_competitor_flank_score: float | None = None
    flank_any_coverage: float | None = None


def classify_projection(
    chain_id: str, query_chrom: str, query_start: int, query_end: int, ref_chrom: str,
    records_by_id: dict[str, ChainRecord], envelopes: dict[str, dict[str, np.ndarray]],
    cache: ChainEnvelopeCache, ref_accessor, query_accessor,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    depth_filter_max: int = DEFAULT_DEPTH_FILTER_MAX,
    flank_any_coverage_threshold: float = DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    cnee_r_threshold: float = DEFAULT_CNEE_R_THRESHOLD,
    match_score: float = DEFAULT_MATCH_SCORE,
    mismatch_score: float = DEFAULT_MISMATCH_SCORE,
    open_gap_score: float = DEFAULT_OPEN_GAP_SCORE,
    extend_gap_score: float = DEFAULT_EXTEND_GAP_SCORE,
) -> ProjectionVerdict:
    """Accept/reject one candidate projection using its own query-side locus."""
    record = records_by_id[chain_id]
    locus_window_start = max(0, query_start - flank_window)
    locus_window_end = min(record.query_size, query_end + flank_window)

    competitor_ids = find_competitor_ids(query_chrom, chain_id, locus_window_start, locus_window_end, envelopes)
    n_competing_chains = len(competitor_ids)
    locus_depth = n_competing_chains + 1

    if depth_filter_max and locus_depth >= depth_filter_max:
        return ProjectionVerdict(
            accepted=False, branch="depth_filtered", outcome=DISCARDED_REPEAT_DRIVEN_DEPTH,
            locus_depth=locus_depth, n_competing_chains=n_competing_chains,
        )

    if n_competing_chains == 0:
        upstream_start, upstream_end = max(0, query_start - flank_window), query_start
        downstream_start, downstream_end = query_end, min(record.query_size, query_end + flank_window)
        coverage = any_chain_flank_coverage(
            query_chrom, chain_id, upstream_start, upstream_end, downstream_start, downstream_end,
            envelopes, cache,
        )
        accepted = coverage > flank_any_coverage_threshold
        return ProjectionVerdict(
            accepted=accepted, branch="isolated",
            outcome=ACCEPTED_ISOLATED_WITH_FLANK_SUPPORT if accepted else REJECTED_ISOLATED_NO_FLANK_SUPPORT,
            locus_depth=locus_depth, n_competing_chains=n_competing_chains, flank_any_coverage=coverage,
        )

    _, own_blocks = cache.get(chain_id)
    own_cnee_score = score_window_for_chain(
        ref_accessor, query_accessor, ref_chrom, query_chrom, record.query_strand,
        query_start, query_end, own_blocks.overlapping(query_start, query_end),
        match_score, mismatch_score, open_gap_score, extend_gap_score,
    )
    competitor_cnee_scores: dict[str, float] = {}
    for competitor_id in competitor_ids:
        comp_record, comp_blocks = cache.get(competitor_id)
        competitor_cnee_scores[competitor_id] = score_window_for_chain(
            ref_accessor, query_accessor, comp_record.ref_chrom, query_chrom, comp_record.query_strand,
            query_start, query_end, comp_blocks.overlapping(query_start, query_end),
            match_score, mismatch_score, open_gap_score, extend_gap_score,
        )
    best_competitor_cnee_score = max(competitor_cnee_scores.values())

    s_max = max(own_cnee_score, best_competitor_cnee_score)
    if s_max > 0:
        r_own = own_cnee_score / s_max
        r_best_competitor = best_competitor_cnee_score / s_max
        if r_own >= 1.0 and r_best_competitor < cnee_r_threshold:
            return ProjectionVerdict(
                accepted=True, branch="competing", outcome=ACCEPTED_CLEAR_CNEE_WINNER,
                locus_depth=locus_depth, n_competing_chains=n_competing_chains,
                own_cnee_score=own_cnee_score, best_competitor_cnee_score=best_competitor_cnee_score,
                r_own_cnee=r_own, r_best_competitor_cnee=r_best_competitor,
            )
        if r_own < cnee_r_threshold:
            return ProjectionVerdict(
                accepted=False, branch="competing", outcome=REJECTED_CNEE_CLEAR_LOSER,
                locus_depth=locus_depth, n_competing_chains=n_competing_chains,
                own_cnee_score=own_cnee_score, best_competitor_cnee_score=best_competitor_cnee_score,
                r_own_cnee=r_own, r_best_competitor_cnee=r_best_competitor,
            )
        # else: too close to call from the CNEE score alone -- fall through
        # to the flank tiebreak below. (r_own/r_best_competitor still
        # reported in the returned verdict for auditability.)
    else:
        r_own = r_best_competitor = None

    upstream_start, upstream_end = max(0, query_start - flank_window), query_start
    downstream_start, downstream_end = query_end, min(record.query_size, query_end + flank_window)

    def flank_score_for(ref_chrom_: str, query_strand: str, blocks: QueryBlocks) -> float:
        return (
            score_window_for_chain(
                ref_accessor, query_accessor, ref_chrom_, query_chrom, query_strand,
                upstream_start, upstream_end, blocks.overlapping(upstream_start, upstream_end),
                match_score, mismatch_score, open_gap_score, extend_gap_score,
            )
            + score_window_for_chain(
                ref_accessor, query_accessor, ref_chrom_, query_chrom, query_strand,
                downstream_start, downstream_end, blocks.overlapping(downstream_start, downstream_end),
                match_score, mismatch_score, open_gap_score, extend_gap_score,
            )
        )

    own_flank_score = flank_score_for(ref_chrom, record.query_strand, own_blocks)
    best_competitor_flank_score = max(
        flank_score_for(records_by_id[cid].ref_chrom, records_by_id[cid].query_strand, cache.get(cid)[1])
        for cid in competitor_ids
    )
    accepted = own_flank_score > best_competitor_flank_score
    return ProjectionVerdict(
        accepted=accepted, branch="competing",
        outcome=ACCEPTED_FLANK_TIEBREAK_WINNER if accepted else REJECTED_FLANK_TIEBREAK_LOSER,
        locus_depth=locus_depth, n_competing_chains=n_competing_chains,
        own_cnee_score=own_cnee_score, best_competitor_cnee_score=best_competitor_cnee_score,
        r_own_cnee=r_own, r_best_competitor_cnee=r_best_competitor,
        own_flank_score=own_flank_score, best_competitor_flank_score=best_competitor_flank_score,
    )
