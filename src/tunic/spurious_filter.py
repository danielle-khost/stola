"""Relative-score spurious-hit filtering for multi-projection ref CNEEs.

Applied to anything ambiguous -- a ref CNEE with more than one candidate
query locus, or more than one ref CNEE sharing a query locus. Works per
query LOCUS (the merged query span of the filtering-eligible projections
landing there -- the "window"), scoring every chain aligned at that locus
once and judging every candidate projection there against the same set:

1. Soft depth filter: chains whose query envelope overlaps the window are
   ranked by their own genome-wide chain score (the chain header score, no
   new computation), and their blocks parsed in that order, keeping the
   first soft_depth_limit (default 20) that actually have aligned bases over
   the window. A projection whose chain doesn't make that cut is rejected
   without being scored (repeat-driven pileup). The scan visits at most
   SCAN_CAP_MULTIPLIER x soft_depth_limit chains, so cost per locus stays
   bounded even at pileups thousands of chains deep.
2. Exactly one kept chain ("isolated"): accept iff the projection's own
   flanks (flank_window each side, core excluded) have aligned-block
   coverage from ANY chain above flank_any_coverage_threshold -- a real
   ortholog's neighborhood is usually alignable even when this exact hit is
   uncontested.
3. Otherwise, each kept chain's local score = sequence identity over the
   window (matches / aligned bases) x aligned fraction over the window plus
   flank_window each side (aligned query bases / slice length; bases the
   chain doesn't reach count as unaligned). Relative score R = local score /
   best local score at the locus; a projection is retained iff
   R > relative_score_threshold (default 0.5). Several chains can clear the
   bar at one locus, so genuine many2one/one2many calls survive.

The vectorized block parsing here is what keeps a genome-wide run fast --
see parse_query_blocks_vectorized's docstring.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Sequence

import numpy as np

from .models import ChainBlock, ChainRecord
from .sequence_compare import AlignedSegment

DEFAULT_SOFT_DEPTH_LIMIT = 20
DEFAULT_FLANK_WINDOW = 10_000
DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD = 0.5
DEFAULT_RELATIVE_SCORE_THRESHOLD = 0.5
# score_locus parses at most this many x soft_depth_limit chains per locus
# (in chain-score order) while looking for soft_depth_limit chains with
# aligned bases over the window.
SCAN_CAP_MULTIPLIER = 10

# Outcome labels (pruned_query_projections.tsv's rejection_reason column).
REJECTED_HIGH_DEPTH_LOW_CHAIN_RANK = "high_depth_low_chain_rank"
ACCEPTED_ISOLATED_WITH_FLANK_SUPPORT = "accepted_isolated_with_flank_support"
REJECTED_ISOLATED_NO_FLANK_SUPPORT = "no_flank_support"
ACCEPTED_RELATIVE_SCORE = "accepted_relative_score"
REJECTED_LOW_RELATIVE_SCORE = "low_relative_score"
REJECTED_ZERO_LOCUS_SCORE = "zero_locus_score"


@dataclass(frozen=True)
class QueryBlocks:
    """A chain's aligned blocks, sorted by query position."""

    query_starts: np.ndarray
    query_ends: np.ndarray
    ref_starts: np.ndarray
    ref_ends: np.ndarray

    def _overlap_range(self, start: int, end: int) -> tuple[int, int]:
        if end <= start or self.query_starts.size == 0:
            return 0, 0
        lo = int(np.searchsorted(self.query_ends, start, side="right"))
        hi = int(np.searchsorted(self.query_starts, end, side="left"))
        return lo, max(lo, hi)

    def overlapping(self, start: int, end: int) -> list[ChainBlock]:
        lo, hi = self._overlap_range(start, end)
        return [
            ChainBlock(
                ref_start=int(self.ref_starts[i]), ref_end=int(self.ref_ends[i]),
                query_start=int(self.query_starts[i]), query_end=int(self.query_ends[i]),
            )
            for i in range(lo, hi)
        ]

    def aligned_bases(self, start: int, end: int) -> int:
        """Query bases in [start, end) covered by this chain's aligned blocks."""
        lo, hi = self._overlap_range(start, end)
        if hi <= lo:
            return 0
        clipped_starts = np.maximum(self.query_starts[lo:hi], start)
        clipped_ends = np.minimum(self.query_ends[lo:hi], end)
        return int(np.clip(clipped_ends - clipped_starts, 0, None).sum())


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
    """Per-query-chromosome chain envelope and chain-score arrays, for cheap overlap search."""
    by_chrom: dict[str, list[tuple[int, int, str, int]]] = {}
    for record in records:
        by_chrom.setdefault(record.query_chrom, []).append(
            (record.query_start, record.query_end, record.chain_id, record.score)
        )
    envelopes: dict[str, dict[str, np.ndarray]] = {}
    for chrom, rows in by_chrom.items():
        starts, ends, chain_ids, scores = zip(*rows)
        envelopes[chrom] = {
            "starts": np.array(starts, dtype=np.int64),
            "ends": np.array(ends, dtype=np.int64),
            "chain_ids": np.array(chain_ids),
            "scores": np.array(scores, dtype=np.int64),
        }
    return envelopes


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


def window_matches_for_chain(
    ref_accessor, query_accessor, ref_chrom: str, query_chrom: str, query_strand: str,
    blocks: QueryBlocks, window_start: int, window_end: int,
) -> tuple[int, int]:
    """(matches, compared bases) of one chain's aligned blocks over a query window."""
    segments = aligned_segments_for_query_interval(
        blocks.overlapping(window_start, window_end), window_start, window_end, query_strand
    )
    if not segments:
        return 0, 0

    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required for sequence comparison. Activate the biopython conda environment."
        ) from exc

    total_matches = 0
    total_compared = 0
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
        total_matches += int(np.count_nonzero(ref_arr == query_arr))
        total_compared += compared
    return total_matches, total_compared


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


class LocusProjection(Protocol):
    """What the filter needs from a candidate projection (see orthology.Projection)."""

    chain_id: str
    query_chrom: str
    query_start: int
    query_end: int


@dataclass(frozen=True)
class SliceScore:
    """One kept chain's alignment over a locus window and its flanks."""

    chain_id: str
    chain_score: int
    window_aligned_bases: int
    window_matches: int
    slice_aligned_bases: int


@dataclass
class LocusScores:
    query_chrom: str
    window_start: int
    window_end: int
    slice_start: int
    slice_end: int
    envelope_depth: int  # chains whose query envelope overlaps the window
    slices: list[SliceScore]  # kept chains, best chain score first
    # Only for isolated loci (exactly one slice): per locus projection (same
    # order), its any-chain flank coverage, or None if it isn't on the kept chain.
    isolated_flank_coverage: list[float | None] | None = None

    def local_score(self, slice_score: SliceScore) -> float:
        identity = (
            slice_score.window_matches / slice_score.window_aligned_bases
            if slice_score.window_aligned_bases else 0.0
        )
        slice_len = self.slice_end - self.slice_start
        aligned_fraction = slice_score.slice_aligned_bases / slice_len if slice_len > 0 else 0.0
        return identity * aligned_fraction


def score_locus(
    locus_projections: Sequence[LocusProjection],
    records_by_id: dict[str, ChainRecord],
    envelopes: dict[str, dict[str, np.ndarray]],
    cache: ChainEnvelopeCache,
    ref_accessor,
    query_accessor,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    soft_depth_limit: int = DEFAULT_SOFT_DEPTH_LIMIT,
) -> LocusScores:
    """Score every chain aligned over one query locus (see module docstring, steps 1-3)."""
    query_chrom = locus_projections[0].query_chrom
    window_start = min(p.query_start for p in locus_projections)
    window_end = max(p.query_end for p in locus_projections)
    query_size = records_by_id[locus_projections[0].chain_id].query_size
    slice_start = max(0, window_start - flank_window)
    slice_end = min(query_size, window_end + flank_window)

    envelope = envelopes[query_chrom]
    mask = (envelope["starts"] < window_end) & (envelope["ends"] > window_start)
    hit_ids = envelope["chain_ids"][mask]
    neg_scores = -envelope["scores"][mask]
    envelope_depth = int(hit_ids.size)

    # Only the top scan_cap by chain score can ever be visited, so partition
    # first instead of fully sorting a possibly ~10^5-long array; then order
    # by descending chain score, ties broken by envelope order.
    scan_cap = SCAN_CAP_MULTIPLIER * soft_depth_limit
    if envelope_depth > scan_cap:
        candidates = np.argpartition(neg_scores, scan_cap - 1)[:scan_cap]
    else:
        candidates = np.arange(envelope_depth)
    order = candidates[np.lexsort((candidates, neg_scores[candidates]))]

    slices: list[SliceScore] = []
    for idx in order:
        if len(slices) >= soft_depth_limit:
            break
        record, blocks = cache.get(str(hit_ids[idx]))
        if blocks.aligned_bases(window_start, window_end) == 0:
            continue
        matches, compared = window_matches_for_chain(
            ref_accessor, query_accessor, record.ref_chrom, query_chrom, record.query_strand,
            blocks, window_start, window_end,
        )
        slices.append(SliceScore(
            chain_id=record.chain_id, chain_score=record.score,
            window_aligned_bases=compared, window_matches=matches,
            slice_aligned_bases=blocks.aligned_bases(slice_start, slice_end),
        ))

    scores = LocusScores(
        query_chrom=query_chrom, window_start=window_start, window_end=window_end,
        slice_start=slice_start, slice_end=slice_end, envelope_depth=envelope_depth, slices=slices,
    )
    if len(slices) == 1:
        only_chain_id = slices[0].chain_id
        coverages: list[float | None] = []
        for p in locus_projections:
            if p.chain_id != only_chain_id:
                coverages.append(None)
                continue
            coverages.append(any_chain_flank_coverage(
                query_chrom, p.chain_id,
                max(0, p.query_start - flank_window), p.query_start,
                p.query_end, min(query_size, p.query_end + flank_window),
                envelopes, cache,
            ))
        scores.isolated_flank_coverage = coverages
    return scores


@dataclass(frozen=True)
class ProjectionVerdict:
    accepted: bool
    branch: str  # "soft_depth_filtered" | "isolated" | "competing"
    outcome: str
    locus_depth: int
    n_competing_chains: int
    local_score: float | None = None
    best_competitor_local_score: float | None = None
    relative_score: float | None = None
    flank_any_coverage: float | None = None


def decide_locus(
    locus_projections: Sequence[LocusProjection],
    scores: LocusScores,
    relative_score_threshold: float = DEFAULT_RELATIVE_SCORE_THRESHOLD,
    flank_any_coverage_threshold: float = DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
) -> list[ProjectionVerdict]:
    """One verdict per entry in locus_projections (same order), from score_locus's output."""
    n_kept = len(scores.slices)
    local_by_chain = {s.chain_id: scores.local_score(s) for s in scores.slices}
    s_max = max(local_by_chain.values(), default=0.0)

    verdicts: list[ProjectionVerdict] = []
    for i, p in enumerate(locus_projections):
        if p.chain_id not in local_by_chain:
            verdicts.append(ProjectionVerdict(
                accepted=False, branch="soft_depth_filtered", outcome=REJECTED_HIGH_DEPTH_LOW_CHAIN_RANK,
                locus_depth=scores.envelope_depth, n_competing_chains=n_kept,
            ))
            continue

        if n_kept == 1:
            coverage = scores.isolated_flank_coverage[i]
            accepted = coverage > flank_any_coverage_threshold
            verdicts.append(ProjectionVerdict(
                accepted=accepted, branch="isolated",
                outcome=ACCEPTED_ISOLATED_WITH_FLANK_SUPPORT if accepted else REJECTED_ISOLATED_NO_FLANK_SUPPORT,
                locus_depth=scores.envelope_depth, n_competing_chains=0, flank_any_coverage=coverage,
            ))
            continue

        own = local_by_chain[p.chain_id]
        best_other = max(v for cid, v in local_by_chain.items() if cid != p.chain_id)
        if s_max <= 0:
            verdicts.append(ProjectionVerdict(
                accepted=False, branch="competing", outcome=REJECTED_ZERO_LOCUS_SCORE,
                locus_depth=scores.envelope_depth, n_competing_chains=n_kept - 1,
                local_score=own, best_competitor_local_score=best_other,
            ))
            continue
        relative = own / s_max
        accepted = relative > relative_score_threshold
        verdicts.append(ProjectionVerdict(
            accepted=accepted, branch="competing",
            outcome=ACCEPTED_RELATIVE_SCORE if accepted else REJECTED_LOW_RELATIVE_SCORE,
            locus_depth=scores.envelope_depth, n_competing_chains=n_kept - 1,
            local_score=own, best_competitor_local_score=best_other, relative_score=relative,
        ))
    return verdicts
