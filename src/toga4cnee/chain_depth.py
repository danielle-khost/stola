"""Per-base ALIGNED-BLOCK depth across a chain file's reference and query
genomes.

"Aligned-block depth" at a base = the number of chains whose own *aligned
block* (the ungapped intervals chain.py's iter_chain_blocks extracts from
chain body lines, honoring the dt/dq gap fields between blocks -- NOT the
whole-chain [start, end) bounding span, which would count a chain's internal
gaps as "covered") overlaps that base. Computed independently for the
reference and query side of the same chain file in one streaming pass over
chain bodies (mirroring iter_chain_blocks' coordinate walk for both sides at
once, since a chain's ref and query chromosome/contig sets are independent
axes but their block boundaries advance together per body line).

Needs numpy for the per-chromosome diff/cumsum sweep -- the only place in
the toga4cnee package that does, since every other summary here operates on
already-aggregated per-candidate rows rather than per-base arrays spanning
whole chromosomes (order 10^8-10^9 positions each).
"""

from __future__ import annotations

import gzip
from array import array
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DepthStats:
    n_total_bases: int
    mean_depth: float
    q1_depth: float
    median_depth: float
    q3_depth: float
    min_depth: int
    max_depth: int


@dataclass(frozen=True)
class AlignedBlockDepthResult:
    ref_histogram: dict[int, int]  # depth -> # reference bases at that depth
    query_histogram: dict[int, int]  # depth -> # query bases at that depth
    n_chains: int
    n_blocks: int


def _open_chain_text(chain_path: str | Path):
    chain_path = Path(chain_path)
    if chain_path.suffix == ".gz":
        return gzip.open(chain_path, "rt", encoding="utf-8")
    return chain_path.open("r", encoding="utf-8")


def _both_side_blocks(
    header_fields: list[str], body_lines: list[str]
) -> tuple[str, int, list[int], list[int], str, int, list[int], list[int]]:
    # Mirrors chain.py's iter_chain_blocks coordinate walk, computing both
    # the reference-side AND query-side aligned-block [start, end) intervals
    # (normalized to plus-strand coordinates) in the same pass over body
    # lines, since both progress counters advance together per block.
    ref_chrom = header_fields[2]
    ref_size = int(header_fields[3])
    ref_strand = header_fields[4]
    ref_raw_start = int(header_fields[5])
    query_chrom = header_fields[7]
    query_size = int(header_fields[8])
    query_strand = header_fields[9]
    query_raw_start = int(header_fields[10])

    ref_progress = ref_raw_start if ref_strand == "+" else ref_size - ref_raw_start
    query_progress = query_raw_start if query_strand == "+" else query_size - query_raw_start

    ref_starts: list[int] = []
    ref_ends: list[int] = []
    query_starts: list[int] = []
    query_ends: list[int] = []
    for line in body_lines:
        parts = line.split()
        size = int(parts[0])
        dt = int(parts[1]) if len(parts) == 3 else 0
        dq = int(parts[2]) if len(parts) == 3 else 0

        if ref_strand == "+":
            r_start = ref_progress
            r_end = r_start + size
            ref_progress = r_end + dt
        else:
            r_end = ref_progress
            r_start = r_end - size
            ref_progress = r_start - dt
        ref_starts.append(r_start)
        ref_ends.append(r_end)

        if query_strand == "+":
            q_start = query_progress
            q_end = q_start + size
            query_progress = q_end + dq
        else:
            q_end = query_progress
            q_start = q_end - size
            query_progress = q_start - dq
        query_starts.append(q_start)
        query_ends.append(q_end)

    return ref_chrom, ref_size, ref_starts, ref_ends, query_chrom, query_size, query_starts, query_ends


def _histogram_from_spans(
    chrom_sizes: dict[str, int],
    starts_by_chrom: dict[str, array],
    ends_by_chrom: dict[str, array],
) -> dict[int, int]:
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "numpy is required for aligned-block depth computation. "
            "Activate the biopython conda environment or install numpy."
        ) from exc

    histogram: dict[int, int] = defaultdict(int)
    for chrom, size in chrom_sizes.items():
        starts = np.frombuffer(starts_by_chrom[chrom], dtype=np.int64)
        ends = np.frombuffer(ends_by_chrom[chrom], dtype=np.int64)

        diff = np.zeros(size + 1, dtype=np.int32)
        np.add.at(diff, starts, 1)
        np.add.at(diff, ends, -1)
        depth = np.cumsum(diff, out=diff)[:size]

        counts = np.bincount(depth)
        for d, c in enumerate(counts):
            if c:
                histogram[int(d)] += int(c)
    return dict(histogram)


def compute_aligned_block_depth_histograms(chain_path: str | Path) -> AlignedBlockDepthResult:
    """Stream a chain file once, returning per-base aligned-block depth
    histograms for both its reference and query genomes."""
    ref_chrom_sizes: dict[str, int] = {}
    ref_starts_by_chrom: defaultdict[str, array] = defaultdict(lambda: array("q"))
    ref_ends_by_chrom: defaultdict[str, array] = defaultdict(lambda: array("q"))
    query_chrom_sizes: dict[str, int] = {}
    query_starts_by_chrom: defaultdict[str, array] = defaultdict(lambda: array("q"))
    query_ends_by_chrom: defaultdict[str, array] = defaultdict(lambda: array("q"))
    n_chains = 0
    n_blocks = 0

    header_fields: list[str] | None = None
    body_lines: list[str] = []

    def flush() -> None:
        nonlocal n_chains, n_blocks
        if header_fields is None:
            return
        ref_chrom, ref_size, ref_starts, ref_ends, query_chrom, query_size, query_starts, query_ends = (
            _both_side_blocks(header_fields, body_lines)
        )
        ref_chrom_sizes.setdefault(ref_chrom, ref_size)
        ref_starts_by_chrom[ref_chrom].extend(ref_starts)
        ref_ends_by_chrom[ref_chrom].extend(ref_ends)
        query_chrom_sizes.setdefault(query_chrom, query_size)
        query_starts_by_chrom[query_chrom].extend(query_starts)
        query_ends_by_chrom[query_chrom].extend(query_ends)
        n_chains += 1
        n_blocks += len(ref_starts)

    with _open_chain_text(chain_path) as handle:
        for line in handle:
            if line.startswith("chain "):
                flush()
                header_fields = line.split()
                body_lines = []
            elif line.strip():
                body_lines.append(line)
        flush()

    ref_histogram = _histogram_from_spans(ref_chrom_sizes, ref_starts_by_chrom, ref_ends_by_chrom)
    query_histogram = _histogram_from_spans(query_chrom_sizes, query_starts_by_chrom, query_ends_by_chrom)
    return AlignedBlockDepthResult(ref_histogram, query_histogram, n_chains, n_blocks)


def _weighted_quantile(sorted_items: list[tuple[int, int]], total: int, q: float) -> float:
    threshold = total * q
    cumulative = 0
    for depth, count in sorted_items:
        cumulative += count
        if cumulative >= threshold:
            return float(depth)
    return float(sorted_items[-1][0])


def compute_depth_stats(histogram: dict[int, int]) -> DepthStats:
    sorted_items = sorted(histogram.items())
    total = sum(count for _, count in sorted_items)
    mean = sum(depth * count for depth, count in sorted_items) / total
    return DepthStats(
        n_total_bases=total,
        mean_depth=mean,
        q1_depth=_weighted_quantile(sorted_items, total, 0.25),
        median_depth=_weighted_quantile(sorted_items, total, 0.5),
        q3_depth=_weighted_quantile(sorted_items, total, 0.75),
        min_depth=sorted_items[0][0],
        max_depth=sorted_items[-1][0],
    )
