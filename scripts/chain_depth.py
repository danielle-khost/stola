#!/usr/bin/env python
"""Per-base ALIGNED-BLOCK depth across a chain file's reference and query
genomes, as a 2-panel density SVG plus mean/quartile stats.

    python scripts/chain_depth.py alignment.chain.gz --out chain_depth.svg [--bin-max 100]

"Aligned-block depth" at a base = the number of chains whose own *aligned
block* (the ungapped intervals chain.py's iter_chain_blocks extracts from
chain body lines, honoring the dt/dq gap fields between blocks -- NOT the
whole-chain [start, end) bounding span, which would count a chain's internal
gaps as "covered") overlaps that base. Computed independently for the
reference and query side of the same chain file in one streaming pass over
chain bodies (mirroring iter_chain_blocks' coordinate walk for both sides at
once, since a chain's ref and query chromosome/contig sets are independent
axes but their block boundaries advance together per body line).

Needs numpy for the per-chromosome diff/cumsum sweep over per-base arrays
spanning whole chromosomes (order 10^8-10^9 positions each).
"""

from __future__ import annotations

import argparse
import gzip
import html
import math
import sys
from array import array
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from tunic.plots import _svg_header  # noqa: E402


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


@dataclass(frozen=True)
class AlignedBlockDepthPlotResult:
    path: Path


def _log_scale_y(value: float, min_exp: float, max_exp: float, plot_top: float, plot_height: float) -> float:
    # Like _scale_y but in log10 space; a non-positive value has no bar at
    # all (same convention as a log-scale axis silently dropping it).
    if value <= 0:
        return plot_top + plot_height
    exponent = min(max_exp, max(min_exp, math.log10(value)))
    frac = (exponent - min_exp) / (max_exp - min_exp) if max_exp > min_exp else 1.0
    return plot_top + plot_height - frac * plot_height


def _depth_panel_svg(
    histogram: dict[int, int],
    stats: dict[str, float],
    bin_max: int,
    title: str,
    x0: float,
    y0: float,
    width: float,
    height: float,
) -> str:
    # One log-scale depth-density bar chart (depths 0..bin_max plus one
    # ">bin_max" overflow bin) with a small inset mean/quartile table,
    # translated to (x0, y0) within the caller's combined SVG canvas.
    left, right, top, bottom = 60, 20, 36, 46
    plot_width = width - left - right
    plot_height = height - top - bottom

    total = sum(histogram.values())
    fractions = [histogram.get(d, 0) / total for d in range(bin_max + 1)]
    overflow_fraction = sum(c for d, c in histogram.items() if d > bin_max) / total

    positive = [f for f in fractions + [overflow_fraction] if f > 0]
    max_exp = 0.0
    min_exp = math.floor(math.log10(min(positive))) if positive else -1.0

    slot = plot_width / (bin_max + 2)
    bar_width = slot * 0.85

    parts = [f'<g transform="translate({x0} {y0})">\n']
    parts.append(f'<text class="label" x="{left}" y="20" font-weight="700">{html.escape(title)}</text>\n')

    for exponent in range(int(min_exp), 1):
        value = 10.0**exponent
        y = _log_scale_y(value, min_exp, max_exp, top, plot_height)
        parts.append(
            f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}"/>\n'
        )
        label = f"1e{exponent}" if exponent != 0 else "1"
        parts.append(
            f'<text class="tick" x="{left - 8}" y="{y + 4:.2f}" text-anchor="end">{label}</text>\n'
        )

    baseline = top + plot_height
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{baseline}"/>\n')
    parts.append(f'<line class="axis" x1="{left}" y1="{baseline}" x2="{width - right}" y2="{baseline}"/>\n')

    tick_step = max(bin_max // 10, 1)
    for depth in list(range(0, bin_max + 1)) + [bin_max + 1]:
        value = overflow_fraction if depth == bin_max + 1 else fractions[depth]
        x = left + depth * slot + (slot - bar_width) / 2
        y = _log_scale_y(value, min_exp, max_exp, top, plot_height)
        if y < baseline:
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" '
                f'height="{baseline - y:.2f}" fill="#1f6feb"/>\n'
            )
        # Skip bin_max itself (sits right next to the ">bin_max" overflow
        # tick and would collide with it); every other step gets a label.
        if depth == bin_max + 1:
            label = f">{bin_max}"
        elif depth % tick_step == 0 and depth != bin_max:
            label = str(depth)
        else:
            continue
        parts.append(
            f'<text class="tick" x="{x + bar_width / 2:.2f}" y="{baseline + 16}" '
            f'text-anchor="middle">{label}</text>\n'
        )

    parts.append(
        f'<text class="label" x="{left + plot_width / 2}" y="{height - 6}" '
        'text-anchor="middle">Aligned-block depth (# overlapping aligned chain blocks)</text>\n'
    )

    table_rows = [
        ("Mean", f"{stats['mean_depth']:.3f}"),
        ("Q1", f"{stats['q1_depth']:.0f}"),
        ("Median", f"{stats['median_depth']:.0f}"),
        ("Q3", f"{stats['q3_depth']:.0f}"),
    ]
    table_x = width - right - 150
    table_y = top + 6
    parts.append(
        f'<rect x="{table_x}" y="{table_y}" width="150" height="{18 * (len(table_rows) + 1)}" '
        'fill="#f8fafc" stroke="#d9e2ec" stroke-width="1"/>\n'
    )
    parts.append(f'<text class="table-header" x="{table_x + 8}" y="{table_y + 15}">Stat</text>\n')
    parts.append(f'<text class="table-header" x="{table_x + 142}" y="{table_y + 15}" text-anchor="end">Depth</text>\n')
    for row_index, (label, value) in enumerate(table_rows, start=1):
        row_y = table_y + 15 + row_index * 18
        parts.append(f'<text class="table-text" x="{table_x + 8}" y="{row_y}">{label}</text>\n')
        parts.append(
            f'<text class="table-text" x="{table_x + 142}" y="{row_y}" text-anchor="end">{value}</text>\n'
        )
    parts.append("</g>\n")
    return "".join(parts)


def write_aligned_block_depth_density(
    ref_histogram: dict[int, int],
    ref_stats: dict[str, float],
    query_histogram: dict[int, int],
    query_stats: dict[str, float],
    out_path: str | Path,
    chain_label: str = "",
    bin_max: int = 100,
) -> AlignedBlockDepthPlotResult:
    """Write a 2-panel (reference top, query bottom) per-base aligned-block
    depth density plot, each depth above bin_max collapsed into one overflow
    bin, with an inset mean/Q1/median/Q3 table. ref_stats/query_stats are
    plain dicts with those four keys (asdict of compute_depth_stats)."""
    width = 1100
    panel_height = 420
    title_height = 44
    gap = 24
    height = title_height + 2 * panel_height + gap

    parts = [_svg_header(width, height)]
    if chain_label:
        parts.append(
            f'<text class="title" x="30" y="30" font-size="18">'
            f'Per-base aligned-block depth -- {html.escape(chain_label)}</text>\n'
        )
    ref_total = ref_stats.get("n_total_bases", sum(ref_histogram.values()))
    query_total = query_stats.get("n_total_bases", sum(query_histogram.values()))
    parts.append(
        _depth_panel_svg(
            ref_histogram, ref_stats, bin_max, f"Reference genome (n={ref_total:,} bp)",
            0, title_height, width, panel_height,
        )
    )
    parts.append(
        _depth_panel_svg(
            query_histogram, query_stats, bin_max, f"Query genome (n={query_total:,} bp)",
            0, title_height + panel_height + gap, width, panel_height,
        )
    )
    parts.append("</svg>\n")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(parts), encoding="utf-8")
    return AlignedBlockDepthPlotResult(path=out_path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Per-base aligned-block depth density plot for a chain file's ref/query genomes"
    )
    parser.add_argument("chain_file", type=Path, help="Input .chain or .chain.gz file")
    parser.add_argument("--out", required=True, type=Path, help="Output SVG path")
    parser.add_argument(
        "--bin-max",
        type=int,
        default=100,
        help="Collapse every depth above this into one overflow bin [default: 100]",
    )
    args = parser.parse_args(argv)

    result = compute_aligned_block_depth_histograms(args.chain_file)
    print(f"{result.n_chains} chains, {result.n_blocks} aligned blocks")
    ref_stats = compute_depth_stats(result.ref_histogram)
    query_stats = compute_depth_stats(result.query_histogram)
    print(f"Reference: {ref_stats}")
    print(f"Query: {query_stats}")

    plot_result = write_aligned_block_depth_density(
        result.ref_histogram,
        asdict(ref_stats),
        result.query_histogram,
        asdict(query_stats),
        args.out,
        chain_label=args.chain_file.name,
        bin_max=args.bin_max,
    )
    print(f"Wrote aligned-block depth density plot -> {plot_result.path}")


if __name__ == "__main__":
    main()
