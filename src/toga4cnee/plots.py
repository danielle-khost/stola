"""Dependency-free SVG QC plots for TOGA4CNEE outputs."""

from __future__ import annotations

import csv
import html
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PlotResult:
    histogram_path: Path
    overlap_density_path: Path
    sequence_identity_density_path: Path
    histogram_n: int
    overlap_density_n: int
    sequence_identity_density_n: int


def write_qc_plots(
    summary_path: str | Path,
    stats_path: str | Path,
    histogram_path: str | Path,
    overlap_density_path: str | Path,
    sequence_identity_density_path: str | Path,
) -> PlotResult:
    histogram_n = write_alignment_histogram(summary_path, histogram_path)
    overlap_density_n = write_overlap_density(stats_path, overlap_density_path)
    sequence_identity_density_n = write_sequence_identity_density(
        stats_path, sequence_identity_density_path
    )
    return PlotResult(
        histogram_path=Path(histogram_path),
        overlap_density_path=Path(overlap_density_path),
        sequence_identity_density_path=Path(sequence_identity_density_path),
        histogram_n=histogram_n,
        overlap_density_n=overlap_density_n,
        sequence_identity_density_n=sequence_identity_density_n,
    )


def write_alignment_histogram(summary_path: str | Path, out_path: str | Path) -> int:
    labels = [str(i) for i in range(11)] + [">10"]
    counts = [0 for _ in labels]
    total = 0
    with Path(summary_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "num_query_alignments" not in (reader.fieldnames or []):
            raise ValueError("Summary file is missing column: num_query_alignments")
        for row in reader:
            value = int(row["num_query_alignments"])
            counts[value if value <= 10 else 11] += 1
            total += 1

    svg = _bar_svg(
        title="Query alignments per reference CNEE",
        subtitle=f"Reference CNEEs: {total}",
        x_label="Number of query alignments",
        y_label="Reference CNEEs",
        labels=labels,
        counts=counts,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return total


def write_overlap_density(
    stats_path: str | Path, out_path: str | Path, bins: int = 100
) -> int:
    counts = [0 for _ in range(bins)]
    total = 0
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"aligned_bases", "covered_fraction"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            if int(row["aligned_bases"]) <= 0:
                continue
            value = float(row["covered_fraction"])
            value = min(1.0, max(0.0, value))
            index = min(bins - 1, int(value * bins))
            counts[index] += 1
            total += 1

    bin_width = 1.0 / bins
    densities = [count / (total * bin_width) if total else 0.0 for count in counts]
    smoothed = _moving_average(densities, radius=2)
    svg = _density_svg(
        title="Overlap fraction density",
        subtitle=f"Aligned candidate projections: {total}",
        x_label="CNEE overlap fraction",
        y_label="Density",
        values=smoothed,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return total


def write_sequence_identity_density(
    stats_path: str | Path, out_path: str | Path, bins: int = 100
) -> int:
    counts = [0 for _ in range(bins)]
    total = 0
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "seq_identity" not in (reader.fieldnames or []):
            total = 0
        else:
            for row in reader:
                value_text = row["seq_identity"]
                if value_text == ".":
                    continue
                value = min(1.0, max(0.0, float(value_text)))
                index = min(bins - 1, int(value * bins))
                counts[index] += 1
                total += 1

    bin_width = 1.0 / bins
    densities = [count / (total * bin_width) if total else 0.0 for count in counts]
    smoothed = _moving_average(densities, radius=2)
    svg = _density_svg(
        title="Sequence identity density",
        subtitle=f"Sequence-compared projections: {total}",
        x_label="Sequence identity",
        y_label="Density",
        values=smoothed,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return total


def _moving_average(values: list[float], radius: int) -> list[float]:
    if not values:
        return []
    output = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        output.append(sum(values[start:end]) / (end - start))
    return output


def _scale_y(value: float, max_value: float, plot_top: float, plot_height: float) -> float:
    if max_value <= 0:
        return plot_top + plot_height
    return plot_top + plot_height - (value / max_value) * plot_height


def _svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        "<style>"
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}"
        ".title{font-size:24px;font-weight:700}"
        ".subtitle{font-size:13px;fill:#52606d}"
        ".axis{stroke:#1f2933;stroke-width:1.2}"
        ".grid{stroke:#d9e2ec;stroke-width:1}"
        ".label{font-size:13px}"
        ".tick{font-size:11px;fill:#52606d}"
        "</style>\n"
    )


def _bar_svg(
    title: str,
    subtitle: str,
    x_label: str,
    y_label: str,
    labels: list[str],
    counts: list[int],
) -> str:
    width, height = 900, 520
    left, right, top, bottom = 80, 35, 75, 80
    plot_width = width - left - right
    plot_height = height - top - bottom
    max_count = max(counts) if counts else 0
    tick_count = 5
    slot = plot_width / len(labels)
    bar_width = slot * 0.72

    parts = [_svg_header(width, height)]
    parts.append(f'<text class="title" x="{left}" y="34">{html.escape(title)}</text>\n')
    parts.append(
        f'<text class="subtitle" x="{left}" y="55">{html.escape(subtitle)}</text>\n'
    )
    for tick in range(tick_count + 1):
        value = max_count * tick / tick_count if tick_count else 0
        y = _scale_y(value, max_count, top, plot_height)
        parts.append(
            f'<line class="grid" x1="{left}" y1="{y:.2f}" '
            f'x2="{width - right}" y2="{y:.2f}"/>\n'
        )
        parts.append(
            f'<text class="tick" x="{left - 8}" y="{y + 4:.2f}" '
            f'text-anchor="end">{int(round(value))}</text>\n'
        )
    parts.append(
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}"/>\n'
    )
    parts.append(
        f'<line class="axis" x1="{left}" y1="{top + plot_height}" '
        f'x2="{width - right}" y2="{top + plot_height}"/>\n'
    )
    for index, (label, count) in enumerate(zip(labels, counts)):
        x = left + index * slot + (slot - bar_width) / 2
        y = _scale_y(count, max_count, top, plot_height)
        bar_height = top + plot_height - y
        parts.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" '
            f'height="{bar_height:.2f}" fill="#2f80ed"/>\n'
        )
        parts.append(
            f'<text class="tick" x="{x + bar_width / 2:.2f}" y="{top + plot_height + 22}" '
            f'text-anchor="middle">{html.escape(label)}</text>\n'
        )
    parts.append(
        f'<text class="label" x="{left + plot_width / 2}" y="{height - 24}" '
        f'text-anchor="middle">{html.escape(x_label)}</text>\n'
    )
    parts.append(
        f'<text class="label" transform="translate(22 {top + plot_height / 2}) rotate(-90)" '
        f'text-anchor="middle">{html.escape(y_label)}</text>\n'
    )
    parts.append("</svg>\n")
    return "".join(parts)


def _density_svg(
    title: str,
    subtitle: str,
    x_label: str,
    y_label: str,
    values: list[float],
) -> str:
    width, height = 900, 520
    left, right, top, bottom = 80, 35, 75, 80
    plot_width = width - left - right
    plot_height = height - top - bottom
    max_value = max(values) if values else 0.0
    tick_count = 5

    parts = [_svg_header(width, height)]
    parts.append(f'<text class="title" x="{left}" y="34">{html.escape(title)}</text>\n')
    parts.append(
        f'<text class="subtitle" x="{left}" y="55">{html.escape(subtitle)}</text>\n'
    )
    for tick in range(tick_count + 1):
        value = max_value * tick / tick_count if tick_count else 0.0
        y = _scale_y(value, max_value, top, plot_height)
        parts.append(
            f'<line class="grid" x1="{left}" y1="{y:.2f}" '
            f'x2="{width - right}" y2="{y:.2f}"/>\n'
        )
        parts.append(
            f'<text class="tick" x="{left - 8}" y="{y + 4:.2f}" '
            f'text-anchor="end">{value:.2g}</text>\n'
        )
    for tick in range(tick_count + 1):
        value = tick / tick_count
        x = left + value * plot_width
        parts.append(
            f'<text class="tick" x="{x:.2f}" y="{top + plot_height + 22}" '
            f'text-anchor="middle">{value:.1f}</text>\n'
        )
    parts.append(
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}"/>\n'
    )
    parts.append(
        f'<line class="axis" x1="{left}" y1="{top + plot_height}" '
        f'x2="{width - right}" y2="{top + plot_height}"/>\n'
    )
    if values:
        points = []
        for index, value in enumerate(values):
            denom = max(1, len(values) - 1)
            x = left + (index / denom) * plot_width
            y = _scale_y(value, max_value, top, plot_height)
            points.append((x, y))
        polyline = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        baseline_y = top + plot_height
        polygon = (
            f"{left:.2f},{baseline_y:.2f} "
            + polyline
            + f" {width - right:.2f},{baseline_y:.2f}"
        )
        parts.append(f'<polygon points="{polygon}" fill="#b9dcff" opacity="0.7"/>\n')
        parts.append(
            f'<polyline points="{polyline}" fill="none" stroke="#1f6feb" '
            'stroke-width="2.5"/>\n'
        )
    parts.append(
        f'<text class="label" x="{left + plot_width / 2}" y="{height - 24}" '
        f'text-anchor="middle">{html.escape(x_label)}</text>\n'
    )
    parts.append(
        f'<text class="label" transform="translate(22 {top + plot_height / 2}) rotate(-90)" '
        f'text-anchor="middle">{html.escape(y_label)}</text>\n'
    )
    parts.append("</svg>\n")
    return "".join(parts)
