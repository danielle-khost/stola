"""Dependency-free SVG QC plots for TOGA4CNEE outputs."""

from __future__ import annotations

import csv
import html
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

CATEGORY_ORDER = [
    "one2one",
    "one2many",
    "many2one",
    "complex/many2many",
    "uncertain_loss",
    "lost",
    "missing",
]

POSITIVE_CATEGORIES = {
    "one2one",
    "one2many",
    "many2one",
    "complex/many2many",
}


@dataclass(frozen=True)
class PlotResult:
    histogram_path: Path
    overlap_density_path: Path
    sequence_identity_density_path: Path
    histogram_n: int
    overlap_density_n: int
    sequence_identity_density_n: int


@dataclass(frozen=True)
class OrthologyDiagnosticPlotResult:
    category_counts_path: Path
    sequence_identity_by_category_path: Path
    gap_fraction_by_category_path: Path
    category_counts: dict[str, int]
    sequence_identity_counts: dict[str, int]
    gap_fraction_counts: dict[str, int]


@dataclass(frozen=True)
class MultiAlignmentMetricPlotResult:
    sequence_identity_density_path: Path
    gap_fraction_density_path: Path
    multi_alignment_ref_count: int
    sequence_identity_density_n: int
    gap_fraction_density_n: int


def write_qc_plots(
    summary_path: str | Path,
    stats_path: str | Path,
    histogram_path: str | Path,
    overlap_density_path: str | Path,
    sequence_identity_density_path: str | Path,
) -> PlotResult:
    # Write the standard trio of SVG QC plots for a completed workflow run.
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


def write_orthology_diagnostic_plots(
    stats_path: str | Path,
    orthology_path: str | Path,
    category_counts_path: str | Path,
    sequence_identity_by_category_path: str | Path,
    gap_fraction_by_category_path: str | Path,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
) -> OrthologyDiagnosticPlotResult:
    # Write category counts plus grouped sequence-identity/gap-fraction densities.
    category_counts = write_orthology_category_counts(
        orthology_path, category_counts_path
    )
    grouped_values = _collect_orthology_metric_values(
        stats_path,
        orthology_path,
        min_sequence_identity=min_sequence_identity,
        max_fraction_gap=max_fraction_gap,
    )
    seq_counts = {
        category: len(values)
        for category, values in grouped_values["sequence_identity"].items()
    }
    gap_counts = {
        category: len(values)
        for category, values in grouped_values["gap_fraction"].items()
    }
    _write_grouped_density(
        grouped_values["sequence_identity"],
        sequence_identity_by_category_path,
        title="Sequence identity by orthology category",
        subtitle="Solid orthology projections plus lost/uncertain-loss candidates",
        x_label="Sequence identity",
    )
    _write_grouped_density(
        grouped_values["gap_fraction"],
        gap_fraction_by_category_path,
        title="Gap fraction by orthology category",
        subtitle="Solid orthology projections plus lost/uncertain-loss candidates",
        x_label="Gap fraction",
    )
    return OrthologyDiagnosticPlotResult(
        category_counts_path=Path(category_counts_path),
        sequence_identity_by_category_path=Path(sequence_identity_by_category_path),
        gap_fraction_by_category_path=Path(gap_fraction_by_category_path),
        category_counts=category_counts,
        sequence_identity_counts=seq_counts,
        gap_fraction_counts=gap_counts,
    )


def write_orthology_category_counts(
    orthology_path: str | Path,
    out_path: str | Path,
) -> dict[str, int]:
    # Plot orthology TSV row counts by category.
    counts: Counter[str] = Counter()
    with Path(orthology_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "orthology_category" not in (reader.fieldnames or []):
            raise ValueError("Orthology file is missing column: orthology_category")
        for row in reader:
            counts[row["orthology_category"]] += 1

    labels = [category for category in CATEGORY_ORDER if counts[category] > 0]
    labels.extend(sorted(category for category in counts if category not in labels))
    svg = _bar_svg(
        title="Orthology category counts",
        subtitle=f"Orthology rows: {sum(counts.values())}",
        x_label="Orthology category",
        y_label="Rows",
        labels=labels,
        counts=[counts[label] for label in labels],
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return dict(counts)


def write_alignment_histogram(summary_path: str | Path, out_path: str | Path) -> int:
    # Plot the distribution of query-alignment counts per reference CNEE.
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
    # Plot a smoothed density of aligned projection overlap fractions.
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
    stats_path: str | Path,
    out_path: str | Path,
    bins: int = 100,
    ref_ids: set[str] | None = None,
    title: str = "Sequence identity density",
    subtitle_label: str = "Sequence-compared projections",
) -> int:
    # Plot a smoothed density of sequence identities when sequence metrics exist.
    counts = [0 for _ in range(bins)]
    total = 0
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cnee_id", "seq_identity"} if ref_ids is not None else {"seq_identity"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            if ref_ids is not None and row["cnee_id"] not in ref_ids:
                continue
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
        title=title,
        subtitle=f"{subtitle_label}: {total}",
        x_label="Sequence identity",
        y_label="Density",
        values=smoothed,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return total


def write_gap_fraction_density(
    stats_path: str | Path,
    out_path: str | Path,
    bins: int = 100,
    ref_ids: set[str] | None = None,
    title: str = "Gap fraction density",
    subtitle_label: str = "Candidate projections",
) -> int:
    # Plot a smoothed density of chain-span gap fractions for candidate rows.
    counts = [0 for _ in range(bins)]
    total = 0
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = (
            {"cnee_id", "gap_fraction_within_chain_span"}
            if ref_ids is not None
            else {"gap_fraction_within_chain_span"}
        )
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            if ref_ids is not None and row["cnee_id"] not in ref_ids:
                continue
            value = min(1.0, max(0.0, float(row["gap_fraction_within_chain_span"])))
            index = min(bins - 1, int(value * bins))
            counts[index] += 1
            total += 1

    bin_width = 1.0 / bins
    densities = [count / (total * bin_width) if total else 0.0 for count in counts]
    smoothed = _moving_average(densities, radius=2)
    svg = _density_svg(
        title=title,
        subtitle=f"{subtitle_label}: {total}",
        x_label="Gap fraction",
        y_label="Density",
        values=smoothed,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return total


def write_multi_alignment_metric_densities(
    summary_path: str | Path,
    stats_path: str | Path,
    sequence_identity_path: str | Path,
    gap_fraction_path: str | Path,
) -> MultiAlignmentMetricPlotResult:
    # Plot sequence/gap metrics only for refs with more than one query alignment.
    ref_ids = _read_multi_alignment_ref_ids(summary_path)
    seq_n = write_sequence_identity_density(
        stats_path,
        sequence_identity_path,
        ref_ids=ref_ids,
        title="Sequence identity density for multi-query CNEEs",
        subtitle_label="Sequence-compared projections from refs with >1 query alignment",
    )
    gap_n = write_gap_fraction_density(
        stats_path,
        gap_fraction_path,
        ref_ids=ref_ids,
        title="Gap fraction density for multi-query CNEEs",
        subtitle_label="Candidate projections from refs with >1 query alignment",
    )
    return MultiAlignmentMetricPlotResult(
        sequence_identity_density_path=Path(sequence_identity_path),
        gap_fraction_density_path=Path(gap_fraction_path),
        multi_alignment_ref_count=len(ref_ids),
        sequence_identity_density_n=seq_n,
        gap_fraction_density_n=gap_n,
    )


def _read_multi_alignment_ref_ids(summary_path: str | Path) -> set[str]:
    # Read summary rows whose reference CNEE has more than one query alignment.
    ref_ids: set[str] = set()
    with Path(summary_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cnee_id", "num_query_alignments"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Summary file is missing columns: {', '.join(missing)}")
        for row in reader:
            if int(row["num_query_alignments"]) > 1:
                ref_ids.add(row["cnee_id"])
    return ref_ids


def _collect_orthology_metric_values(
    stats_path: str | Path,
    orthology_path: str | Path,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> dict[str, dict[str, list[float]]]:
    # Assign candidate projection metrics to ref-level orthology categories.
    ref_categories = _read_ref_categories(orthology_path)
    sequence_identity_values: dict[str, list[float]] = defaultdict(list)
    gap_fraction_values: dict[str, list[float]] = defaultdict(list)

    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_id",
            "aligned_bases",
            "gap_fraction_within_chain_span",
            "seq_identity",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            category = ref_categories.get(row["cnee_id"])
            if category is None or category == "missing":
                continue
            if category in POSITIVE_CATEGORIES and not _stats_row_is_solid(
                row,
                min_sequence_identity=min_sequence_identity,
                max_fraction_gap=max_fraction_gap,
            ):
                continue

            gap_fraction_values[category].append(
                min(1.0, max(0.0, float(row["gap_fraction_within_chain_span"])))
            )
            seq_identity = row.get("seq_identity", ".")
            if seq_identity != "." and seq_identity != "":
                sequence_identity_values[category].append(
                    min(1.0, max(0.0, float(seq_identity)))
                )

    return {
        "sequence_identity": dict(sequence_identity_values),
        "gap_fraction": dict(gap_fraction_values),
    }


def _read_ref_categories(orthology_path: str | Path) -> dict[str, str]:
    # Map each reference CNEE to the category assigned in the orthology TSV.
    categories_by_ref: dict[str, set[str]] = defaultdict(set)
    with Path(orthology_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"partner_ref_cnees", "orthology_category"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Orthology file is missing columns: {', '.join(missing)}")
        for row in reader:
            category = row["orthology_category"]
            for ref_id in row["partner_ref_cnees"].split(","):
                if ref_id:
                    categories_by_ref[ref_id].add(category)
    return {
        ref_id: _choose_ref_category(categories)
        for ref_id, categories in categories_by_ref.items()
    }


def _choose_ref_category(categories: set[str]) -> str:
    # Collapse repeated one2many rows or complex components to one plotting group.
    if "complex/many2many" in categories:
        return "complex/many2many"
    for category in CATEGORY_ORDER:
        if category in categories:
            return category
    return sorted(categories)[0]


def _stats_row_is_solid(
    row: dict[str, str],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # Reapply solid-edge criteria so positive-category densities use solid rows.
    if int(row["aligned_bases"]) <= 0:
        return False
    if float(row["gap_fraction_within_chain_span"]) >= max_fraction_gap:
        return False
    seq_identity = row.get("seq_identity", ".")
    if seq_identity == "." or seq_identity == "":
        return True
    return float(seq_identity) > min_sequence_identity


def _write_grouped_density(
    grouped_values: dict[str, list[float]],
    out_path: str | Path,
    title: str,
    subtitle: str,
    x_label: str,
    bins: int = 100,
) -> dict[str, int]:
    # Write a multi-category density plot and return plotted counts by category.
    series: list[tuple[str, list[float], int]] = []
    for category in CATEGORY_ORDER:
        values = grouped_values.get(category, [])
        if values:
            series.append((category, _density_values(values, bins), len(values)))
    for category in sorted(set(grouped_values).difference(CATEGORY_ORDER)):
        values = grouped_values[category]
        if values:
            series.append((category, _density_values(values, bins), len(values)))

    svg = _grouped_density_svg(
        title=title,
        subtitle=subtitle,
        x_label=x_label,
        y_label="Density",
        series=series,
    )
    Path(out_path).write_text(svg, encoding="utf-8")
    return {category: count for category, _, count in series}


def _density_values(values: list[float], bins: int) -> list[float]:
    # Convert values from [0, 1] into a smoothed density series.
    counts = [0 for _ in range(bins)]
    for value in values:
        value = min(1.0, max(0.0, value))
        index = min(bins - 1, int(value * bins))
        counts[index] += 1
    bin_width = 1.0 / bins
    densities = [count / (len(values) * bin_width) if values else 0.0 for count in counts]
    return _moving_average(densities, radius=2)


def _moving_average(values: list[float], radius: int) -> list[float]:
    # Smooth binned densities with a centered moving average.
    if not values:
        return []
    output = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        output.append(sum(values[start:end]) / (end - start))
    return output


def _scale_y(value: float, max_value: float, plot_top: float, plot_height: float) -> float:
    # Convert a data value into an SVG y coordinate within the plot area.
    if max_value <= 0:
        return plot_top + plot_height
    return plot_top + plot_height - (value / max_value) * plot_height


def _svg_header(width: int, height: int) -> str:
    # Create the common SVG shell and inline styles for QC plots.
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
        ".table-header{font-size:12px;font-weight:700;fill:#1f2933}"
        ".table-text{font-size:11px;fill:#1f2933}"
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
    # Render a labeled bar chart and a compact count table as standalone SVG.
    width, height = 1120, 540
    left, right, top, bottom = 80, 255, 75, 92
    table_left = width - right + 32
    table_top = top + 2
    table_width = right - 55
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
    parts.append(
        f'<rect x="{table_left - 14}" y="{table_top - 22}" '
        f'width="{table_width}" height="{plot_height + 48}" '
        'fill="#f8fafc" stroke="#d9e2ec" stroke-width="1"/>\n'
    )
    parts.append(
        f'<text class="table-header" x="{table_left}" y="{table_top}">Counts</text>\n'
    )
    row_y = table_top + 22
    for label, count in zip(labels, counts):
        parts.append(
            f'<text class="table-text" x="{table_left}" y="{row_y}">'
            f'{html.escape(label)}</text>\n'
        )
        parts.append(
            f'<text class="table-text" x="{width - 36}" y="{row_y}" text-anchor="end">'
            f'{count:,}</text>\n'
        )
        row_y += 18
    parts.append(
        f'<line x1="{table_left}" y1="{row_y - 8}" x2="{width - 36}" '
        f'y2="{row_y - 8}" stroke="#d9e2ec" stroke-width="1"/>\n'
    )
    parts.append(
        f'<text class="table-header" x="{table_left}" y="{row_y + 8}">Total</text>\n'
    )
    parts.append(
        f'<text class="table-header" x="{width - 36}" y="{row_y + 8}" '
        f'text-anchor="end">{sum(counts):,}</text>\n'
    )
    parts.append("</svg>\n")
    return "".join(parts)


def _grouped_density_svg(
    title: str,
    subtitle: str,
    x_label: str,
    y_label: str,
    series: list[tuple[str, list[float], int]],
) -> str:
    # Render several smoothed density curves in one SVG.
    width, height = 980, 560
    left, right, top, bottom = 80, 205, 75, 80
    plot_width = width - left - right
    plot_height = height - top - bottom
    max_value = max((max(values) for _, values, _ in series if values), default=0.0)
    tick_count = 5
    colors = {
        "one2one": "#1f6feb",
        "one2many": "#0f766e",
        "many2one": "#a16207",
        "complex/many2many": "#7c3aed",
        "uncertain_loss": "#c2410c",
        "lost": "#64748b",
        "missing": "#111827",
    }
    fallback_colors = ["#be123c", "#0369a1", "#4d7c0f", "#9333ea"]

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

    for index, (label, values, _) in enumerate(series):
        if not values:
            continue
        denom = max(1, len(values) - 1)
        points = []
        for value_index, value in enumerate(values):
            x = left + (value_index / denom) * plot_width
            y = _scale_y(value, max_value, top, plot_height)
            points.append((x, y))
        polyline = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        color = colors.get(label, fallback_colors[index % len(fallback_colors)])
        parts.append(
            f'<polyline points="{polyline}" fill="none" stroke="{color}" '
            'stroke-width="2.4"/>\n'
        )

    legend_x = width - right + 28
    legend_y = top + 6
    for index, (label, _, count) in enumerate(series):
        y = legend_y + index * 24
        color = colors.get(label, fallback_colors[index % len(fallback_colors)])
        parts.append(
            f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 22}" y2="{y}" '
            f'stroke="{color}" stroke-width="3"/>\n'
        )
        parts.append(
            f'<text class="tick" x="{legend_x + 30}" y="{y + 4}">'
            f'{html.escape(label)} ({count})</text>\n'
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
    # Render a smoothed density curve as a standalone SVG string.
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
