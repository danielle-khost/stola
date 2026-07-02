#!/usr/bin/env python3
"""Generate TOGA4CNEE diagnostic figures from completed test outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from toga4cnee.plots import write_orthology_diagnostic_plots, write_qc_plots


def build_parser() -> argparse.ArgumentParser:
    # Define a test-oriented plotting interface separate from the main workflow.
    parser = argparse.ArgumentParser(
        description=(
            "Generate SVG diagnostic plots from completed TOGA4CNEE output files."
        )
    )
    parser.add_argument("--summary", required=True, type=Path, help="Summary TSV")
    parser.add_argument("--stats", required=True, type=Path, help="Overlap stats TSV")
    parser.add_argument(
        "--orthology", required=True, type=Path, help="Orthology classification TSV"
    )
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix [default: overlap stats stem]",
    )
    parser.add_argument(
        "--min-sequence-identity",
        type=float,
        default=0.7,
        help="Minimum sequence identity used to define solid projections",
    )
    parser.add_argument(
        "--max-fraction-gap",
        type=float,
        default=0.9,
        help="Maximum gap fraction used to define solid projections",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    # Write the previous QC plots plus orthology-category diagnostic plots.
    args = build_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix if args.prefix is not None else args.stats.stem

    alignment_histogram = args.out_dir / f"{prefix}.num_query_alignments.svg"
    overlap_density = args.out_dir / f"{prefix}.overlap_fraction_density.svg"
    sequence_identity_density = args.out_dir / f"{prefix}.sequence_identity_density.svg"
    category_counts = args.out_dir / f"{prefix}.orthology_category_counts.svg"
    sequence_by_category = (
        args.out_dir / f"{prefix}.orthology_sequence_identity_density.svg"
    )
    gap_by_category = args.out_dir / f"{prefix}.orthology_gap_fraction_density.svg"

    qc_result = write_qc_plots(
        args.summary,
        args.stats,
        alignment_histogram,
        overlap_density,
        sequence_identity_density,
    )
    diagnostic_result = write_orthology_diagnostic_plots(
        args.stats,
        args.orthology,
        category_counts,
        sequence_by_category,
        gap_by_category,
        min_sequence_identity=args.min_sequence_identity,
        max_fraction_gap=args.max_fraction_gap,
    )

    print(f"Wrote alignment histogram -> {qc_result.histogram_path}")
    print(f"Wrote overlap density -> {qc_result.overlap_density_path}")
    print(f"Wrote sequence identity density -> {qc_result.sequence_identity_density_path}")
    print(f"Wrote orthology category counts -> {diagnostic_result.category_counts_path}")
    print(
        "Wrote sequence identity by category -> "
        f"{diagnostic_result.sequence_identity_by_category_path}"
    )
    print(
        "Wrote gap fraction by category -> "
        f"{diagnostic_result.gap_fraction_by_category_path}"
    )
    print("Orthology category counts:")
    for category, count in sorted(diagnostic_result.category_counts.items()):
        print(f"  {category}\t{count}")


if __name__ == "__main__":
    main()
