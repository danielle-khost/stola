"""Command line interface for TOGA4CNEE."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .plots import write_qc_plots
from .summary import summarize_overlap_stats
from .workflow import run_workflow


def _run(args: argparse.Namespace) -> None:
    # Resolve default output paths and launch the full one-command workflow.
    out_summary = (
        args.out_summary
        if args.out_summary is not None
        else args.out_stats.with_suffix(".summary.tsv")
    )
    out_alignment_histogram = (
        args.out_alignment_histogram
        if args.out_alignment_histogram is not None
        else args.out_stats.with_suffix(".num_query_alignments.svg")
    )
    out_overlap_density = (
        args.out_overlap_density
        if args.out_overlap_density is not None
        else args.out_stats.with_suffix(".overlap_fraction_density.svg")
    )
    out_sequence_identity_density = (
        args.out_sequence_identity_density
        if args.out_sequence_identity_density is not None
        else args.out_stats.with_suffix(".sequence_identity_density.svg")
    )
    out_orthology = (
        args.out_orthology
        if args.out_orthology is not None
        else args.out_stats.with_suffix(".orthology.tsv")
    )
    result = run_workflow(
        args.chain_file,
        args.cnee_bed,
        args.out_bed,
        args.out_stats,
        out_summary=out_summary,
        out_alignment_histogram=out_alignment_histogram,
        out_overlap_density=out_overlap_density,
        out_sequence_identity_density=out_sequence_identity_density,
        out_orthology=out_orthology,
        ref_fasta=args.ref_fasta,
        query_fasta=args.query_fasta,
        flank_size=args.flank_size,
        work_dir=args.work_dir,
    )
    print(f"Indexed {result.chain_count} chains")
    print(f"Wrote {result.stats_count} overlap-stat rows -> {args.out_stats}")
    print(f"Wrote {result.query_bed_count} query CNEE BED rows -> {args.out_bed}")
    print(
        f"Wrote {result.orthology_result.orthology_rows} orthology rows -> "
        f"{out_orthology}"
    )
    print(f"Wrote {result.summary_count} CNEE summary rows -> {out_summary}")
    print(f"Wrote alignment histogram -> {out_alignment_histogram}")
    print(f"Wrote overlap density plot -> {out_overlap_density}")
    print(f"Wrote sequence identity density plot -> {out_sequence_identity_density}")
    if args.work_dir is not None:
        print(f"Kept intermediates in {args.work_dir}")


def _summarize(args: argparse.Namespace) -> None:
    # Rebuild only the compact per-reference-CNEE summary from existing stats.
    count = summarize_overlap_stats(
        args.overlap_stats_tsv,
        args.out,
        cnee_bed=args.cnee_bed,
        aligned_only=not args.include_gap_only,
    )
    print(f"Wrote {count} CNEE summary rows -> {args.out}")


def _plot(args: argparse.Namespace) -> None:
    # Rebuild only the SVG QC plots from existing summary and stats files.
    result = write_qc_plots(
        args.summary_tsv,
        args.overlap_stats_tsv,
        args.out_alignment_histogram,
        args.out_overlap_density,
        args.out_sequence_identity_density,
    )
    print(f"Wrote alignment histogram for {result.histogram_n} CNEEs -> {result.histogram_path}")
    print(
        "Wrote overlap density for "
        f"{result.overlap_density_n} aligned candidate projections -> "
        f"{result.overlap_density_path}"
    )
    print(
        "Wrote sequence identity density for "
        f"{result.sequence_identity_density_n} sequence-compared projections -> "
        f"{result.sequence_identity_density_path}"
    )


def build_parser() -> argparse.ArgumentParser:
    # Define the run, summarize, and plot command-line interfaces.
    parser = argparse.ArgumentParser(
        prog="toga4cnee",
        description=(
            "Lift reference CNEEs through a chain file and report raw overlap "
            "statistics for every candidate chain/CNEE pair."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="Run the full chain/CNEE workflow",
        description=parser.description,
    )
    run_parser.add_argument("chain_file", type=Path, help="Input .chain or .chain.gz file")
    run_parser.add_argument("cnee_bed", type=Path, help="Reference CNEE BED3+ file")
    run_parser.add_argument(
        "--out-bed", required=True, type=Path, help="Output graph-assigned query CNEE BED"
    )
    run_parser.add_argument(
        "--out-stats", required=True, type=Path, help="Output overlap stats TSV"
    )
    run_parser.add_argument(
        "--out-summary",
        type=Path,
        default=None,
        help="Concise per-reference-CNEE summary TSV [default: OUT_STATS with .summary.tsv suffix]",
    )
    run_parser.add_argument(
        "--out-alignment-histogram",
        type=Path,
        default=None,
        help="SVG histogram of query alignment counts [default: OUT_STATS with .num_query_alignments.svg suffix]",
    )
    run_parser.add_argument(
        "--out-overlap-density",
        type=Path,
        default=None,
        help="SVG density plot of overlap fractions [default: OUT_STATS with .overlap_fraction_density.svg suffix]",
    )
    run_parser.add_argument(
        "--out-sequence-identity-density",
        type=Path,
        default=None,
        help="SVG density plot of sequence identity [default: OUT_STATS with .sequence_identity_density.svg suffix]",
    )
    run_parser.add_argument(
        "--out-orthology",
        type=Path,
        default=None,
        help="Graph-based orthology classification TSV [default: OUT_STATS with .orthology.tsv suffix]",
    )
    run_parser.add_argument("--flank-size", type=int, default=1000)
    run_parser.add_argument(
        "--ref-fasta",
        type=Path,
        default=None,
        help="Reference genome FASTA for sequence identity/mismatch metrics",
    )
    run_parser.add_argument(
        "--query-fasta",
        type=Path,
        default=None,
        help="Query genome FASTA for sequence identity/mismatch metrics",
    )
    run_parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory for intermediate index/candidate files. If omitted, "
            "intermediates are created in a temporary directory."
        ),
    )
    run_parser.set_defaults(func=_run)

    summary_parser = subparsers.add_parser(
        "summarize",
        help="Create the concise CNEE summary from an overlap stats TSV",
    )
    summary_parser.add_argument("overlap_stats_tsv", type=Path)
    summary_parser.add_argument("--out", required=True, type=Path)
    summary_parser.add_argument(
        "--cnee-bed",
        type=Path,
        default=None,
        help="Reference CNEE BED3+ file; when supplied, include CNEEs with zero alignments",
    )
    summary_parser.add_argument(
        "--include-gap-only",
        action="store_true",
        help="Include fully gapped spanning candidates as query alignments",
    )
    summary_parser.set_defaults(func=_summarize)

    plot_parser = subparsers.add_parser(
        "plot",
        help="Create QC SVG plots from existing summary and stats TSVs",
    )
    plot_parser.add_argument("summary_tsv", type=Path)
    plot_parser.add_argument("overlap_stats_tsv", type=Path)
    plot_parser.add_argument("--out-alignment-histogram", required=True, type=Path)
    plot_parser.add_argument("--out-overlap-density", required=True, type=Path)
    plot_parser.add_argument("--out-sequence-identity-density", required=True, type=Path)
    plot_parser.set_defaults(func=_plot)

    return parser


def main(argv: list[str] | None = None) -> None:
    # Parse CLI arguments, supporting bare run arguments as shorthand for run.
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] not in {"run", "summarize", "plot", "-h", "--help"}:
        argv = ["run", *argv]
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
