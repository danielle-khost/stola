"""Command line interface for TOGA4CNEE."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .plots import write_qc_plots
from .summary import summarize_overlap_stats
from .synteny import DEFAULT_FLANK_SIZES
from .workflow import run_workflow


def _run(args: argparse.Namespace) -> None:
    # Resolve default output paths and launch the full one-command workflow.
    out_summary = (
        args.out_summary
        if args.out_summary is not None
        else args.out_stats.with_suffix(".summary.tsv")
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
        out_orthology=out_orthology,
        out_strict_stats=args.out_strict_stats,
        out_strict_summary=args.out_strict_summary,
        out_synteny=args.out_synteny,
        out_synteny_neighbors=args.out_synteny_neighbors,
        ref_fasta=args.ref_fasta,
        query_fasta=args.query_fasta,
        synteny_flank_sizes=args.synteny_flank_sizes,
        synteny_max_one2one_controls=args.synteny_max_one2one_controls,
        min_sequence_identity=args.min_sequence_identity,
        max_fraction_gap=args.max_fraction_gap,
        synteny_check=not args.no_synteny_check,
        work_dir=args.work_dir,
    )
    print(f"Indexed {result.chain_count} chains")
    print(f"Wrote {result.stats_count} all overlap-stat rows -> {args.out_stats}")
    print(
        f"Wrote {result.strict_stats_count} strict overlap-stat rows -> "
        f"{result.strict_stats_path}"
    )
    print(f"Wrote {result.query_bed_count} query CNEE BED rows -> {args.out_bed}")
    print(
        f"Wrote {result.orthology_result.orthology_rows} orthology rows -> "
        f"{out_orthology}"
    )
    print(f"Wrote {result.summary_count} all CNEE summary rows -> {out_summary}")
    print(
        f"Wrote {result.strict_summary_count} strict CNEE summary rows -> "
        f"{result.strict_summary_path}"
    )
    print(
        f"Wrote {result.synteny_result.rows} synteny verification rows -> "
        f"{result.synteny_path}"
    )
    print(
        f"Wrote {result.synteny_neighbor_result.rows} synteny neighbor-check rows -> "
        f"{result.synteny_neighbors_path}"
    )
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
        "--out-orthology",
        type=Path,
        default=None,
        help="Graph-based orthology classification TSV [default: OUT_STATS with .orthology.tsv suffix]",
    )
    run_parser.add_argument(
        "--out-strict-stats",
        type=Path,
        default=None,
        help="Strict solid-only overlap stats TSV [default: OUT_STATS with .strict before extension]",
    )
    run_parser.add_argument(
        "--out-strict-summary",
        type=Path,
        default=None,
        help="Strict per-reference-CNEE summary TSV [default: OUT_STATS with .strict.summary before extension]",
    )
    run_parser.add_argument(
        "--out-synteny",
        type=Path,
        default=None,
        help="Flanking-synteny verification TSV [default: OUT_STATS with .synteny before extension]",
    )
    run_parser.add_argument(
        "--out-synteny-neighbors",
        type=Path,
        default=None,
        help=(
            "Left/right flanking-CNEE synteny check TSV [default: OUT_STATS with "
            ".synteny_neighbors before extension]"
        ),
    )
    run_parser.add_argument(
        "--synteny-flank-sizes",
        type=int,
        nargs="+",
        default=list(DEFAULT_FLANK_SIZES),
        help=(
            "One or more reference flank sizes (bases upstream/downstream) for "
            "flanking-synteny verification; the output TSV gets one row per "
            f"candidate per size [default: {list(DEFAULT_FLANK_SIZES)}]"
        ),
    )
    run_parser.add_argument(
        "--synteny-max-one2one-controls",
        type=int,
        default=50000,
        help="Maximum one2one positive-control projections in synteny output [default: 50000]",
    )
    run_parser.add_argument(
        "--min-sequence-identity",
        type=float,
        default=0.7,
        help=(
            "Minimum sequence identity for a solid orthology edge when FASTAs "
            "are supplied [default: 0.7]"
        ),
    )
    run_parser.add_argument(
        "--max-fraction-gap",
        type=float,
        default=0.1,
        help="Maximum gap fraction for a solid orthology edge [default: 0.1]",
    )
    run_parser.add_argument(
        "--no-synteny-check",
        action="store_true",
        help=(
            "Skip synteny verification/pruning of ambiguous orthology edges "
            "(a non-solid single ref/query edge, or any edge in a one2many/"
            "many2one/complex component). Categories are still assigned from "
            "the unpruned graph, just without the extra quality gate."
        ),
    )
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
