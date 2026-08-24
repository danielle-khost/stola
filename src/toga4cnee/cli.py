"""Command line interface for TOGA4CNEE."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .orthology import (
    DEFAULT_MIN_ALIGN_SCORE,
    DEFAULT_MIN_GAP_LENGTH,
    DEFAULT_MIN_NEIGHBOR_FRACTION,
    DEFAULT_QUERY_ID_PREFIX,
    DEFAULT_VERIFICATION_NUM_NEIGHBORS,
)
from .plots import write_qc_plots
from .summary import summarize_overlap_stats
from .synteny import DEFAULT_FLANK_SIZES
from .workflow import run_workflow


def _run(args: argparse.Namespace) -> None:
    # Launch the full one-command workflow; every output lands under args.outdir.
    result = run_workflow(
        args.chain_file,
        args.cnee_bed,
        args.ref_fasta,
        args.query_fasta,
        outdir=args.outdir,
        query_id_prefix=args.query_prefix,
        synteny_flank_sizes=args.synteny_flank_sizes,
        synteny_max_one2one_controls=args.synteny_max_one2one_controls,
        min_sequence_identity=args.min_sequence_identity,
        max_fraction_gap=args.max_fraction_gap,
        synteny_check=not args.no_synteny_check,
        min_align_score=args.min_align_score,
        min_neighbor_fraction=args.min_neighbor_fraction,
        verification_num_neighbors=args.verification_num_neighbors,
        min_gap_length=args.min_gap_length,
    )
    print(f"Indexed {result.chain_count} chains")
    print(f"Wrote {result.stats_count} all overlap-stat rows -> {result.stats_path}")
    print(f"Wrote {result.query_bed_count} query CNEE BED rows -> {result.query_bed_path}")
    print(
        f"Wrote {result.orthology_result.orthology_rows} query orthology rows -> "
        f"{result.orthology_path}"
    )
    print(f"Wrote {result.summary_count} all CNEE summary rows -> {result.summary_path}")
    print(
        f"Wrote {result.synteny_result.rows} synteny verification rows -> "
        f"{result.synteny_path}"
    )
    print(
        f"Wrote {result.synteny_neighbor_result.rows} synteny neighbor-check rows -> "
        f"{result.synteny_neighbors_path}"
    )
    print(
        f"Wrote {result.orthology_result.culled_projection_rows} pruned-projection rows -> "
        f"{result.culled_projections_path}"
    )
    print(
        f"Wrote {result.orthology_result.ref_cnee_rows} ref orthology rows -> "
        f"{result.ref_cnees_path}"
    )
    print(f"All outputs and intermediates written under {result.outdir}")


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
    run_parser.add_argument("ref_fasta", type=Path, help="Reference genome, FASTA or .2bit")
    run_parser.add_argument("query_fasta", type=Path, help="Query genome, FASTA or .2bit")
    run_parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=".",
        help=(
            "Directory for all outputs and intermediates (created if missing) "
            "[default: current directory]. See output_files.md for the fixed "
            "filenames written underneath it."
        ),
    )
    run_parser.add_argument(
        "--query-prefix",
        type=str,
        default=DEFAULT_QUERY_ID_PREFIX,
        help=(
            "Prefix for every graph-assigned query CNEE id, e.g. 'q_ce1234' "
            f"[default: {DEFAULT_QUERY_ID_PREFIX!r}]"
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
        help="Minimum sequence identity for a solid orthology edge [default: 0.7]",
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
        "--min-align-score",
        type=float,
        default=DEFAULT_MIN_ALIGN_SCORE,
        help=(
            "Minimum align_score (this projection's sequence identity relative "
            f"to the best in its connected component) to survive synteny "
            f"verification [default: {DEFAULT_MIN_ALIGN_SCORE}]"
        ),
    )
    run_parser.add_argument(
        "--min-neighbor-fraction",
        type=float,
        default=DEFAULT_MIN_NEIGHBOR_FRACTION,
        help=(
            "Minimum fraction of trusted neighboring reference CNEEs whose own "
            "query ortholog corroborates this projection's synteny, to survive "
            f"verification [default: {DEFAULT_MIN_NEIGHBOR_FRACTION}]"
        ),
    )
    run_parser.add_argument(
        "--verification-num-neighbors",
        type=int,
        default=DEFAULT_VERIFICATION_NUM_NEIGHBORS,
        help=(
            "Number of trusted (solid one2one) reference-CNEE neighbors checked "
            f"per flank during synteny verification [default: {DEFAULT_VERIFICATION_NUM_NEIGHBORS}]"
        ),
    )
    run_parser.add_argument(
        "--min-gap-length",
        type=int,
        default=DEFAULT_MIN_GAP_LENGTH,
        help=(
            "Minimum N-run length in the query FASTA considered a candidate "
            "deliberate assembly-gap marker, when inferring contig boundaries "
            f"for lost-ref sub-categorization [default: {DEFAULT_MIN_GAP_LENGTH}]"
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
