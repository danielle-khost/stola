"""Command line interface for TOGA4CNEE."""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

from .orthology import (
    DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    DEFAULT_FLANK_WINDOW,
    DEFAULT_MIN_GAP_LENGTH,
    DEFAULT_QUERY_ID_PREFIX,
    DEFAULT_RELATIVE_SCORE_THRESHOLD,
    DEFAULT_SOFT_DEPTH_LIMIT,
)
from .plots import write_aligned_block_depth_density, write_qc_plots
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
        filter_multi_projections=not args.no_filter_multi_projections,
        soft_depth_limit=args.soft_depth_limit,
        relative_score_threshold=args.relative_score_threshold,
        flank_window=args.flank_window,
        flank_any_coverage_threshold=args.flank_any_coverage_threshold,
        min_gap_length=args.min_gap_length,
    )
    print(f"Indexed {result.chain_count} chains")
    if result.named_ref_bed_path is not None:
        print(
            "Reference BED had no name column; wrote generated CNEE ids -> "
            f"{result.named_ref_bed_path}"
        )
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


def _chain_depth(args: argparse.Namespace) -> None:
    # Compute per-base aligned-block depth across both genomes of a chain
    # file and write the combined reference/query density SVG.
    from .chain_depth import compute_aligned_block_depth_histograms, compute_depth_stats

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
            "Prefix for every graph-assigned query CNEE id, e.g. 'q_ce1234'. "
            "A trailing '_' is added automatically if not already present, "
            "so 'aphWoo' and 'aphWoo_' both give ids like 'aphWoo_ce1234' "
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
        help=(
            "Minimum sequence identity for a 'solid' candidate in the synteny "
            "diagnostic reports and the overlap-stats summary (does not gate "
            "orthology categorization) [default: 0.7]"
        ),
    )
    run_parser.add_argument(
        "--max-fraction-gap",
        type=float,
        default=0.1,
        help=(
            "Maximum gap fraction for a 'solid' candidate in the synteny "
            "diagnostic reports and the overlap-stats summary (does not gate "
            "orthology categorization) [default: 0.1]"
        ),
    )
    run_parser.add_argument(
        "--no-filter-multi-projections",
        action="store_true",
        help=(
            "Skip relative-score filtering of multi-projection orthology edges "
            "(any edge in a one2many/many2one/complex component). A lone "
            "ref/query edge is never filtered either way. Categories are "
            "still assigned from the unfiltered graph, just without the "
            "extra quality gate."
        ),
    )
    run_parser.add_argument(
        "--soft-depth-limit",
        type=int,
        default=DEFAULT_SOFT_DEPTH_LIMIT,
        help=(
            "At each query locus, score only this many chains with aligned "
            "bases there, taking the best by chain score; a projection on any "
            "other chain is rejected as part of a repeat-driven pileup "
            f"[default: {DEFAULT_SOFT_DEPTH_LIMIT}]"
        ),
    )
    run_parser.add_argument(
        "--relative-score-threshold",
        type=float,
        default=DEFAULT_RELATIVE_SCORE_THRESHOLD,
        help=(
            "Keep a projection if its chain's local score (sequence identity "
            "over the locus x aligned fraction over the locus and its flanks) "
            "is above this fraction of the best chain's local score at that "
            f"locus [default: {DEFAULT_RELATIVE_SCORE_THRESHOLD}]"
        ),
    )
    run_parser.add_argument(
        "--flank-window",
        type=int,
        default=DEFAULT_FLANK_WINDOW,
        help=(
            "Query-side flank size (bases each side) included in each chain's "
            f"aligned fraction, and checked at single-chain loci [default: {DEFAULT_FLANK_WINDOW}]"
        ),
    )
    run_parser.add_argument(
        "--flank-any-coverage-threshold",
        type=float,
        default=DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
        help=(
            "At a locus where only one chain aligns: minimum fraction of the "
            "projection's flanks that must have aligned-block coverage from "
            f"ANY chain to accept it [default: {DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD}]"
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

    chain_depth_parser = subparsers.add_parser(
        "chain-depth",
        help="Per-base aligned-block depth density plot for a chain file's ref/query genomes",
    )
    chain_depth_parser.add_argument("chain_file", type=Path, help="Input .chain or .chain.gz file")
    chain_depth_parser.add_argument("--out", required=True, type=Path, help="Output SVG path")
    chain_depth_parser.add_argument(
        "--bin-max",
        type=int,
        default=100,
        help="Collapse every depth above this into one overflow bin [default: 100]",
    )
    chain_depth_parser.set_defaults(func=_chain_depth)

    return parser


def main(argv: list[str] | None = None) -> None:
    # Parse CLI arguments, supporting bare run arguments as shorthand for run.
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] not in {"run", "summarize", "plot", "chain-depth", "-h", "--help"}:
        argv = ["run", *argv]
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
