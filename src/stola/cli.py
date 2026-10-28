"""Command line interface for STOLA."""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
from typing import Callable

from . import __version__
from .orthology import (
    DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    DEFAULT_FLANK_WINDOW,
    DEFAULT_MIN_GAP_LENGTH,
    DEFAULT_QUERY_ID_PREFIX,
    DEFAULT_RELATIVE_SCORE_THRESHOLD,
    DEFAULT_SOFT_DEPTH_LIMIT,
)
from .workflow import run_workflow


def bounded(
    kind: type,
    minimum: float | None = None,
    maximum: float | None = None,
    below_maximum: bool = False,
) -> Callable[[str], float]:
    """argparse type for a number parsed with kind (int or float) that must be
    >= minimum and <= maximum (or < maximum, with below_maximum)."""

    def parse(text: str):
        try:
            value = kind(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"invalid {kind.__name__} value: {text!r}") from None
        if not math.isfinite(value):
            raise argparse.ArgumentTypeError(f"must be a finite number (got {text})")
        too_low = minimum is not None and value < minimum
        too_high = maximum is not None and (value >= maximum if below_maximum else value > maximum)
        if too_low or too_high:
            if maximum is None:
                allowed = f">= {minimum}"
            elif below_maximum:
                allowed = f">= {minimum} and < {maximum}"
            else:
                allowed = f"between {minimum} and {maximum}"
            raise argparse.ArgumentTypeError(f"must be {allowed} (got {text})")
        return value

    return parse


def existing_file(text: str) -> Path:
    """argparse type for an input file that must already exist."""
    path = Path(text)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {text}")
    return path


def _print_progress(message: str) -> None:
    # Timestamp each line and flush it at once, so progress shows up as it
    # happens even when stdout is redirected to a log file.
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def _run(args: argparse.Namespace) -> None:
    # Launch the full one-command workflow; every output lands under args.outdir.
    started = time.monotonic()
    _print_progress(f"stola {__version__}: writing outputs under {args.outdir}")
    result = run_workflow(
        args.chain_file,
        args.cnee_bed,
        args.ref_fasta,
        args.query_fasta,
        outdir=args.outdir,
        query_id_prefix=args.query_prefix,
        filter_multi_projections=not args.no_filter_multi_projections,
        soft_depth_limit=args.soft_depth_limit,
        relative_score_threshold=args.relative_score_threshold,
        flank_window=args.flank_window,
        flank_any_coverage_threshold=args.flank_any_coverage_threshold,
        min_gap_length=args.min_gap_length,
        keep_intermediates=args.keep_intermediates,
        progress=_print_progress,
    )
    minutes, seconds = divmod(round(time.monotonic() - started), 60)
    _print_progress(
        f"All outputs and intermediates written under {result.outdir} ({minutes}m {seconds:02d}s)"
    )


def build_parser() -> argparse.ArgumentParser:
    # Define the single-command workflow interface.
    parser = argparse.ArgumentParser(
        prog="stola",
        description=(
            "Lift reference CNEEs through a chain file, classify their query "
            "orthology, and write summary tables and figures."
        ),
    )
    parser.add_argument("chain_file", type=existing_file, help="Input .chain or .chain.gz file")
    parser.add_argument("cnee_bed", type=existing_file, help="Reference CNEE BED3+ file")
    parser.add_argument("ref_fasta", type=existing_file, help="Reference genome, FASTA or .2bit")
    parser.add_argument("query_fasta", type=existing_file, help="Query genome, FASTA or .2bit")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=".",
        help=(
            "Directory for all outputs and intermediates (created if missing) "
            "[default: current directory]. See the README's Output section for "
            "the fixed filenames written underneath it."
        ),
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help=(
            "Keep the large intermediates under OUTDIR/work (the chain index, "
            "the decompressed copy of a .gz chain file, and the candidate "
            "chain/CNEE pairs) instead of deleting them after a successful "
            "run. Only needed for scripts/plot_chain_pileup.py."
        ),
    )
    parser.add_argument(
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
    parser.add_argument(
        "--no-filter-multi-projections",
        action="store_true",
        help=(
            "Skip the projection filtering step for multi-projection orthology "
            "edges (any edge in a one2many/many2one/complex component). A lone "
            "ref/query edge is never filtered either way. Categories are "
            "still assigned from the unfiltered graph, just without the "
            "extra quality gate."
        ),
    )
    parser.add_argument(
        "--soft-depth-limit",
        type=bounded(int, minimum=1),
        default=DEFAULT_SOFT_DEPTH_LIMIT,
        help=(
            "At each query locus, score only this many chains with aligned "
            "bases there, taking the best by chain score; a projection on any "
            "other chain is rejected as part of a repeat-driven pileup "
            f"[default: {DEFAULT_SOFT_DEPTH_LIMIT}]"
        ),
    )
    parser.add_argument(
        "--relative-score-threshold",
        type=bounded(float, minimum=0, maximum=1, below_maximum=True),
        default=DEFAULT_RELATIVE_SCORE_THRESHOLD,
        help=(
            "Keep a projection if its chain's local score (sequence identity "
            "over the locus x aligned fraction over the locus and its flanks) "
            "is above this fraction of the best chain's local score at that "
            f"locus [default: {DEFAULT_RELATIVE_SCORE_THRESHOLD}]"
        ),
    )
    parser.add_argument(
        "--flank-window",
        type=bounded(int, minimum=1),
        default=DEFAULT_FLANK_WINDOW,
        help=(
            "Query-side flank size (bases each side) included in each chain's "
            f"aligned fraction, and checked at single-chain loci [default: {DEFAULT_FLANK_WINDOW}]"
        ),
    )
    parser.add_argument(
        "--flank-any-coverage-threshold",
        type=bounded(float, minimum=0, maximum=1, below_maximum=True),
        default=DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
        help=(
            "At a locus where only one chain aligns: accept the projection if "
            "the fraction of its flanks with aligned-block coverage from ANY "
            f"chain is above this value [default: {DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD}]"
        ),
    )
    parser.add_argument(
        "--min-gap-length",
        type=bounded(int, minimum=1),
        default=DEFAULT_MIN_GAP_LENGTH,
        help=(
            "Minimum N-run length in the query genome considered a candidate "
            "deliberate assembly-gap marker, when inferring contig boundaries "
            f"for lost-ref sub-categorization [default: {DEFAULT_MIN_GAP_LENGTH}]"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    _run(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
