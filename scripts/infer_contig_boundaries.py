#!/usr/bin/env python
"""CLI wrapper for toga4cnee.assembly_gaps -- report a genome's inferred gap convention.

See toga4cnee/assembly_gaps.py for the algorithm. This script scans a FASTA,
writes every N-run found to a TSV, and prints the length distribution plus
the confidence verdict.

Example:
    PYTHONPATH=src python scripts/infer_contig_boundaries.py \\
        --fasta tests/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna \\
        --out tests/test_runs/colLivi/query_n_runs.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from toga4cnee.assembly_gaps import (
    DEFAULT_MIN_GAP_LENGTH,
    MIN_RUNS_FOR_INFERENCE,
    infer_gap_convention,
    iter_n_runs,
)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="Per-N-run TSV output path")
    parser.add_argument(
        "--min-gap-length",
        type=int,
        default=DEFAULT_MIN_GAP_LENGTH,
        help=(
            "Minimum N-run length considered a candidate deliberate gap marker "
            f"(rather than an isolated ambiguous basecall) [default: {DEFAULT_MIN_GAP_LENGTH}]"
        ),
    )
    args = parser.parse_args(argv)

    print(f"Scanning {args.fasta} for N-runs...", file=sys.stderr)
    runs = list(iter_n_runs(args.fasta))
    print(f"  found {len(runs)} N-run(s) across the genome", file=sys.stderr)

    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "length"])
        for run in runs:
            writer.writerow([run.chrom, run.start, run.end, run.length])
    print(f"Wrote {len(runs)} rows -> {args.out}")

    convention = infer_gap_convention(runs, args.min_gap_length)

    print("\nLength distribution (top 10 by count, all run lengths):")
    label_width = max((len(str(length)) for length in convention.length_counts), default=1)
    for length, count in convention.length_counts.most_common(10):
        print(f"  length={str(length).rjust(label_width)}  count={count}")
    if len(convention.length_counts) > 10:
        print(f"  ... {len(convention.length_counts) - 10} more distinct length(s)")

    total_filtered = sum(
        c for length, c in convention.length_counts.items() if length >= args.min_gap_length
    )
    print(
        f"\nInference (restricted to runs >= --min-gap-length={args.min_gap_length}, "
        f"n={total_filtered}):"
    )
    if convention.confidence == "high":
        (length,) = convention.boundary_lengths
        count = convention.length_counts[length]
        print(
            f"  HIGH confidence: {length}bp N-runs look like this assembly's fixed "
            f"gap-of-unknown-size marker ({count} occurrences, {count / total_filtered:.1%} "
            "of qualifying runs)."
        )
    elif convention.confidence == "medium":
        print(
            "  MEDIUM confidence: no single dominant length, but this small set covers "
            "most qualifying runs -- plausibly a few distinct gap-type conventions:"
        )
        for length in sorted(convention.boundary_lengths, key=lambda l: -convention.length_counts[l]):
            count = convention.length_counts[length]
            print(f"    {length}bp: {count} occurrences ({count / total_filtered:.1%})")
    elif total_filtered < MIN_RUNS_FOR_INFERENCE:
        print(
            f"  Could not confidently infer a gap-length convention -- only "
            f"{total_filtered} qualifying N-run(s) found, too few to draw a "
            "conclusion either way. This assembly may be effectively gapless, or "
            "use isolated ambiguous bases rather than deliberate gap padding."
        )
    else:
        print(
            "  Could not confidently infer a fixed gap-length convention -- run "
            "lengths are too spread out (no length or small set of lengths "
            "dominates). This assembly likely records estimated gap sizes rather "
            "than a fixed placeholder, or mixes conventions. See the full "
            "distribution above/in --out to judge by eye."
        )


if __name__ == "__main__":
    main()
