#!/usr/bin/env python
"""Compare TOGA4CNEE query CNEEs against UCSC liftOver outputs.

This is a standalone test/diagnostic script, not part of the main package CLI.

Example:
    conda activate biopython
    python tests/compare_liftover_results.py \
        --toga-bed /private/tmp/toga4cnee_pigeon_chicken_orthology/query_cnees.bed \
        --toga-summary /private/tmp/toga4cnee_pigeon_chicken_orthology/overlap_stats.summary.tsv \
        --toga-orthology /private/tmp/toga4cnee_pigeon_chicken_orthology/overlap_stats.orthology.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_UCSC_BED = Path(__file__).resolve().parent / "HLcolLivi3.ucscLiftover.bed"
DEFAULT_UCSC_UNMAPPED = Path(__file__).resolve().parent / "unmapped"


def read_bed(path: str | Path, source: str) -> pd.DataFrame:
    # Load BED3+ rows, keeping only coordinates and the column-4 annotation ID.
    rows = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number} has fewer than 3 BED columns")
            rows.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "name": fields[3] if len(fields) > 3 and fields[3] else ".",
                    "source": source,
                }
            )
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "source"])


def read_orthology_map(path: str | Path | None) -> dict[str, set[str]]:
    # Map graph-assigned query IDs back to partner reference CNEE IDs.
    if path is None:
        return {}
    orthology_path = Path(path)
    if not orthology_path.exists():
        return {}
    df = pd.read_csv(orthology_path, sep="\t", dtype=str).fillna("")
    required = {"query_cnee_id", "partner_ref_cnees"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"{orthology_path} is missing columns: {', '.join(sorted(missing))}"
        )
    mapping: dict[str, set[str]] = {}
    for row in df.itertuples(index=False):
        query_id = getattr(row, "query_cnee_id")
        if query_id == "NA":
            continue
        partners = {
            value
            for value in getattr(row, "partner_ref_cnees").split(",")
            if value
        }
        mapping.setdefault(query_id, set()).update(partners)
    return mapping


def infer_toga_ref_ids(query_id: str, orthology_map: dict[str, set[str]]) -> set[str]:
    # Prefer orthology TSV mappings, then fall back to simple q_ref or q_ref_N IDs.
    if query_id in orthology_map:
        return orthology_map[query_id]
    if not query_id.startswith("q_"):
        return {query_id}
    stripped = query_id[2:]
    match = re.match(r"(.+)_\d+$", stripped)
    if match and not stripped.startswith(("complex_", "many2one_")):
        return {match.group(1)}
    if not stripped.startswith(("complex_", "many2one_")):
        return {stripped}
    return set()


def add_comparison_keys(
    df: pd.DataFrame, orthology_map: dict[str, set[str]] | None = None
) -> pd.DataFrame:
    # Add exact locus keys and comparable reference CNEE ID sets.
    output = df.copy()
    output["locus_key"] = list(zip(output["chrom"], output["start"], output["end"]))
    if orthology_map is None:
        output["id_set"] = output["name"].map(lambda value: {str(value)})
    else:
        output["id_set"] = output["name"].map(
            lambda value: infer_toga_ref_ids(str(value), orthology_map)
        )
    return output


def summarize_toga_missing(path: str | Path) -> int:
    # Count TOGA4CNEE reference CNEEs with no aligned query BED row.
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "projection_status" in df.columns:
        return int(df["projection_status"].isin(["lost", "missing"]).sum())
    if "num_query_alignments" not in df.columns:
        raise ValueError(
            f"{path} needs projection_status or num_query_alignments to count missing rows"
        )
    return int((df["num_query_alignments"].astype(int) == 0).sum())


def read_unmapped_ids(path: str | Path) -> set[str]:
    # Read UCSC liftOver unmapped BED-like records, skipping reason/comment lines.
    ids: set[str] = set()
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 4:
                ids.add(fields[3])
    return ids


def build_locus_index(df: pd.DataFrame) -> dict[tuple[str, int, int], set[str]]:
    # Index rows by exact output locus, merging all comparable IDs at that locus.
    grouped: dict[tuple[str, int, int], set[str]] = {}
    for row in df.itertuples(index=False):
        grouped.setdefault(row.locus_key, set()).update(row.id_set)
    return grouped


def classify_rows(
    df: pd.DataFrame, other_by_locus: dict[tuple[str, int, int], set[str]]
) -> dict[str, int]:
    # Classify each BED row once relative to the other method's exact loci.
    counts = {"perfect": 0, "different_ids": 0, "unique": 0}
    for row in df.itertuples(index=False):
        other_ids = other_by_locus.get(row.locus_key)
        if other_ids is None:
            counts["unique"] += 1
        elif row.id_set.intersection(other_ids):
            counts["perfect"] += 1
        else:
            counts["different_ids"] += 1
    return counts


def compare_outputs(
    toga_bed: str | Path,
    toga_summary: str | Path,
    ucsc_bed: str | Path,
    ucsc_unmapped: str | Path,
    toga_orthology: str | Path | None = None,
) -> pd.DataFrame:
    # Build the requested method-level comparison summary table.
    orthology_map = read_orthology_map(toga_orthology)
    toga = add_comparison_keys(read_bed(toga_bed, "toga4cnee"), orthology_map)
    ucsc = add_comparison_keys(read_bed(ucsc_bed, "ucsc_liftover"))

    toga_by_locus = build_locus_index(toga)
    ucsc_by_locus = build_locus_index(ucsc)
    toga_counts = classify_rows(toga, ucsc_by_locus)
    ucsc_counts = classify_rows(ucsc, toga_by_locus)

    rows = [
        {
            "category": "total annotations",
            "toga4cnee": len(toga),
            "ucsc_liftover": len(ucsc),
        },
        {
            "category": "total missing",
            "toga4cnee": summarize_toga_missing(toga_summary),
            "ucsc_liftover": len(read_unmapped_ids(ucsc_unmapped)),
        },
        {
            "category": "perfect match",
            "toga4cnee": toga_counts["perfect"],
            "ucsc_liftover": ucsc_counts["perfect"],
        },
        {
            "category": "different IDs",
            "toga4cnee": toga_counts["different_ids"],
            "ucsc_liftover": ucsc_counts["different_ids"],
        },
        {
            "category": "unique annotations",
            "toga4cnee": toga_counts["unique"],
            "ucsc_liftover": ucsc_counts["unique"],
        },
    ]
    return pd.DataFrame(rows, columns=["category", "toga4cnee", "ucsc_liftover"])


def build_parser() -> argparse.ArgumentParser:
    # Define file paths for standalone comparison runs.
    parser = argparse.ArgumentParser(
        description="Compare TOGA4CNEE query BED output to UCSC liftOver results."
    )
    parser.add_argument("--toga-bed", required=True, type=Path)
    parser.add_argument("--toga-summary", required=True, type=Path)
    parser.add_argument("--toga-orthology", default=None, type=Path)
    parser.add_argument("--ucsc-bed", default=DEFAULT_UCSC_BED, type=Path)
    parser.add_argument("--ucsc-unmapped", default=DEFAULT_UCSC_UNMAPPED, type=Path)
    parser.add_argument("--out", default=None, type=Path)
    return parser


def main() -> None:
    # Run the comparison and print or write the summary table.
    args = build_parser().parse_args()
    summary = compare_outputs(
        toga_bed=args.toga_bed,
        toga_summary=args.toga_summary,
        toga_orthology=args.toga_orthology,
        ucsc_bed=args.ucsc_bed,
        ucsc_unmapped=args.ucsc_unmapped,
    )
    if args.out is None:
        print(summary.to_csv(sep="\t", index=False), end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote comparison summary -> {args.out}")


if __name__ == "__main__":
    main()
