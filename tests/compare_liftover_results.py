#!/usr/bin/env python
"""Compare reference CNEE recovery between TOGA4CNEE and UCSC liftOver.

This is a standalone test/diagnostic script, not part of the main package CLI.

Example:
    conda activate biopython
    python tests/compare_liftover_results.py \
        --toga-bed /private/tmp/toga4cnee_pigeon_chicken_orthology/query_cnees.bed \
        --toga-orthology /private/tmp/toga4cnee_pigeon_chicken_orthology/overlap_stats.orthology.tsv \
        --out /private/tmp/pigeon_liftover_ref_comparison.tsv
"""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

import pandas as pd


DEFAULT_UCSC_BED = Path(__file__).resolve().parent / "HLcolLivi3.ucscLiftover.bed"
DEFAULT_UCSC_UNMAPPED = Path(__file__).resolve().parent / "unmapped"
DEFAULT_REF_BED = Path(__file__).resolve().parent / "gGal7_cnee_final_data_filtered.fixed.bed"


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


def read_ref_ids(path: str | Path) -> set[str]:
    # Read the source reference BED and return its reference CNEE IDs.
    return set(read_bed(path, "reference")["name"])


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


def collect_ref_ids(df: pd.DataFrame) -> set[str]:
    # Collapse BED rows into the set of source reference CNEEs they represent.
    ref_ids: set[str] = set()
    for row in df.itertuples(index=False):
        ref_ids.update(row.id_set)
    ref_ids.discard("")
    ref_ids.discard(".")
    return ref_ids


def write_venn_svg(
    out_path: str | Path,
    species: str,
    left_only: int,
    shared: int,
    right_only: int,
) -> None:
    # Draw a small dependency-free Venn diagram summarizing reference-ID overlap.
    width, height = 760, 440
    left_cx, right_cx, cy, radius = 310, 450, 210, 135
    title = html.escape(f"{species} reference CNEE recovery")
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text{{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}}
.title{{font-size:24px;font-weight:700}}
.label{{font-size:17px;font-weight:700}}
.count{{font-size:28px;font-weight:700}}
.note{{font-size:13px;fill:#52606d}}
</style>
<rect width="100%" height="100%" fill="#ffffff"/>
<text class="title" x="40" y="44">{title}</text>
<circle cx="{left_cx}" cy="{cy}" r="{radius}" fill="#3b82f6" fill-opacity="0.32" stroke="#2563eb" stroke-width="2"/>
<circle cx="{right_cx}" cy="{cy}" r="{radius}" fill="#f59e0b" fill-opacity="0.34" stroke="#d97706" stroke-width="2"/>
<text class="label" x="{left_cx - 75}" y="{cy - radius - 22}">TOGA4CNEE</text>
<text class="label" x="{right_cx - 52}" y="{cy - radius - 22}">UCSC liftOver</text>
<text class="count" x="{left_cx - 82}" y="{cy + 8}" text-anchor="middle">{left_only}</text>
<text class="count" x="{(left_cx + right_cx) / 2:.0f}" y="{cy + 8}" text-anchor="middle">{shared}</text>
<text class="count" x="{right_cx + 82}" y="{cy + 8}" text-anchor="middle">{right_only}</text>
<text class="note" x="{left_cx - 82}" y="{cy + 36}" text-anchor="middle">unique</text>
<text class="note" x="{(left_cx + right_cx) / 2:.0f}" y="{cy + 36}" text-anchor="middle">shared</text>
<text class="note" x="{right_cx + 82}" y="{cy + 36}" text-anchor="middle">unique</text>
</svg>
"""
    Path(out_path).write_text(svg, encoding="utf-8")


def compare_outputs(
    toga_bed: str | Path,
    ucsc_bed: str | Path,
    ref_bed: str | Path,
    toga_orthology: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    # Build a reference-ID recovery summary for both methods.
    ref_ids = read_ref_ids(ref_bed)
    orthology_map = read_orthology_map(toga_orthology)
    toga = add_comparison_keys(read_bed(toga_bed, "toga4cnee"), orthology_map)
    ucsc = add_comparison_keys(read_bed(ucsc_bed, "ucsc_liftover"))

    toga_ref_ids = collect_ref_ids(toga).intersection(ref_ids)
    ucsc_ref_ids = collect_ref_ids(ucsc).intersection(ref_ids)
    shared = toga_ref_ids.intersection(ucsc_ref_ids)
    toga_unique = toga_ref_ids.difference(ucsc_ref_ids)
    ucsc_unique = ucsc_ref_ids.difference(toga_ref_ids)

    rows = [
        {
            "category": "total ref annotations",
            "toga4cnee": len(ref_ids),
            "ucsc_liftover": len(ref_ids),
        },
        {
            "category": "total query annotations",
            "toga4cnee": len(toga),
            "ucsc_liftover": len(ucsc),
        },
        {
            "category": "total shared annotations",
            "toga4cnee": len(shared),
            "ucsc_liftover": len(shared),
        },
        {
            "category": "total unique annotations",
            "toga4cnee": len(toga_unique),
            "ucsc_liftover": len(ucsc_unique),
        },
    ]
    venn_counts = {
        "toga4cnee_unique": len(toga_unique),
        "shared": len(shared),
        "ucsc_liftover_unique": len(ucsc_unique),
    }
    return pd.DataFrame(rows, columns=["category", "toga4cnee", "ucsc_liftover"]), venn_counts


def build_parser() -> argparse.ArgumentParser:
    # Define file paths for standalone comparison runs.
    parser = argparse.ArgumentParser(
        description="Compare TOGA4CNEE query BED output to UCSC liftOver results."
    )
    parser.add_argument("--toga-bed", required=True, type=Path)
    parser.add_argument("--ref-bed", default=DEFAULT_REF_BED, type=Path)
    parser.add_argument("--toga-orthology", default=None, type=Path)
    parser.add_argument("--ucsc-bed", default=DEFAULT_UCSC_BED, type=Path)
    parser.add_argument("--out", default=None, type=Path)
    parser.add_argument("--out-venn", default=None, type=Path)
    parser.add_argument(
        "--toga-summary",
        default=None,
        type=Path,
        help="Accepted for compatibility with older comparison commands; unused.",
    )
    parser.add_argument(
        "--ucsc-unmapped",
        default=DEFAULT_UCSC_UNMAPPED,
        type=Path,
        help="Accepted for compatibility with older comparison commands; unused.",
    )
    return parser


def main() -> None:
    # Run the comparison and print or write the summary table.
    args = build_parser().parse_args()
    species = Path(args.toga_bed).resolve().parent.name
    summary, venn_counts = compare_outputs(
        toga_bed=args.toga_bed,
        toga_orthology=args.toga_orthology,
        ucsc_bed=args.ucsc_bed,
        ref_bed=args.ref_bed,
    )
    out_venn = args.out_venn
    if out_venn is None and args.out is not None:
        out_venn = args.out.with_suffix(".venn.svg")
    if args.out is None:
        print(summary.to_csv(sep="\t", index=False), end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote comparison summary -> {args.out}")
    if out_venn is not None:
        write_venn_svg(
            out_venn,
            species=species,
            left_only=venn_counts["toga4cnee_unique"],
            shared=venn_counts["shared"],
            right_only=venn_counts["ucsc_liftover_unique"],
        )
        print(f"Wrote Venn diagram -> {out_venn}")


if __name__ == "__main__":
    main()
