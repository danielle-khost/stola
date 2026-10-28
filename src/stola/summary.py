"""Concise per-reference-CNEE summaries from overlap statistics."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed

SUMMARY_COLUMNS = [
    "cnee_uid",
    "num_query_alignments",
    "overlap_fractions",
    "gap_bases",
    "sequence_identities",
    "sequence_mismatches",
]
# Decimal places for the overlap_fractions/sequence_identities lists, which
# hold one value per projection; overlap_stats.tsv keeps full precision.
SUMMARY_DECIMALS = 3


@dataclass
class CneeSummary:
    cnee_uid: str
    overlap_fractions: list[str] = field(default_factory=list)
    gap_bases: list[str] = field(default_factory=list)
    sequence_identities: list[str] = field(default_factory=list)
    sequence_mismatches: list[str] = field(default_factory=list)

    @property
    def num_query_alignments(self) -> int:
        # Count aligned query projections represented in this summary row.
        return len(self.overlap_fractions)

    def row(self) -> dict[str, str | int]:
        # Convert collected per-CNEE values into the output TSV row shape.
        return {
            "cnee_uid": self.cnee_uid,
            "num_query_alignments": self.num_query_alignments,
            "overlap_fractions": ",".join(self.overlap_fractions),
            "gap_bases": ",".join(self.gap_bases),
            "sequence_identities": ",".join(self.sequence_identities),
            "sequence_mismatches": ",".join(self.sequence_mismatches),
        }


def summarize_overlap_stats(
    stats_path: str | Path,
    out_path: str | Path,
    cnee_bed: str | Path | None = None,
) -> int:
    """Write a compact per-CNEE summary from an overlap stats TSV.

    Each reference CNEE's row lists the metrics of its candidate projections
    with at least one aligned CNEE base (its "query alignments"). With
    cnee_bed, every CNEE in it gets a row, even one with no candidates.
    """
    # Aggregate overlap, gap, and sequence metrics by reference CNEE ID.
    summaries: dict[str, CneeSummary] = {}
    if cnee_bed is not None:
        for interval in parse_bed(cnee_bed):
            summaries.setdefault(interval.cnee_uid, CneeSummary(interval.cnee_uid))

    with Path(stats_path).open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        required = {
            "cnee_uid",
            "aligned_bases",
            "covered_fraction",
            "gap_bases_within_chain_span",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            cnee_uid = row["cnee_uid"]
            summary = summaries.setdefault(cnee_uid, CneeSummary(cnee_uid))
            if int(row["aligned_bases"]) <= 0:
                continue
            summary.overlap_fractions.append(_rounded(row["covered_fraction"]))
            summary.gap_bases.append(row["gap_bases_within_chain_span"])
            seq_identity = row.get("seq_identity", ".")
            if seq_identity != ".":
                summary.sequence_identities.append(_rounded(seq_identity))
            if "seq_mismatches" in row and seq_identity != ".":
                summary.sequence_mismatches.append(row["seq_mismatches"])

    with Path(out_path).open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for summary in summaries.values():
            writer.writerow(summary.row())
    return len(summaries)


def _rounded(value: str) -> str:
    # Format one overlap-stats float to SUMMARY_DECIMALS places.
    return f"{float(value):.{SUMMARY_DECIMALS}f}"
