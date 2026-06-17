"""Concise per-reference-CNEE summaries from overlap statistics."""

from __future__ import annotations

import csv
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed

SUMMARY_COLUMNS = [
    "cnee_id",
    "num_query_alignments",
    "overlap_fractions",
    "gap_bases",
]


@dataclass
class CneeSummary:
    cnee_id: str
    overlap_fractions: list[str] = field(default_factory=list)
    gap_bases: list[str] = field(default_factory=list)

    @property
    def num_query_alignments(self) -> int:
        return len(self.overlap_fractions)

    def row(self) -> dict[str, str | int]:
        return {
            "cnee_id": self.cnee_id,
            "num_query_alignments": self.num_query_alignments,
            "overlap_fractions": ",".join(self.overlap_fractions),
            "gap_bases": ",".join(self.gap_bases),
        }


def summarize_overlap_stats(
    stats_path: str | Path,
    out_path: str | Path,
    cnee_bed: str | Path | None = None,
    aligned_only: bool = True,
) -> int:
    """Write a compact per-CNEE summary from an overlap stats TSV.

    By default, only rows with at least one aligned CNEE base are counted as
    query alignments. This matches the query BED output, which omits fully
    gapped spanning candidates.
    """
    summaries: OrderedDict[str, CneeSummary] = OrderedDict()
    if cnee_bed is not None:
        for interval in parse_bed(cnee_bed):
            summaries.setdefault(interval.cnee_id, CneeSummary(cnee_id=interval.cnee_id))

    with Path(stats_path).open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        required = {
            "cnee_id",
            "aligned_bases",
            "covered_fraction",
            "gap_bases_within_chain_span",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            if aligned_only and int(row["aligned_bases"]) <= 0:
                continue
            cnee_id = row["cnee_id"]
            summary = summaries.setdefault(cnee_id, CneeSummary(cnee_id=cnee_id))
            summary.overlap_fractions.append(row["covered_fraction"])
            summary.gap_bases.append(row["gap_bases_within_chain_span"])

    with Path(out_path).open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t")
        writer.writeheader()
        for summary in summaries.values():
            writer.writerow(summary.row())
    return len(summaries)
