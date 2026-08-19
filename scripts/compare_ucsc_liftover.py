#!/usr/bin/env python
"""Compare our pipeline's query CNEE annotations against UCSC's liftOver set.

Overlap is purely geometric (share any overlap at all on the same query
chromosome, not necessarily an exact-match interval, and independent of the
two files' different naming schemes -- UCSC's bed names by the *reference*
cnee_id, ours by the graph-assigned query_cnee_id).

Reports:
  - total overlapping annotations: our query CNEEs that overlap >=1 UCSC entry
  - unique to UCSC: UCSC entries with no overlap in our query bed at all
  - unique to our method: our query CNEEs with no overlap in the UCSC bed
  - a breakdown of the overlapping set by orthology category (joined from
    the corresponding overlap_stats.orthology.tsv)

Example:
    PYTHONPATH=src python scripts/compare_ucsc_liftover.py \\
        --ucsc-bed tests/HLcolLivi3.ucscLiftover.bed \\
        --query-bed tests/test_runs/colLivi/query_cnees.bed \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv \\
        --out tests/test_runs/colLivi/ucsc_liftover_comparison.tsv
"""

from __future__ import annotations

import argparse
import csv
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

CATEGORY_ORDER = ["one2one", "one2one(weak)", "one2many", "many2one", "complex/many2many"]


@dataclass(frozen=True)
class BedEntry:
    chrom: str
    start: int
    end: int
    name: str


def _read_bed(path: str | Path) -> list[BedEntry]:
    entries: list[BedEntry] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            name = fields[3] if len(fields) > 3 and fields[3] else ""
            entries.append(BedEntry(fields[0], int(fields[1]), int(fields[2]), name))
    return entries


class IntervalIndex:
    """Overlap-only lookup, assuming intervals within one file don't overlap
    each other (same non-overlapping-annotation assumption used elsewhere in
    this codebase) -- sort by start and bisect on end, same pattern as
    synteny.py's BlockIndex.overlapping().
    """

    def __init__(self, entries: list[BedEntry]) -> None:
        by_chrom: dict[str, list[BedEntry]] = defaultdict(list)
        for entry in entries:
            by_chrom[entry.chrom].append(entry)
        self.by_chrom: dict[str, list[BedEntry]] = {}
        self.ends_by_chrom: dict[str, list[int]] = {}
        for chrom, chrom_entries in by_chrom.items():
            ordered = sorted(chrom_entries, key=lambda e: (e.start, e.end))
            self.by_chrom[chrom] = ordered
            self.ends_by_chrom[chrom] = [e.end for e in ordered]

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        entries = self.by_chrom.get(chrom)
        if not entries:
            return False
        ends = self.ends_by_chrom[chrom]
        index = bisect_right(ends, start)
        while index < len(entries):
            entry = entries[index]
            if entry.start >= end:
                break
            if entry.end > start:
                return True
            index += 1
        return False


def _read_query_categories(orthology_path: str | Path) -> dict[str, str]:
    category_by_query_id: dict[str, str] = {}
    with Path(orthology_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["query_cnee_id"] != "NA":
                category_by_query_id[row["query_cnee_id"]] = row["orthology_category"]
    return category_by_query_id


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucsc-bed", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path, help="Our pipeline's query_cnees.bed")
    parser.add_argument("--orthology", required=True, type=Path, help="overlap_stats.orthology.tsv")
    parser.add_argument("--out", type=Path, default=None, help="Optional TSV output path")
    args = parser.parse_args(argv)

    ucsc_entries = _read_bed(args.ucsc_bed)
    ours_entries = _read_bed(args.query_bed)
    ucsc_index = IntervalIndex(ucsc_entries)
    ours_index = IntervalIndex(ours_entries)

    overlapping_ours = [e for e in ours_entries if ucsc_index.overlaps(e.chrom, e.start, e.end)]
    unique_to_ours = [e for e in ours_entries if not ucsc_index.overlaps(e.chrom, e.start, e.end)]
    unique_to_ucsc = [e for e in ucsc_entries if not ours_index.overlaps(e.chrom, e.start, e.end)]

    category_by_query_id = _read_query_categories(args.orthology)
    category_counts = Counter(category_by_query_id.get(e.name, "unknown") for e in overlapping_ours)

    rows: list[tuple[str, int]] = [
        ("Total overlapping annotations", len(overlapping_ours)),
        ("Unique to UCSC", len(unique_to_ucsc)),
        ("Unique to our method", len(unique_to_ours)),
    ]
    seen_categories = set(category_counts)
    for category in CATEGORY_ORDER:
        if category in seen_categories:
            rows.append((f"  {category}", category_counts[category]))
            seen_categories.discard(category)
    for category in sorted(seen_categories):
        rows.append((f"  {category}", category_counts[category]))

    label_width = max(len(label) for label, _ in rows)
    for label, count in rows:
        print(f"{label.ljust(label_width)}\t{count}")

    if args.out is not None:
        with args.out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for label, count in rows:
                writer.writerow([label.strip(), count])
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
