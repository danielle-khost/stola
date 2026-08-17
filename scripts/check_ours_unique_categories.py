#!/usr/bin/env python
"""For our query CNEEs with no overlap in UCSC's liftOver bed, break down
which final orthology category they fall into.

Complements check_ucsc_unique_pruning.py, which looks at the opposite
direction (UCSC-unique annotations and what happened to them in our
pipeline). Here the join is on query_cnee_id, since our query bed names
each entry by the graph-assigned query_cnee_id -- the same key used in
overlap_stats.orthology.tsv.

Example:
    PYTHONPATH=src python scripts/check_ours_unique_categories.py \\
        --ucsc-bed tests/HLcolLivi3.ucscLiftover.bed \\
        --query-bed tests/test_runs/colLivi/query_cnees.bed \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv \\
        --out tests/test_runs/colLivi/ours_unique_categories.tsv
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
    parser.add_argument("--out", type=Path, default=None, help="Optional per-CNEE TSV output path")
    args = parser.parse_args(argv)

    ucsc_entries = _read_bed(args.ucsc_bed)
    ours_entries = _read_bed(args.query_bed)
    ucsc_index = IntervalIndex(ucsc_entries)

    unique_to_ours = [e for e in ours_entries if not ucsc_index.overlaps(e.chrom, e.start, e.end)]

    category_by_query_id = _read_query_categories(args.orthology)

    rows: list[tuple[str, str]] = []
    for entry in unique_to_ours:
        category = category_by_query_id.get(entry.name, "not_in_orthology_output")
        rows.append((entry.name, category))

    counts = Counter(category for _, category in rows)

    print(f"Our query CNEEs with no overlap in UCSC's liftOver bed: {len(rows)}\n")
    label_width = max(len(cat) for cat in counts) if counts else 0
    seen = set(counts)
    ordered_categories = [c for c in CATEGORY_ORDER if c in seen] + sorted(seen.difference(CATEGORY_ORDER))
    for category in ordered_categories:
        count = counts[category]
        pct = 100 * count / len(rows) if rows else 0.0
        print(f"  {category.ljust(label_width)}\t{count}\t({pct:.1f}%)")

    if args.out is not None:
        with args.out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["query_cnee_id", "orthology_category"])
            writer.writerows(rows)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
