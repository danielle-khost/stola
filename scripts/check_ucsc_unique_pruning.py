#!/usr/bin/env python
"""For UCSC liftOver entries with no overlap in our final query bed, check
whether our pipeline actually produced projection(s) for that reference CNEE
and, if so, what final orthology category it landed in -- in particular,
whether it was pruned to "uncertain_loss" during synteny verification.

UCSC's liftOver bed names each entry by the *reference* cnee_id, so the
"unique to UCSC" set from compare_ucsc_liftover.py can be joined directly
against overlap_stats.orthology.tsv's partner_ref_cnees column.

Example:
    PYTHONPATH=src python scripts/check_ucsc_unique_pruning.py \\
        --ucsc-bed tests/HLcolLivi3.ucscLiftover.bed \\
        --query-bed tests/test_runs/colLivi/query_cnees.bed \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv \\
        --out tests/test_runs/colLivi/ucsc_unique_pruning.tsv
"""

from __future__ import annotations

import argparse
import csv
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


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


def _final_category_by_ref(orthology_path: str | Path) -> dict[str, str]:
    category_by_ref: dict[str, str] = {}
    with Path(orthology_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"partner_ref_cnees", "orthology_category"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Orthology file is missing columns: {', '.join(missing)}")
        for row in reader:
            category = row["orthology_category"]
            for ref_id in row["partner_ref_cnees"].split(","):
                if ref_id:
                    category_by_ref[ref_id] = category
    return category_by_ref


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucsc-bed", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path, help="Our pipeline's query_cnees.bed")
    parser.add_argument("--orthology", required=True, type=Path, help="overlap_stats.orthology.tsv")
    parser.add_argument("--out", type=Path, default=None, help="Optional per-CNEE TSV output path")
    args = parser.parse_args(argv)

    ucsc_entries = _read_bed(args.ucsc_bed)
    ours_entries = _read_bed(args.query_bed)
    ours_index = IntervalIndex(ours_entries)

    unique_to_ucsc = [e for e in ucsc_entries if not ours_index.overlaps(e.chrom, e.start, e.end)]

    category_by_ref = _final_category_by_ref(args.orthology)

    rows: list[tuple[str, str]] = []
    for entry in unique_to_ucsc:
        category = category_by_ref.get(entry.name, "not_in_orthology_output")
        rows.append((entry.name, category))

    counts = Counter(category for _, category in rows)

    print(f"UCSC entries with no overlap in our final query bed: {len(rows)}\n")
    label_width = max(len(cat) for cat in counts) if counts else 0
    for category, count in counts.most_common():
        pct = 100 * count / len(rows) if rows else 0.0
        print(f"  {category.ljust(label_width)}\t{count}\t({pct:.1f}%)")

    uncertain = counts.get("uncertain_loss", 0)
    if rows:
        print(
            f"\n{uncertain}/{len(rows)} ({100 * uncertain / len(rows):.1f}%) of UCSC-unique CNEEs "
            "were pruned to uncertain_loss during synteny verification -- i.e. our pipeline "
            "did produce a projection for these, but it failed the synteny checks and got dropped."
        )

    if args.out is not None:
        with args.out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["ref_cnee_id", "final_orthology_category"])
            writer.writerows(rows)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
