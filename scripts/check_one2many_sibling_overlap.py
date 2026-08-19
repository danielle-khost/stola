#!/usr/bin/env python
"""For one2many clusters, check whether the "unique to our method" copies are
genuinely novel loci, or just extra paralogous projections sitting alongside
a sibling copy that UCSC's liftOver did call.

UCSC's liftOver only ever reports one locus per reference CNEE, so a ref CNEE
that our pipeline resolved as one2many (multiple query copies) can have at
most one of its query copies overlap the UCSC bed. This groups our one2many
query CNEEs by their shared partner ref CNEE and reports, per cluster,
whether any sibling copy overlaps UCSC at all.

Example:
    PYTHONPATH=src python scripts/check_one2many_sibling_overlap.py \\
        --ucsc-bed tests/HLcolLivi3.ucscLiftover.bed \\
        --query-bed tests/test_runs/colLivi/query_cnees.bed \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv
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


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucsc-bed", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path, help="Our pipeline's query_cnees.bed")
    parser.add_argument("--orthology", required=True, type=Path, help="overlap_stats.orthology.tsv")
    parser.add_argument("--category", default="one2many", help="Orthology category to cluster by ref (default: one2many)")
    args = parser.parse_args(argv)

    ucsc_entries = _read_bed(args.ucsc_bed)
    ours_entries = _read_bed(args.query_bed)
    ucsc_index = IntervalIndex(ucsc_entries)
    ours_by_name = {e.name: e for e in ours_entries}

    query_ids_by_ref: dict[str, list[str]] = defaultdict(list)
    with Path(args.orthology).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["orthology_category"] != args.category:
                continue
            ref_ids = [r for r in row["partner_ref_cnees"].split(",") if r]
            for ref_id in ref_ids:
                query_ids_by_ref[ref_id].append(row["query_cnee_id"])

    cluster_status = Counter()
    query_count_by_status = Counter()
    missing_bed_entries = 0

    for ref_id, query_ids in query_ids_by_ref.items():
        overlap_flags = []
        for qid in query_ids:
            entry = ours_by_name.get(qid)
            if entry is None:
                missing_bed_entries += 1
                continue
            overlap_flags.append(ucsc_index.overlaps(entry.chrom, entry.start, entry.end))
        any_overlap = any(overlap_flags)
        status = "has_ucsc_sibling" if any_overlap else "fully_novel_cluster"
        cluster_status[status] += 1
        query_count_by_status[status] += len(query_ids)

    total_clusters = sum(cluster_status.values())
    total_queries = sum(query_count_by_status.values())

    print(f"'{args.category}' ref clusters: {total_clusters} (covering {total_queries} query CNEEs)\n")
    for status, count in cluster_status.most_common():
        qcount = query_count_by_status[status]
        print(
            f"  {status.ljust(20)}\t{count} clusters ({100*count/total_clusters:.1f}%)"
            f"\t-> {qcount} query CNEEs ({100*qcount/total_queries:.1f}%)"
        )

    has_sibling_queries = query_count_by_status["has_ucsc_sibling"]
    novel_cluster_queries = query_count_by_status["fully_novel_cluster"]
    print(
        f"\nOf {total_queries} query CNEEs in '{args.category}' clusters: {has_sibling_queries} sit in a "
        f"cluster where a sibling copy overlaps UCSC (so most of those {has_sibling_queries} are the 'extra' "
        f"paralog copies beyond the one UCSC could call), while {novel_cluster_queries} are in clusters where "
        "NO copy overlaps UCSC at all (genuinely novel loci, not just extra paralogs)."
    )
    if missing_bed_entries:
        print(f"\nwarning: {missing_bed_entries} query_cnee_id(s) from --orthology had no entry in --query-bed", flush=True)


if __name__ == "__main__":
    main()
