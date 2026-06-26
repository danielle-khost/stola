"""BED parsing for CNEE intervals."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .models import BedInterval


def parse_bed(path: str | Path) -> list[BedInterval]:
    """Parse BED3+ records as single CNEE intervals.

    Column 4 is used as the display CNEE ID when present. Otherwise, a stable
    coordinate-based ID is synthesized. The cnee_uid is unique within the file.
    """
    # Load BED rows into stable CNEE interval objects used across the workflow.
    intervals: list[BedInterval] = []
    seen_ids: defaultdict[str, int] = defaultdict(int)
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {line_number} has fewer than 3 columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end < start:
                raise ValueError(f"BED line {line_number} has invalid coordinates")
            cnee_id = (
                fields[3]
                if len(fields) >= 4 and fields[3]
                else f"cnee_{line_number:06d}_{chrom}_{start}_{end}"
            )
            seen_ids[cnee_id] += 1
            cnee_uid = (
                cnee_id if seen_ids[cnee_id] == 1 else f"{cnee_id}__line{line_number}"
            )
            intervals.append(
                BedInterval(
                    cnee_uid=cnee_uid,
                    cnee_id=cnee_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    line_number=line_number,
                )
            )
    return intervals


def group_bed_by_chrom(intervals: list[BedInterval]) -> dict[str, list[BedInterval]]:
    # Group and sort CNEEs by reference chromosome for interval scanning.
    grouped: dict[str, list[BedInterval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.chrom].append(interval)
    for chrom in grouped:
        grouped[chrom].sort(key=lambda item: (item.start, item.end, item.cnee_uid))
    return dict(grouped)
