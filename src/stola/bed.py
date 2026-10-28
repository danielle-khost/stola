"""BED parsing for CNEE intervals."""

from __future__ import annotations

import warnings
from collections import defaultdict
from pathlib import Path

from .models import BedInterval


def parse_bed(path: str | Path) -> list[BedInterval]:
    """Parse BED3+ records as single CNEE intervals.

    Column 4 is used as the display CNEE ID when present. Otherwise, a stable
    coordinate-based ID is synthesized. The cnee_uid is unique within the file;
    downstream orthology/summary logic keys on cnee_uid specifically so
    that two intervals sharing the same display cnee_id are never conflated
    into one reference locus.
    """
    # Load BED rows into stable CNEE interval objects used across the workflow.
    intervals: list[BedInterval] = []
    seen_ids: defaultdict[str, int] = defaultdict(int)
    duplicated_ids: set[str] = set()
    zero_length_ids: list[str] = []
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
            if seen_ids[cnee_id] > 1:
                duplicated_ids.add(cnee_id)
            cnee_uid = (
                cnee_id if seen_ids[cnee_id] == 1 else f"{cnee_id}__line{line_number}"
            )
            if end == start:
                zero_length_ids.append(cnee_uid)
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
    if duplicated_ids:
        shown = ", ".join(sorted(duplicated_ids)[:10])
        more = ", ..." if len(duplicated_ids) > 10 else ""
        warnings.warn(
            f"{path}: {len(duplicated_ids)} cnee_id value(s) are not unique "
            f"({shown}{more}). These CNEEs are distinct genomic intervals that "
            "happen to share a display name; downstream output keys on the "
            "unique cnee_uid (e.g. 'name__line123') to keep them from being "
            "treated as the same reference locus.",
            stacklevel=2,
        )
    if zero_length_ids:
        shown = ", ".join(sorted(zero_length_ids)[:10])
        more = ", ..." if len(zero_length_ids) > 10 else ""
        warnings.warn(
            f"{path}: {len(zero_length_ids)} CNEE(s) have zero-length intervals "
            f"({shown}{more}). A zero-length interval can never overlap a chain "
            "and will always be reported as missing.",
            stacklevel=2,
        )
    return intervals


def bed_missing_name_column(path: str | Path) -> bool:
    """Report whether any row of a BED file lacks a usable column 4.

    parse_bed synthesizes a coordinate-based cnee_id for such rows, which
    stays stable across separate pipeline runs on the same reference BED, but
    is otherwise invisible to the user -- callers use this to decide whether
    to write out the synthesized names as a companion "_named.bed" file.
    """
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4 or not fields[3]:
                return True
    return False


def write_named_bed(intervals: list[BedInterval], out_path: str | Path) -> int:
    """Write reference CNEEs as BED4, naming each row with its cnee_uid.

    Used to hand the user a stable, named copy of a reference BED that
    lacked a name column, so that the same synthesized ids can be recognized
    across separate runs (e.g. the same reference CNEEs lifted to different
    query species).
    """
    count = 0
    with Path(out_path).open("w", encoding="utf-8") as handle:
        for interval in sorted(intervals, key=lambda item: item.line_number):
            handle.write(
                "\t".join(
                    [interval.chrom, str(interval.start), str(interval.end), interval.cnee_uid]
                )
                + "\n"
            )
            count += 1
    return count


def group_bed_by_chrom(intervals: list[BedInterval]) -> dict[str, list[BedInterval]]:
    # Group and sort CNEEs by reference chromosome for interval scanning.
    grouped: dict[str, list[BedInterval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.chrom].append(interval)
    for chrom in grouped:
        grouped[chrom].sort(key=lambda item: (item.start, item.end, item.cnee_uid))
    return dict(grouped)
