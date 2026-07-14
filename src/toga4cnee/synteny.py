"""Flanking-synteny verification for weak, lost, and uncertain CNEE calls."""

from __future__ import annotations

import csv
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

from .chain import iter_chain_blocks
from .index import load_index
from .models import ChainBlock, ChainRecord
from .sequence_compare import SequenceComparator, aligned_segments_for_interval
from .summary import stats_row_is_solid

# cnee_uid holds the reference CNEE's cnee_uid for lost/uncertain_loss rows,
# but a query-side ID for weak/one2one rows -- see SyntenyTarget.output_id.
SYNTENY_COLUMNS = [
    "cnee_uid",
    "orthology_status",
    "flank_size",
    "flank_seq_identity",
    "flank_fraction_gap",
]

# Verify synteny at several scales at once so results across scales can be
# compared without rerunning the workflow.
DEFAULT_FLANK_SIZES: tuple[int, ...] = (10000, 1000, 100)


@dataclass(frozen=True)
class QueryBedRecord:
    chrom: str
    start: int
    end: int
    query_id: str


@dataclass(frozen=True)
class SyntenyTarget:
    row: dict[str, str]
    status: str
    output_id: str


@dataclass(frozen=True)
class SyntenyResult:
    rows: int


@dataclass(frozen=True)
class BlockIndex:
    blocks: list[ChainBlock]
    ref_ends: list[int]

    @classmethod
    def from_blocks(cls, blocks: list[ChainBlock]) -> "BlockIndex":
        # Sort blocks and cache reference ends for overlap queries.
        ordered = sorted(blocks, key=lambda block: (block.ref_start, block.ref_end))
        return cls(ordered, [block.ref_end for block in ordered])

    def overlapping(self, start: int, end: int) -> Iterator[ChainBlock]:
        # Yield blocks overlapping a half-open reference interval.
        index = bisect_right(self.ref_ends, start)
        while index < len(self.blocks):
            block = self.blocks[index]
            if block.ref_start >= end:
                break
            yield block
            index += 1


def write_synteny_verification(
    index_path: str | Path,
    stats_path: str | Path,
    query_bed_path: str | Path,
    orthology_path: str | Path,
    out_path: str | Path,
    flank_sizes: Sequence[int] = DEFAULT_FLANK_SIZES,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    max_one2one_controls: int = 50000,
    sequence_comparator: SequenceComparator | None = None,
) -> SyntenyResult:
    """Write flanking synteny metrics for weak/lost/uncertain projections.

    One row is written per candidate per requested flank size (long format),
    so synteny support can be compared across scales in one output file.
    """
    # Select rows to verify, calculate flank metrics at each scale, and write
    # a long-format TSV (one row per candidate x flank size).
    if not flank_sizes:
        raise ValueError("flank_sizes must include at least one value")
    for flank_size in flank_sizes:
        if flank_size < 0:
            raise ValueError("flank_size must be non-negative")
    rows = _read_stats(stats_path)
    targets = _select_synteny_targets(
        rows,
        query_bed_path,
        orthology_path,
        min_sequence_identity=min_sequence_identity,
        max_fraction_gap=max_fraction_gap,
        max_one2one_controls=max_one2one_controls,
    )
    payload, records = load_index(index_path)
    chain_path = payload["indexed_chain_path"]
    record_by_id = {record.chain_id: record for record in records}
    block_cache: dict[str, BlockIndex] = {}
    count = 0

    with Path(chain_path).open("rb") as chain_handle:
        with Path(out_path).open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=SYNTENY_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for target in targets:
                chain_id = target.row["chain_id"]
                record = record_by_id[chain_id]
                if chain_id not in block_cache:
                    chain_handle.seek(record.start_byte)
                    chain_text = chain_handle.read(record.byte_length).decode("utf-8")
                    block_cache[chain_id] = BlockIndex.from_blocks(
                        iter_chain_blocks(chain_text)
                    )
                for flank_size in flank_sizes:
                    flank_metrics = _flank_metrics(
                        target.row,
                        record,
                        block_cache[chain_id],
                        flank_size=flank_size,
                        sequence_comparator=sequence_comparator,
                    )
                    writer.writerow(
                        {
                            "cnee_uid": target.output_id,
                            "orthology_status": target.status,
                            "flank_size": flank_size,
                            "flank_seq_identity": flank_metrics["flank_seq_identity"],
                            "flank_fraction_gap": flank_metrics["flank_fraction_gap"],
                        }
                    )
                    count += 1
    return SyntenyResult(rows=count)


def _read_stats(stats_path: str | Path) -> list[dict[str, str]]:
    # Load overlap stats and validate columns needed for synteny verification.
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_id",
            "cnee_uid",
            "chain_id",
            "ref_chrom",
            "ref_start",
            "ref_end",
            "aligned_bases",
            "query_chrom",
            "query_aligned_start",
            "query_aligned_end",
            "gap_fraction_within_chain_span",
            "seq_identity",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        return list(reader)


def _select_synteny_targets(
    rows: list[dict[str, str]],
    query_bed_path: str | Path,
    orthology_path: str | Path,
    min_sequence_identity: float,
    max_fraction_gap: float,
    max_one2one_controls: int,
) -> list[SyntenyTarget]:
    # Select weak, uncertain_loss, lost, and capped one2one control rows.
    # Reference identity is keyed on cnee_uid (guaranteed unique) rather than
    # cnee_id, so two distinct loci sharing a display name are never conflated.
    refs_with_aligned = {
        row["cnee_uid"] for row in rows if int(row["aligned_bases"]) > 0
    }
    refs_with_solid = {
        row["cnee_uid"]
        for row in rows
        if stats_row_is_solid(
            row,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
    }
    query_lookup = _read_query_bed(query_bed_path)
    one2one_ids = _read_one2one_query_ids(orthology_path)
    targets: list[SyntenyTarget] = []
    one2one_count = 0
    for row in rows:
        cnee_uid = row["cnee_uid"]
        aligned_bases = int(row["aligned_bases"])
        row_is_solid = stats_row_is_solid(
            row,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
        if aligned_bases <= 0:
            if cnee_uid in refs_with_aligned:
                continue
            targets.append(SyntenyTarget(row=row, status="lost", output_id=cnee_uid))
            continue
        if cnee_uid not in refs_with_solid:
            targets.append(
                SyntenyTarget(row=row, status="uncertain_loss", output_id=cnee_uid)
            )
            continue
        if row_is_solid and cnee_uid in one2one_ids and one2one_count < max_one2one_controls:
            targets.append(
                SyntenyTarget(
                    row=row,
                    status="one2one",
                    output_id=one2one_ids[cnee_uid],
                )
            )
            one2one_count += 1
            continue
        if not row_is_solid:
            query_id = _query_id_for_projection(row, query_lookup)
            targets.append(SyntenyTarget(row=row, status="weak", output_id=query_id))
    return targets


def _read_one2one_query_ids(orthology_path: str | Path) -> dict[str, str]:
    # Map one2one reference IDs to their all-mode query CNEE IDs.
    one2one: dict[str, str] = {}
    path = Path(orthology_path)
    if not path.exists():
        return one2one
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_cnee_id", "partner_ref_cnees", "orthology_category"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Orthology file is missing columns: {', '.join(missing)}")
        for row in reader:
            if row["orthology_category"] != "one2one":
                continue
            if row["query_cnee_id"] == "NA":
                continue
            ref_ids = [ref_id for ref_id in row["partner_ref_cnees"].split(",") if ref_id]
            if len(ref_ids) == 1:
                one2one[ref_ids[0]] = row["query_cnee_id"]
    return one2one


def _read_query_bed(path: str | Path) -> dict[str, list[QueryBedRecord]]:
    # Load all-mode query BED records for weak projection ID lookup.
    by_chrom: dict[str, list[QueryBedRecord]] = defaultdict(list)
    bed_path = Path(path)
    if not bed_path.exists():
        return {}
    with bed_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            by_chrom[fields[0]].append(
                QueryBedRecord(
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    query_id=fields[3],
                )
            )
    for records in by_chrom.values():
        records.sort(key=lambda item: (item.start, item.end, item.query_id))
    return dict(by_chrom)


def _query_id_for_projection(
    row: dict[str, str], query_lookup: dict[str, list[QueryBedRecord]]
) -> str:
    # Map a projected query interval back to the all-mode query BED ID.
    query_start = row["query_aligned_start"]
    query_end = row["query_aligned_end"]
    if query_start == "." or query_end == ".":
        return row["cnee_id"]
    start = int(query_start)
    end = int(query_end)
    for record in query_lookup.get(row["query_chrom"], []):
        if start < record.end and end > record.start:
            return record.query_id
    return row["query_projection_name"] if "query_projection_name" in row else row["cnee_id"]


def _flank_metrics(
    row: dict[str, str],
    record: ChainRecord,
    blocks: BlockIndex,
    flank_size: int,
    sequence_comparator: SequenceComparator | None,
) -> dict[str, str | float]:
    # Calculate combined left/right flank identity and gap fraction.
    ref_start = int(row["ref_start"])
    ref_end = int(row["ref_end"])
    intervals = [
        (max(0, ref_start - flank_size), ref_start),
        (ref_end, min(record.ref_size, ref_end + flank_size)),
    ]
    flank_len = sum(end - start for start, end in intervals if end > start)
    if flank_len == 0:
        return {"flank_seq_identity": 0, "flank_fraction_gap": "."}

    aligned = 0
    compared = 0
    matches = 0
    for start, end in intervals:
        if end <= start:
            continue
        overlapping_blocks = list(blocks.overlapping(start, end))
        aligned += _aligned_bases(overlapping_blocks, start, end)
        if sequence_comparator is None or not overlapping_blocks:
            continue
        metrics = sequence_comparator.compare(
            str(row["ref_chrom"]),
            start,
            end,
            record.query_chrom,
            record.query_strand,
            aligned_segments_for_interval(
                overlapping_blocks, start, end, record.query_strand
            ),
        )
        compared += int(metrics["seq_compared_bases"])
        matches += int(metrics["seq_matches"])

    identity = (matches / compared) if compared else 0
    fraction_gap = (flank_len - aligned) / flank_len
    return {
        "flank_seq_identity": identity,
        "flank_fraction_gap": fraction_gap,
    }


def _aligned_bases(blocks: list[ChainBlock], start: int, end: int) -> int:
    # Count reference bases in an interval covered by aligned chain blocks.
    aligned = 0
    for block in blocks:
        aligned += max(0, min(end, block.ref_end) - max(start, block.ref_start))
    return aligned
