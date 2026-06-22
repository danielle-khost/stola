"""Block-level overlap statistics for candidate chain/CNEE pairs."""

from __future__ import annotations

import csv
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterator
from pathlib import Path

from .chain import iter_chain_blocks
from .candidates import CANDIDATE_COLUMNS
from .index import load_index
from .intervals import fraction, overlap_length
from .models import ChainBlock, ChainRecord
from .sequence_compare import (
    MISSING_SEQUENCE_METRICS,
    SEQUENCE_COLUMNS,
    SequenceComparator,
    aligned_segments_for_interval,
)

OVERLAP_COLUMNS = CANDIDATE_COLUMNS + [
    "query_projection_name",
    "cnee_len",
    "aligned_bases",
    "covered_fraction",
    "gap_bases_within_chain_span",
    "gap_fraction_within_chain_span",
    "query_aligned_start",
    "query_aligned_end",
    "query_block_count",
    "query_block_sizes",
    "query_block_starts",
    "left_flank_len",
    "left_flank_aligned_bases",
    "left_flank_covered_fraction",
    "right_flank_len",
    "right_flank_aligned_bases",
    "right_flank_covered_fraction",
    "num_candidate_chains_for_cnee",
    "rank_by_chain_score",
] + SEQUENCE_COLUMNS


@dataclass(frozen=True)
class BlockIndex:
    blocks: list[ChainBlock]
    ref_ends: list[int]

    @classmethod
    def from_blocks(cls, blocks: list[ChainBlock]) -> "BlockIndex":
        ordered = sorted(blocks, key=lambda block: (block.ref_start, block.ref_end))
        return cls(ordered, [block.ref_end for block in ordered])

    def overlapping(self, start: int, end: int) -> Iterator[ChainBlock]:
        index = bisect_right(self.ref_ends, start)
        while index < len(self.blocks):
            block = self.blocks[index]
            if block.ref_start >= end:
                break
            yield block
            index += 1


def read_candidates(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in CANDIDATE_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"Candidate file is missing columns: {', '.join(missing)}")
        return list(reader)


def _project_overlap_to_query(
    block: ChainBlock,
    query_strand: str,
    ref_overlap_start: int,
    ref_overlap_end: int,
) -> tuple[int, int]:
    offset_start = ref_overlap_start - block.ref_start
    offset_end = ref_overlap_end - block.ref_start
    if query_strand == "+":
        query_start = block.query_start + offset_start
        query_end = block.query_start + offset_end
    else:
        query_start = block.query_end - offset_end
        query_end = block.query_end - offset_start
    return min(query_start, query_end), max(query_start, query_end)


def _aligned_bases(
    blocks: BlockIndex,
    start: int,
    end: int,
    query_strand: str | None = None,
) -> tuple[int, int | None, int | None, list[tuple[int, int]]]:
    aligned = 0
    query_segments: list[tuple[int, int]] = []
    for block in blocks.overlapping(start, end):
        overlap = overlap_length(start, end, block.ref_start, block.ref_end)
        if overlap == 0:
            continue
        ref_overlap_start = max(start, block.ref_start)
        ref_overlap_end = min(end, block.ref_end)
        aligned += overlap
        if query_strand is not None:
            query_segments.append(
                _project_overlap_to_query(
                    block, query_strand, ref_overlap_start, ref_overlap_end
                )
            )
    if not query_segments:
        return aligned, None, None, []
    merged_segments = _merge_segments(query_segments)
    return (
        aligned,
        min(segment[0] for segment in merged_segments),
        max(segment[1] for segment in merged_segments),
        merged_segments,
    )


def _merge_segments(segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not segments:
        return []
    ordered = sorted(segments)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _block_fields(
    query_start: int | None, query_segments: list[tuple[int, int]]
) -> tuple[int, str, str]:
    if query_start is None or not query_segments:
        return 0, ".", "."
    sizes = ",".join(str(end - start) for start, end in query_segments) + ","
    starts = ",".join(str(start - query_start) for start, _ in query_segments) + ","
    return len(query_segments), sizes, starts


def _rank_candidates(rows: list[dict[str, str]]) -> dict[tuple[str, str], int]:
    grouped: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["cnee_uid"]].append(row)

    ranks: dict[tuple[str, str], int] = {}
    for cnee_uid, group in grouped.items():
        ordered = sorted(
            group,
            key=lambda row: (-int(row["chain_score"]), str(row["chain_id"])),
        )
        for rank, row in enumerate(ordered, start=1):
            ranks[(cnee_uid, row["chain_id"])] = rank
    return ranks


def iter_overlap_rows(
    index_path: str | Path,
    candidate_path: str | Path,
    flank_size: int = 1000,
    sequence_comparator: SequenceComparator | None = None,
) -> Iterator[dict[str, str | int | float]]:
    if flank_size < 0:
        raise ValueError("flank_size must be non-negative")
    payload, records = load_index(index_path)
    chain_path = payload["indexed_chain_path"]
    record_by_id: dict[str, ChainRecord] = {record.chain_id: record for record in records}
    candidates = read_candidates(candidate_path)
    candidate_counts = Counter(row["cnee_uid"] for row in candidates)
    ranks = _rank_candidates(candidates)
    block_cache: dict[str, BlockIndex] = {}

    with Path(chain_path).open("rb") as chain_handle:
        for candidate in candidates:
            chain_id = candidate["chain_id"]
            record = record_by_id[chain_id]
            if chain_id not in block_cache:
                chain_handle.seek(record.start_byte)
                chain_text = chain_handle.read(record.byte_length).decode("utf-8")
                block_cache[chain_id] = BlockIndex.from_blocks(
                    iter_chain_blocks(chain_text)
                )
            blocks = block_cache[chain_id]

            ref_start = int(candidate["ref_start"])
            ref_end = int(candidate["ref_end"])
            cnee_len = ref_end - ref_start
            span_overlap = int(candidate["span_overlap_bases"])
            overlapping_blocks = list(blocks.overlapping(ref_start, ref_end))
            cnee_block_index = BlockIndex.from_blocks(overlapping_blocks)
            aligned, query_start, query_end, query_segments = _aligned_bases(
                cnee_block_index, ref_start, ref_end, record.query_strand
            )
            gap_bases = max(0, span_overlap - aligned)
            query_block_count, query_block_sizes, query_block_starts = _block_fields(
                query_start, query_segments
            )
            sequence_metrics: dict[str, str | int | float]
            if sequence_comparator is None:
                sequence_metrics = dict(MISSING_SEQUENCE_METRICS)
            else:
                sequence_metrics = sequence_comparator.compare(
                    str(candidate["ref_chrom"]),
                    ref_start,
                    ref_end,
                    record.query_chrom,
                    record.query_strand,
                    aligned_segments_for_interval(
                        overlapping_blocks, ref_start, ref_end, record.query_strand
                    ),
                )

            left_start = max(0, ref_start - flank_size)
            left_end = ref_start
            right_start = ref_end
            right_end = ref_end + flank_size
            left_len = left_end - left_start
            right_len = right_end - right_start
            left_aligned, _, _, _ = _aligned_bases(blocks, left_start, left_end)
            right_aligned, _, _, _ = _aligned_bases(blocks, right_start, right_end)

            output_row: dict[str, str | int | float] = dict(candidate)
            projection_name = f"{candidate['cnee_uid']}#{chain_id}"
            output_row.update(
                {
                    "query_projection_name": projection_name,
                    "cnee_len": cnee_len,
                    "aligned_bases": aligned,
                    "covered_fraction": fraction(aligned, cnee_len),
                    "gap_bases_within_chain_span": gap_bases,
                    "gap_fraction_within_chain_span": fraction(gap_bases, span_overlap),
                    "query_aligned_start": "." if query_start is None else query_start,
                    "query_aligned_end": "." if query_end is None else query_end,
                    "query_block_count": query_block_count,
                    "query_block_sizes": query_block_sizes,
                    "query_block_starts": query_block_starts,
                    "left_flank_len": left_len,
                    "left_flank_aligned_bases": left_aligned,
                    "left_flank_covered_fraction": fraction(left_aligned, left_len),
                    "right_flank_len": right_len,
                    "right_flank_aligned_bases": right_aligned,
                    "right_flank_covered_fraction": fraction(right_aligned, right_len),
                    "num_candidate_chains_for_cnee": candidate_counts[
                        candidate["cnee_uid"]
                    ],
                    "rank_by_chain_score": ranks[(candidate["cnee_uid"], chain_id)],
                    **sequence_metrics,
                }
            )
            yield output_row


def calculate_overlap_rows(
    index_path: str | Path,
    candidate_path: str | Path,
    flank_size: int = 1000,
) -> list[dict[str, str | int | float]]:
    return list(iter_overlap_rows(index_path, candidate_path, flank_size))


def write_overlap_stats(
    index_path: str | Path,
    bed_path: str | Path,
    candidate_path: str | Path,
    out_path: str | Path,
    flank_size: int = 1000,
) -> int:
    # BED is intentionally accepted here to keep the CLI contract explicit. The
    # candidate table carries the CNEE coordinates used for v1 stats.
    _ = bed_path
    rows = iter_overlap_rows(index_path, candidate_path, flank_size)
    count = 0
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OVERLAP_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count
