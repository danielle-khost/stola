"""Block-level overlap statistics for candidate chain/CNEE pairs."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from typing import Iterator
from pathlib import Path

from .chain import RefBlockIndex, parse_chain_blocks, read_chain_text
from .candidates import CANDIDATE_COLUMNS
from .index import ChainIndex, as_chain_index
from .intervals import fraction
from .models import ChainRecord
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
    "num_candidate_chains_for_cnee",
    "rank_by_chain_score",
] + SEQUENCE_COLUMNS


def read_candidates(path: str | Path) -> list[dict[str, str]]:
    # Load the candidate TSV and verify it has the expected columns.
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in CANDIDATE_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"Candidate file is missing columns: {', '.join(missing)}")
        return list(reader)


def _merge_segments(segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    # Collapse overlapping or touching query segments into nonredundant spans.
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
    # Format merged query segments as BED-style block counts, sizes, and starts.
    if query_start is None or not query_segments:
        return 0, ".", "."
    sizes = ",".join(str(end - start) for start, end in query_segments) + ","
    starts = ",".join(str(start - query_start) for start, _ in query_segments) + ","
    return len(query_segments), sizes, starts


def _rank_candidates(rows: list[dict[str, str]]) -> dict[tuple[str, str], int]:
    # Rank each CNEE's candidate chains by score, breaking ties by chain ID.
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
    chain_index: ChainIndex | str | Path,
    candidate_path: str | Path,
    sequence_comparator: SequenceComparator | None = None,
) -> Iterator[dict[str, str | int | float]]:
    # Stream one overlap-stat row per candidate chain/CNEE pair.
    chain_index = as_chain_index(chain_index)
    record_by_id: dict[str, ChainRecord] = {
        record.chain_id: record for record in chain_index.records
    }
    candidates = read_candidates(candidate_path)
    candidate_counts = Counter(row["cnee_uid"] for row in candidates)
    ranks = _rank_candidates(candidates)
    block_cache: dict[str, RefBlockIndex] = {}
    cached_chrom: str | None = None

    with chain_index.indexed_chain_path.open("rb") as chain_handle:
        for candidate in candidates:
            if candidate["ref_chrom"] != cached_chrom:
                # Candidates are sorted by reference chromosome, and a chain
                # aligns to a single reference chromosome, so the previous
                # chromosome's chains are never needed again.
                block_cache.clear()
                cached_chrom = candidate["ref_chrom"]
            chain_id = candidate["chain_id"]
            record = record_by_id[chain_id]
            if chain_id not in block_cache:
                block_cache[chain_id] = RefBlockIndex.from_query_blocks(
                    parse_chain_blocks(read_chain_text(chain_handle, record))
                )
            blocks = block_cache[chain_id]

            ref_start = int(candidate["ref_start"])
            ref_end = int(candidate["ref_end"])
            cnee_len = ref_end - ref_start
            span_overlap = int(candidate["span_overlap_bases"])
            # Each aligned block's overlap with the CNEE, paired with its
            # query-side projection.
            aligned_segments = aligned_segments_for_interval(
                blocks.overlapping(ref_start, ref_end), ref_start, ref_end, record.query_strand
            )
            aligned = sum(segment.ref_end - segment.ref_start for segment in aligned_segments)
            query_segments = _merge_segments(
                [(segment.query_start, segment.query_end) for segment in aligned_segments]
            )
            query_start = query_segments[0][0] if query_segments else None
            query_end = query_segments[-1][1] if query_segments else None
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
                    aligned_segments,
                )

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
                    "num_candidate_chains_for_cnee": candidate_counts[
                        candidate["cnee_uid"]
                    ],
                    "rank_by_chain_score": ranks[(candidate["cnee_uid"], chain_id)],
                    **sequence_metrics,
                }
            )
            yield output_row


def calculate_overlap_rows(
    chain_index: ChainIndex | str | Path,
    candidate_path: str | Path,
    sequence_comparator: SequenceComparator | None = None,
) -> list[dict[str, str | int | float]]:
    # Materialize overlap rows for tests and small programmatic callers.
    return list(iter_overlap_rows(chain_index, candidate_path, sequence_comparator))


def write_overlap_stats(
    chain_index: ChainIndex | str | Path,
    candidate_path: str | Path,
    out_path: str | Path,
    sequence_comparator: SequenceComparator | None = None,
) -> int:
    # Stream overlap-stat rows straight to disk from an index and candidate table.
    rows = iter_overlap_rows(chain_index, candidate_path, sequence_comparator)
    count = 0
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OVERLAP_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count
