"""UCSC chain parsing and byte-offset extraction."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .models import ChainBlock, ChainRecord


def _normalize_start_end(size: int, strand: str, start: int, end: int) -> tuple[int, int]:
    # Convert strand-relative chain coordinates into plus-strand coordinates.
    if strand == "+":
        return start, end
    if strand == "-":
        return size - end, size - start
    raise ValueError(f"Invalid chain strand: {strand}")


def parse_chain_header(header_line: str, start_byte: int = 0, byte_length: int = 0) -> ChainRecord:
    # Parse one UCSC chain header and attach byte-offset extraction metadata.
    fields = header_line.rstrip("\n").split()
    if len(fields) != 13 or fields[0] != "chain":
        raise ValueError(f"Invalid chain header: {header_line.rstrip()}")

    score = int(fields[1])
    ref_chrom = fields[2]
    ref_size = int(fields[3])
    ref_strand = fields[4]
    ref_raw_start = int(fields[5])
    ref_raw_end = int(fields[6])
    query_chrom = fields[7]
    query_size = int(fields[8])
    query_strand = fields[9]
    query_raw_start = int(fields[10])
    query_raw_end = int(fields[11])
    chain_id = fields[12]

    ref_start, ref_end = _normalize_start_end(
        ref_size, ref_strand, ref_raw_start, ref_raw_end
    )
    query_start, query_end = _normalize_start_end(
        query_size, query_strand, query_raw_start, query_raw_end
    )
    return ChainRecord(
        chain_id=chain_id,
        score=score,
        ref_chrom=ref_chrom,
        ref_size=ref_size,
        ref_strand=ref_strand,
        ref_start=ref_start,
        ref_end=ref_end,
        query_chrom=query_chrom,
        query_size=query_size,
        query_strand=query_strand,
        query_start=query_start,
        query_end=query_end,
        start_byte=start_byte,
        byte_length=byte_length,
    )


def extract_chain_text(chain_path: str | Path, record: ChainRecord) -> str:
    # Read one complete chain from the indexed plain chain file by byte offset.
    with Path(chain_path).open("rb") as handle:
        handle.seek(record.start_byte)
        return handle.read(record.byte_length).decode("utf-8")


def iter_chain_blocks(chain_text: str) -> list[ChainBlock]:
    """Return aligned blocks from a single chain as plus-coordinate intervals."""
    # Walk the chain body rows and emit aligned reference/query block intervals.
    lines = [line for line in chain_text.splitlines() if line]
    if not lines:
        raise ValueError("Cannot parse an empty chain")
    header = lines[0].split()
    if len(header) != 13 or header[0] != "chain":
        raise ValueError("Chain text does not begin with a valid chain header")

    ref_size = int(header[3])
    ref_strand = header[4]
    ref_raw_start = int(header[5])
    query_size = int(header[8])
    query_strand = header[9]
    query_raw_start = int(header[10])

    if ref_strand == "+":
        ref_progress = ref_raw_start
    elif ref_strand == "-":
        ref_progress = ref_size - ref_raw_start
    else:
        raise ValueError(f"Invalid chain target strand: {ref_strand}")

    if query_strand == "+":
        query_progress = query_raw_start
    elif query_strand == "-":
        query_progress = query_size - query_raw_start
    else:
        raise ValueError(f"Invalid chain query strand: {query_strand}")

    blocks: list[ChainBlock] = []
    for line_number, line in enumerate(lines[1:], start=2):
        parts = line.split()
        if len(parts) not in {1, 3}:
            raise ValueError(f"Invalid chain block at line {line_number}: {line}")
        size = int(parts[0])
        dt = int(parts[1]) if len(parts) == 3 else 0
        dq = int(parts[2]) if len(parts) == 3 else 0

        if ref_strand == "+":
            ref_start = ref_progress
            ref_end = ref_start + size
            ref_progress = ref_end + dt
        else:
            ref_end = ref_progress
            ref_start = ref_end - size
            ref_progress = ref_start - dt

        if query_strand == "+":
            query_start = query_progress
            query_end = query_start + size
            query_progress = query_end + dq
        else:
            query_end = query_progress
            query_start = query_end - size
            query_progress = query_start - dq

        blocks.append(
            ChainBlock(
                ref_start=min(ref_start, ref_end),
                ref_end=max(ref_start, ref_end),
                query_start=min(query_start, query_end),
                query_end=max(query_start, query_end),
            )
        )
    return blocks


@dataclass(frozen=True)
class BlockIndex:
    blocks: list[ChainBlock]
    ref_ends: list[int]
    ref_starts: list[int]

    @classmethod
    def from_blocks(cls, blocks: list[ChainBlock]) -> "BlockIndex":
        # Sort blocks and cache reference starts/ends for overlap and nearest-
        # block queries.
        ordered = sorted(blocks, key=lambda block: (block.ref_start, block.ref_end))
        return cls(
            ordered,
            [block.ref_end for block in ordered],
            [block.ref_start for block in ordered],
        )

    def overlapping(self, start: int, end: int) -> Iterator[ChainBlock]:
        # Yield blocks overlapping a half-open reference interval.
        index = bisect_right(self.ref_ends, start)
        while index < len(self.blocks):
            block = self.blocks[index]
            if block.ref_start >= end:
                break
            yield block
            index += 1


def _flanking_blocks(
    blocks: BlockIndex, start: int, end: int
) -> tuple[ChainBlock | None, ChainBlock | None]:
    # Find the aligned blocks immediately preceding/following a reference
    # interval that itself overlaps no aligned block -- i.e. a "lost"
    # candidate sitting fully inside a chain gap.
    preceding_index = bisect_right(blocks.ref_ends, start) - 1
    following_index = bisect_left(blocks.ref_starts, end)
    preceding = blocks.blocks[preceding_index] if preceding_index >= 0 else None
    following = (
        blocks.blocks[following_index] if following_index < len(blocks.blocks) else None
    )
    return preceding, following
