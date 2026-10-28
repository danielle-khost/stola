"""UCSC chain parsing and byte-offset extraction."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import BinaryIO, Iterator

import numpy as np

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


def read_chain_text(handle: BinaryIO, record: ChainRecord) -> str:
    # Read one complete chain by byte offset from an open (binary) handle on
    # the indexed plain chain file.
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

    @classmethod
    def from_blocks(cls, blocks: list[ChainBlock]) -> "BlockIndex":
        # Sort blocks and cache reference ends for fast overlap lookups.
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


def _flanking_blocks(
    blocks: BlockIndex, start: int, end: int
) -> tuple[ChainBlock | None, ChainBlock | None]:
    # Find the aligned blocks immediately preceding/following a reference
    # interval that itself overlaps no aligned block -- i.e. a "lost"
    # candidate sitting fully inside a chain gap.
    preceding_index = bisect_right(blocks.ref_ends, start) - 1
    following_index = bisect_left(blocks.blocks, end, key=lambda block: block.ref_start)
    preceding = blocks.blocks[preceding_index] if preceding_index >= 0 else None
    following = (
        blocks.blocks[following_index] if following_index < len(blocks.blocks) else None
    )
    return preceding, following


@dataclass(frozen=True)
class QueryBlocks:
    """A chain's aligned blocks, sorted by query position."""

    query_starts: np.ndarray
    query_ends: np.ndarray
    ref_starts: np.ndarray
    ref_ends: np.ndarray

    def _overlap_range(self, start: int, end: int) -> tuple[int, int]:
        if end <= start or self.query_starts.size == 0:
            return 0, 0
        lo = int(np.searchsorted(self.query_ends, start, side="right"))
        hi = int(np.searchsorted(self.query_starts, end, side="left"))
        return lo, max(lo, hi)

    def overlapping(self, start: int, end: int) -> list[ChainBlock]:
        lo, hi = self._overlap_range(start, end)
        return [
            ChainBlock(
                ref_start=int(self.ref_starts[i]), ref_end=int(self.ref_ends[i]),
                query_start=int(self.query_starts[i]), query_end=int(self.query_ends[i]),
            )
            for i in range(lo, hi)
        ]

    def aligned_bases(self, start: int, end: int) -> int:
        """Query bases in [start, end) covered by this chain's aligned blocks."""
        lo, hi = self._overlap_range(start, end)
        if hi <= lo:
            return 0
        clipped_starts = np.maximum(self.query_starts[lo:hi], start)
        clipped_ends = np.minimum(self.query_ends[lo:hi], end)
        return int(np.clip(clipped_ends - clipped_starts, 0, None).sum())


def parse_query_blocks_vectorized(header_fields: list[str], body_lines: list[str]) -> QueryBlocks:
    """Vectorized equivalent of iter_chain_blocks.

    iter_chain_blocks walks the chain body one line at a time in pure
    Python, building a ChainBlock dataclass instance per line. When the
    multi-projection filter used it, profiling a real genome-wide run showed
    that loop consuming ~80% of total wall time -- almost entirely
    interpreter overhead -- and its blocks take ~8x the memory of these
    arrays. (Lost-ref classification still uses iter_chain_blocks.)

    This computes the exact same (ref_start, ref_end, query_start, query_end)
    per block, but as whole-chain numpy arrays: every body line's (size, dt,
    dq) is tokenized in ONE big str.split() (the ragged last line, which
    omits dt/dq, is parsed separately) and cast to an int64 array in one
    shot, then both coordinate axes are obtained via np.cumsum instead of a
    running Python counter. Blocks are monotonic along each axis (a single
    chain's blocks never reorder or overlap), so getting them into ascending
    query order needs at most one np.flip -- never an O(n log n) sort.
    """
    ref_strand = header_fields[4]
    ref_size = int(header_fields[3])
    ref_raw_start = int(header_fields[5])
    query_strand = header_fields[9]
    query_size = int(header_fields[8])
    query_raw_start = int(header_fields[10])

    n = len(body_lines)
    if n == 0:
        raise ValueError("Cannot parse a chain with no body lines")

    if n > 1:
        tokens = " ".join(body_lines[:-1]).split()
        regular = np.array(tokens, dtype=np.int64).reshape(-1, 3)
        sizes, dts, dqs = regular[:, 0], regular[:, 1], regular[:, 2]
    else:
        sizes = np.empty(0, dtype=np.int64)
        dts = np.empty(0, dtype=np.int64)
        dqs = np.empty(0, dtype=np.int64)

    last_size = int(body_lines[-1].split()[0])
    sizes = np.concatenate([sizes, [last_size]])
    dts = np.concatenate([dts, [0]])
    dqs = np.concatenate([dqs, [0]])

    prefix_ref = np.concatenate(([0], np.cumsum(sizes + dts)[:-1]))
    prefix_query = np.concatenate(([0], np.cumsum(sizes + dqs)[:-1]))

    if ref_strand == "+":
        ref_start = ref_raw_start + prefix_ref
        ref_end = ref_start + sizes
    elif ref_strand == "-":
        ref_end = (ref_size - ref_raw_start) - prefix_ref
        ref_start = ref_end - sizes
    else:
        raise ValueError(f"Invalid chain target strand: {ref_strand}")

    if query_strand == "+":
        query_start = query_raw_start + prefix_query
        query_end = query_start + sizes
    elif query_strand == "-":
        query_end = (query_size - query_raw_start) - prefix_query
        query_start = query_end - sizes
    else:
        raise ValueError(f"Invalid chain query strand: {query_strand}")

    if query_strand == "-":
        ref_start, ref_end, query_start, query_end = (
            arr[::-1] for arr in (ref_start, ref_end, query_start, query_end)
        )
    return QueryBlocks(query_start, query_end, ref_start, ref_end)


def parse_chain_blocks(chain_text: str) -> QueryBlocks:
    """Parse one chain's text (header plus body) with the vectorized parser."""
    lines = [line for line in chain_text.splitlines() if line]
    if not lines:
        raise ValueError("Cannot parse an empty chain")
    return parse_query_blocks_vectorized(lines[0].split(), lines[1:])


@dataclass(frozen=True)
class RefBlockIndex:
    """A chain's aligned blocks as numpy arrays in reference order, for
    reference-interval lookups (QueryBlocks is the query-ordered version)."""

    ref_starts: np.ndarray
    ref_ends: np.ndarray
    query_starts: np.ndarray
    query_ends: np.ndarray

    @classmethod
    def from_query_blocks(cls, blocks: QueryBlocks) -> "RefBlockIndex":
        # A chain's blocks are monotonic along both axes, so reference order
        # is query order, reversed when the two axes run in opposite directions.
        arrays = (blocks.ref_starts, blocks.ref_ends, blocks.query_starts, blocks.query_ends)
        if blocks.ref_starts.size > 1 and blocks.ref_starts[0] > blocks.ref_starts[-1]:
            arrays = tuple(array[::-1] for array in arrays)
        return cls(*arrays)

    def overlapping(self, start: int, end: int) -> list[ChainBlock]:
        # Blocks whose reference intervals overlap [start, end), in reference order.
        lo = int(np.searchsorted(self.ref_ends, start, side="right"))
        hi = int(np.searchsorted(self.ref_starts, end, side="left"))
        return [
            ChainBlock(
                ref_start=int(self.ref_starts[i]), ref_end=int(self.ref_ends[i]),
                query_start=int(self.query_starts[i]), query_end=int(self.query_ends[i]),
            )
            for i in range(lo, hi)
        ]
