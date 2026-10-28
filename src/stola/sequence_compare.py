"""Sequence comparison for chain-projected CNEEs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .genome_io import open_genome_accessor, reverse_complement
from .models import ChainBlock

SEQUENCE_COLUMNS = [
    "seq_ref_len",
    "seq_query_len",
    "seq_compared_bases",
    "seq_matches",
    "seq_mismatches",
    "seq_identity",
]

MISSING_SEQUENCE_METRICS: dict[str, str | int] = {
    "seq_ref_len": ".",
    "seq_query_len": ".",
    "seq_compared_bases": 0,
    "seq_matches": 0,
    "seq_mismatches": 0,
    "seq_identity": ".",
}


@dataclass(frozen=True)
class AlignedSegment:
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int


class SequenceComparator:
    """Compare reference and query sequence implied by chain aligned blocks."""

    def __init__(
        self,
        ref_fasta: str | Path,
        query_fasta: str | Path,
        work_dir: str | Path,
    ) -> None:
        # Prepare random-access readers for the reference and query genome
        # files -- each independently either FASTA or .2bit (see genome_io.py).
        work_path = Path(work_dir)
        work_path.mkdir(parents=True, exist_ok=True)
        self.ref = open_genome_accessor(ref_fasta, work_path / "ref_fasta.idx")
        self.query = open_genome_accessor(query_fasta, work_path / "query_fasta.idx")

    def close(self) -> None:
        # Close both genome readers held by the comparator.
        self.ref.close()
        self.query.close()

    def compare(
        self,
        ref_chrom: str,
        ref_start: int,
        ref_end: int,
        query_chrom: str,
        query_strand: str,
        aligned_segments: list[AlignedSegment],
    ) -> dict[str, str | int | float]:
        # Compare chain-aligned reference/query sequence and summarize identity.
        ref_len = ref_end - ref_start
        if not aligned_segments:
            return {
                "seq_ref_len": ref_len,
                "seq_query_len": 0,
                "seq_compared_bases": 0,
                "seq_matches": 0,
                "seq_mismatches": 0,
                "seq_identity": ".",
            }

        ref_parts: list[str] = []
        query_parts: list[str] = []
        for segment in aligned_segments:
            ref_part = self.ref.fetch(ref_chrom, segment.ref_start, segment.ref_end)
            query_part = self.query.fetch(
                query_chrom, segment.query_start, segment.query_end
            )
            if query_strand == "-":
                query_part = reverse_complement(query_part)
            ref_parts.append(ref_part)
            query_parts.append(query_part)

        ref_seq = "".join(ref_parts).upper()
        query_seq = "".join(query_parts).upper()
        compared = min(len(ref_seq), len(query_seq))
        matches = sum(1 for a, b in zip(ref_seq[:compared], query_seq[:compared]) if a == b)
        mismatches = compared - matches
        return {
            "seq_ref_len": ref_len,
            "seq_query_len": len(query_seq),
            "seq_compared_bases": compared,
            "seq_matches": matches,
            "seq_mismatches": mismatches,
            "seq_identity": (matches / compared) if compared else ".",
        }


def aligned_segments_for_interval(
    blocks: Iterable[ChainBlock],
    ref_start: int,
    ref_end: int,
    query_strand: str,
) -> list[AlignedSegment]:
    # Project each aligned block overlap into paired reference/query segments.
    segments: list[AlignedSegment] = []
    for block in blocks:
        overlap_start = max(ref_start, block.ref_start)
        overlap_end = min(ref_end, block.ref_end)
        if overlap_end <= overlap_start:
            continue
        offset_start = overlap_start - block.ref_start
        offset_end = overlap_end - block.ref_start
        if query_strand == "+":
            query_start = block.query_start + offset_start
            query_end = block.query_start + offset_end
        else:
            query_start = block.query_end - offset_end
            query_end = block.query_end - offset_start
        segments.append(
            AlignedSegment(
                ref_start=overlap_start,
                ref_end=overlap_end,
                query_start=min(query_start, query_end),
                query_end=max(query_start, query_end),
            )
        )
    return segments
