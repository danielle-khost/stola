"""Sequence comparison for chain-projected CNEEs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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


@dataclass(frozen=True)
class FastaRecordIndex:
    seq_id: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class FastaAccessor:
    """Random-access FASTA reader with conservative accession-version aliases."""

    def __init__(self, fasta_path: str | Path, index_path: str | Path) -> None:
        # Open the FASTA and make sure a reusable byte-offset index is available.
        self.fasta_path = Path(fasta_path)
        self.index_path = Path(index_path)
        if not self.index_path.exists():
            self._build_index()
        self.records = self._read_index()
        self.aliases = self._build_aliases(self.records.keys())
        self.handle = self.fasta_path.open("rb")

    @staticmethod
    def _strip_version(seq_id: str) -> str:
        # Remove a terminal numeric accession version from a sequence ID.
        prefix, dot, suffix = seq_id.rpartition(".")
        if dot and suffix.isdigit():
            return prefix
        return seq_id

    @classmethod
    def _build_aliases(cls, keys: Iterable[str]) -> dict[str, str]:
        # Map exact IDs and unambiguous version-stripped IDs to FASTA records.
        exact = {key: key for key in keys}
        stripped: dict[str, list[str]] = {}
        for key in exact:
            stripped.setdefault(cls._strip_version(key), []).append(key)
        aliases = dict(exact)
        for key, values in stripped.items():
            if key in aliases or len(values) != 1:
                continue
            aliases[key] = values[0]
        return aliases

    def _build_index(self) -> None:
        # Scan the FASTA once and save per-record byte layout metadata.
        self.index_path.parent.mkdir(parents=True, exist_ok=True)
        rows: list[FastaRecordIndex] = []
        current_id: str | None = None
        current_length = 0
        current_offset = 0
        current_line_bases = 0
        current_line_width = 0

        with self.fasta_path.open("rb") as handle:
            while True:
                line_start = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if current_id is not None:
                        rows.append(
                            FastaRecordIndex(
                                current_id,
                                current_length,
                                current_offset,
                                current_line_bases,
                                current_line_width,
                            )
                        )
                    current_id = line[1:].split(None, 1)[0].decode("utf-8")
                    current_length = 0
                    current_offset = handle.tell()
                    current_line_bases = 0
                    current_line_width = 0
                    continue
                if current_id is None:
                    raise ValueError(f"FASTA sequence encountered before header in {self.fasta_path}")
                stripped = line.rstrip(b"\r\n")
                if not stripped:
                    continue
                if current_line_bases == 0:
                    current_line_bases = len(stripped)
                    current_line_width = len(line)
                    current_offset = line_start
                current_length += len(stripped)

        if current_id is not None:
            rows.append(
                FastaRecordIndex(
                    current_id,
                    current_length,
                    current_offset,
                    current_line_bases,
                    current_line_width,
                )
            )

        with self.index_path.open("w", encoding="utf-8") as out:
            out.write("#seq_id\tlength\toffset\tline_bases\tline_width\n")
            for row in rows:
                out.write(
                    f"{row.seq_id}\t{row.length}\t{row.offset}\t"
                    f"{row.line_bases}\t{row.line_width}\n"
                )

    def _read_index(self) -> dict[str, FastaRecordIndex]:
        # Load saved FASTA byte-offset metadata into typed records.
        records: dict[str, FastaRecordIndex] = {}
        with self.index_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                seq_id, length, offset, line_bases, line_width = line.split("\t")
                records[seq_id] = FastaRecordIndex(
                    seq_id=seq_id,
                    length=int(length),
                    offset=int(offset),
                    line_bases=int(line_bases),
                    line_width=int(line_width),
                )
        return records

    def fetch(self, seq_id: str, start: int, end: int) -> str:
        # Fetch a half-open FASTA interval using byte offsets instead of scanning.
        resolved = self.aliases.get(seq_id)
        if resolved is None:
            raise KeyError(f"Sequence {seq_id!r} was not found in {self.fasta_path}")
        record = self.records[resolved]
        if start < 0 or end < start or end > record.length:
            raise ValueError(
                f"Invalid interval {seq_id}:{start}-{end} for sequence length {record.length}"
            )
        remaining = end - start
        position = start
        chunks: list[bytes] = []
        while remaining > 0:
            line_index = position // record.line_bases
            line_offset = position % record.line_bases
            take = min(remaining, record.line_bases - line_offset)
            byte_offset = record.offset + line_index * record.line_width + line_offset
            self.handle.seek(byte_offset)
            chunks.append(self.handle.read(take))
            position += take
            remaining -= take
        return b"".join(chunks).decode("ascii").upper()

    def close(self) -> None:
        # Close the underlying FASTA file handle.
        self.handle.close()


class SequenceComparator:
    """Compare reference and query sequence implied by chain aligned blocks."""

    def __init__(
        self,
        ref_fasta: str | Path,
        query_fasta: str | Path,
        work_dir: str | Path,
    ) -> None:
        # Prepare random-access readers for the reference and query FASTA files.
        work_path = Path(work_dir)
        work_path.mkdir(parents=True, exist_ok=True)
        self.ref = FastaAccessor(ref_fasta, work_path / "ref_fasta.idx")
        self.query = FastaAccessor(query_fasta, work_path / "query_fasta.idx")

    def close(self) -> None:
        # Close both FASTA readers held by the comparator.
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

        try:
            from Bio.Seq import Seq
        except ImportError as exc:
            raise RuntimeError(
                "Biopython is required for sequence comparison. "
                "Activate the biopython conda environment or install biopython."
            ) from exc

        ref_parts: list[str] = []
        query_parts: list[str] = []
        for segment in aligned_segments:
            ref_part = self.ref.fetch(ref_chrom, segment.ref_start, segment.ref_end)
            query_part = self.query.fetch(
                query_chrom, segment.query_start, segment.query_end
            )
            if query_strand == "-":
                query_part = str(Seq(query_part).reverse_complement()).upper()
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
