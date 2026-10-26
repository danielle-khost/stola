"""Genome sequence access shared by FASTA- and 2bit-backed inputs.

Format is detected by content (a .2bit file's 4-byte magic signature), not
by filename extension. Everything downstream of this module -- sequence
identity comparison in sequence_compare.py, contig-boundary N-run inference
in assembly_gaps.py -- works through the same two entry points regardless of
which format a given ref/query genome file turns out to be:

  - open_genome_accessor(path, index_path): random-access point lookups
    (fetch a short window by coordinate), for sequence comparison.
  - iter_named_sequences(path): full-sequence iteration, one record at a
    time, for whole-genome scans.

ref and query are sniffed independently, so a FASTA ref paired with a 2bit
query (or vice versa) just works.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

TWOBIT_SIGNATURE = 0x1A412743


def sniff_genome_format(path: str | Path) -> str:
    """Return "twobit" or "fasta" by peeking the file's first 4 bytes."""
    with Path(path).open("rb") as handle:
        data = handle.read(4)
    if len(data) == 4:
        for byteorder in ("little", "big"):
            if int.from_bytes(data, byteorder) == TWOBIT_SIGNATURE:
                return "twobit"
    return "fasta"


def check_twobit_complete(path: str | Path) -> None:
    """Raise ValueError if a .2bit file is shorter than its own index says.

    Biopython's TwoBitIterator doesn't notice truncation: a record whose
    offset lies past end-of-file decodes from empty reads as a zero-length
    sequence, so a partially downloaded .2bit (which keeps its full name
    index, since that's at the front of the file) opens without complaint
    and only fails much later, as an "Invalid interval ... for sequence
    length 0" from whichever fetch() first lands on a missing chromosome.
    This walks the index and checks that every record's header and packed
    bases fit inside the file.
    """
    path = Path(path)
    size = path.stat().st_size

    def _truncated(detail: str) -> ValueError:
        return ValueError(
            f"{path} is truncated ({size:,} bytes): {detail}. This usually "
            "means an interrupted download -- re-download the file (`wget -c "
            "URL` resumes a partial one) and rerun."
        )

    with path.open("rb") as handle:

        def _read(count: int) -> bytes:
            data = handle.read(count)
            if len(data) < count:
                raise _truncated("the file ends inside its sequence-name index")
            return data

        header = _read(16)
        for order in ("<", ">"):
            signature, version, seq_count, _ = struct.unpack(f"{order}IIII", header)
            if signature == TWOBIT_SIGNATURE:
                break
        else:
            return  # Not a .2bit; TwoBitIterator reports the bad signature.
        if version not in (0, 1):
            return  # TwoBitIterator rejects unknown versions with its own message.
        offset_format = f"{order}Q" if version == 1 else f"{order}I"
        offset_size = struct.calcsize(offset_format)

        index: list[tuple[str, int]] = []
        for _ in range(seq_count):
            name = _read(_read(1)[0]).decode("ascii", errors="replace")
            (offset,) = struct.unpack(offset_format, _read(offset_size))
            index.append((name, offset))

        incomplete: list[str] = []
        for name, offset in index:
            # Record layout: dnaSize, nBlockCount, nBlock starts+sizes,
            # maskBlockCount, mask starts+sizes, reserved, then 4 bases/byte.
            handle.seek(offset)
            fixed = handle.read(8)
            if len(fixed) < 8:
                incomplete.append(name)
                continue
            dna_size, n_block_count = struct.unpack(f"{order}II", fixed)
            handle.seek(offset + 8 + 8 * n_block_count)
            mask_count_bytes = handle.read(4)
            if len(mask_count_bytes) < 4:
                incomplete.append(name)
                continue
            (mask_block_count,) = struct.unpack(f"{order}I", mask_count_bytes)
            record_end = (
                offset + 8 + 8 * n_block_count + 4 + 8 * mask_block_count + 4
                + (dna_size + 3) // 4
            )
            if record_end > size:
                incomplete.append(name)

    if incomplete:
        raise _truncated(
            f"{len(incomplete)} of {seq_count} sequences extend past the end "
            f"of the file, starting with {incomplete[0]!r}"
        )


def _strip_version(seq_id: str) -> str:
    # Remove a terminal numeric accession version from a sequence ID.
    prefix, dot, suffix = seq_id.rpartition(".")
    if dot and suffix.isdigit():
        return prefix
    return seq_id


def _build_aliases(keys: Iterable[str]) -> dict[str, str]:
    # Map exact IDs and unambiguous version-stripped IDs to genome records.
    exact = {key: key for key in keys}
    stripped: dict[str, list[str]] = {}
    for key in exact:
        stripped.setdefault(_strip_version(key), []).append(key)
    aliases = dict(exact)
    for key, values in stripped.items():
        if key in aliases or len(values) != 1:
            continue
        aliases[key] = values[0]
    return aliases


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
        self.aliases = _build_aliases(self.records.keys())
        self.handle = self.fasta_path.open("rb")

    def _build_index(self) -> None:
        # Scan the FASTA once and save per-record byte layout metadata.
        self.index_path.parent.mkdir(parents=True, exist_ok=True)
        rows: list[FastaRecordIndex] = []
        current_id: str | None = None
        current_length = 0
        current_offset = 0
        current_line_bases = 0
        current_line_width = 0
        # Byte-offset lookup in fetch() assumes every line of a record except
        # the last has the same width. A line's own length can't be checked
        # until we know it isn't the record's last line, so the check on line
        # N is deferred until line N+1 (or the next header/EOF) arrives.
        previous_line_bases: int | None = None
        previous_line_number = 0
        line_number = 0

        with self.fasta_path.open("rb") as handle:
            while True:
                line_start = handle.tell()
                line = handle.readline()
                line_number += 1
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
                    previous_line_bases = None
                    continue
                if current_id is None:
                    raise ValueError(f"FASTA sequence encountered before header in {self.fasta_path}")
                stripped = line.rstrip(b"\r\n")
                if not stripped:
                    continue
                if previous_line_bases is not None and previous_line_bases != current_line_bases:
                    # Fail loudly rather than silently reflow -- an incorrect
                    # auto-fix here is worse than an error the user can act on.
                    raise ValueError(
                        f"{self.fasta_path}: sequence {current_id!r} has ragged "
                        f"line wrapping -- line {previous_line_number} has "
                        f"{previous_line_bases} bases, but the record's other "
                        f"lines have {current_line_bases}. stola requires "
                        "every line of a FASTA record except the last to be the "
                        "same length. Re-wrap this file to a uniform width (e.g. "
                        "`seqkit seq -w 60 in.fa > out.fa`) and rerun."
                    )
                if current_line_bases == 0:
                    current_line_bases = len(stripped)
                    current_line_width = len(line)
                    current_offset = line_start
                current_length += len(stripped)
                previous_line_bases = len(stripped)
                previous_line_number = line_number

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


class TwoBitAccessor:
    """Random-access .2bit reader, matching FastaAccessor's fetch/close contract.

    Unlike FASTA, .2bit is natively indexed (the format's own header stores
    per-sequence offsets and N-run/mask-run tables), so there's no sidecar
    index file to build or cache -- opening the file is enough.

    Biopython's own per-call lookup (Bio.Seq's lazy .2bit sequence data)
    calls numpy.fromfile fresh on every single fetch(), which profiles at
    ~35-40x slower per call than FastaAccessor's plain seek()/read() even in
    the best case, and far worse under the scattered access pattern this
    pipeline actually produces (confirmed: a real ~600k-chain run went from
    ~2 minutes with FASTA to ~4 hours of CPU time with .2bit). Fully
    decoding a whole chromosome at once, by contrast, runs at ~1.1-1.3 GB/s
    with no sign of superlinear cost at any size tested (a 196Mbp chromosome
    decodes in ~0.15s) -- so each sequence is decoded to a plain str at most
    once, on first access, and every fetch after that is a cheap string
    slice. Worth the memory (up to one full genome's worth of plain-str
    sequence, decoded lazily per chromosome so untouched scaffolds cost
    nothing) for the several-orders-of-magnitude speedup.
    """

    def __init__(self, twobit_path: str | Path) -> None:
        try:
            from Bio.SeqIO.TwoBitIO import TwoBitIterator
        except ImportError as exc:
            raise RuntimeError(
                "Biopython is required to read .2bit genome files. "
                "Activate the biopython conda environment or install biopython."
            ) from exc

        self.twobit_path = Path(twobit_path)
        check_twobit_complete(self.twobit_path)
        self._iterator = TwoBitIterator(str(self.twobit_path))
        self.sequences = self._iterator.sequences
        self.aliases = _build_aliases(self.sequences.keys())
        self._decoded: dict[str, str] = {}

    def _decoded_sequence(self, resolved: str) -> str:
        text = self._decoded.get(resolved)
        if text is None:
            text = str(self.sequences[resolved]).upper()
            self._decoded[resolved] = text
        return text

    def fetch(self, seq_id: str, start: int, end: int) -> str:
        resolved = self.aliases.get(seq_id)
        if resolved is None:
            raise KeyError(f"Sequence {seq_id!r} was not found in {self.twobit_path}")
        text = self._decoded_sequence(resolved)
        if start < 0 or end < start or end > len(text):
            raise ValueError(
                f"Invalid interval {seq_id}:{start}-{end} for sequence length {len(text)}"
            )
        return text[start:end]

    def close(self) -> None:
        self._decoded.clear()
        self._iterator.stream.close()


def open_genome_accessor(
    path: str | Path, index_path: str | Path
) -> FastaAccessor | TwoBitAccessor:
    """Open a genome file for random-access fetch(seq_id, start, end).

    index_path is only used (and only created) for a FASTA input -- a .2bit
    input ignores it entirely, since the format needs no sidecar index.
    """
    if sniff_genome_format(path) == "twobit":
        return TwoBitAccessor(path)
    return FastaAccessor(path, index_path)


def iter_named_sequences(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (name, full sequence as an uppercase str) for every record.

    Dispatches on sniffed format. Each record's full sequence is decoded at
    once (never more than one record in memory at a time), for whole-genome
    scans that need every base -- e.g. finding N-runs.
    """
    if sniff_genome_format(path) == "twobit":
        yield from _iter_twobit_sequences(path)
    else:
        yield from _iter_fasta_sequences(path)


def _iter_fasta_sequences(path: str | Path) -> Iterator[tuple[str, str]]:
    chrom: str | None = None
    buffer: list[bytes] = []

    def _flush() -> Iterator[tuple[str, str]]:
        if chrom is None or not buffer:
            return
        yield chrom, b"".join(buffer).decode("ascii").upper()

    with Path(path).open("rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b">"):
                yield from _flush()
                chrom = raw_line[1:].split(None, 1)[0].decode("utf-8")
                buffer = []
                continue
            buffer.append(raw_line.rstrip(b"\r\n"))
        yield from _flush()


def _iter_twobit_sequences(path: str | Path) -> Iterator[tuple[str, str]]:
    try:
        from Bio.SeqIO.TwoBitIO import TwoBitIterator
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required to read .2bit genome files. "
            "Activate the biopython conda environment or install biopython."
        ) from exc

    check_twobit_complete(path)
    iterator = TwoBitIterator(str(path))
    try:
        for name, seq in iterator.sequences.items():
            yield name, str(seq).upper()
    finally:
        iterator.stream.close()
