"""Flanking-synteny verification for weak, lost, and uncertain CNEE calls."""

from __future__ import annotations

import csv
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

from .bed import group_bed_by_chrom, parse_bed
from .chain import iter_chain_blocks
from .index import load_index
from .models import BedInterval, ChainBlock, ChainRecord
from .sequence_compare import SequenceComparator, aligned_segments_for_interval
from .summary import stats_row_is_solid

# cnee_uid holds the reference CNEE's cnee_uid for lost/uncertain_loss rows,
# but a query-side ID for weak/one2one rows -- see SyntenyTarget.output_id.
# distance_to_nearest_aligned_block and distance_to_nearest_cnee_neighbor are
# properties of the candidate itself, not the flank check, so they repeat
# identically across a candidate's rows for each flank size.
SYNTENY_COLUMNS = [
    "cnee_uid",
    "orthology_status",
    "flank_size",
    "flank_seq_identity",
    "flank_fraction_gap",
    "distance_to_nearest_aligned_block",
    "distance_to_nearest_cnee_neighbor",
]

# Verify synteny at several scales at once so results across scales can be
# compared without rerunning the workflow.
DEFAULT_FLANK_SIZES: tuple[int, ...] = (10000, 1000, 100)

# Long-format columns for write_synteny_neighbor_verification: one row per
# candidate per flank side, reporting the nearest reference/query CNEE
# neighbor on that side and whether they are orthology partners.
SYNTENY_NEIGHBOR_COLUMNS = [
    "cnee_uid",
    "orthology_status",
    "query_projection_id",
    "ref_coordinates",
    "query_projection_coordinates",
    "flank",
    "nearest_ref_cnee_id",
    "distance_to_ref_cnee",
    "nearest_ref_cnee_status",
    "nearest_query_cnee_id",
    "distance_to_query_cnee",
    "nearest_query_cnee_status",
    "synteny_matches",
]


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


def write_synteny_verification(
    index_path: str | Path,
    stats_path: str | Path,
    query_bed_path: str | Path,
    orthology_path: str | Path,
    cnee_bed_path: str | Path,
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
    Each row also reports the candidate's distance to the nearest aligned
    chain block (0 if the candidate already overlaps one) and its distance
    to the nearest neighboring reference CNEE annotation.
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
    neighbor_distances = _nearest_cnee_neighbor_distances(cnee_bed_path)
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
                blocks = block_cache[chain_id]
                block_distance = _distance_to_nearest_block(
                    blocks, int(target.row["ref_start"]), int(target.row["ref_end"])
                )
                neighbor_distance = neighbor_distances.get(target.row["cnee_uid"])
                for flank_size in flank_sizes:
                    flank_metrics = _flank_metrics(
                        target.row,
                        record,
                        blocks,
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
                            "distance_to_nearest_aligned_block": block_distance,
                            "distance_to_nearest_cnee_neighbor": (
                                "." if neighbor_distance is None else neighbor_distance
                            ),
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


def _distance_to_nearest_block(blocks: BlockIndex, start: int, end: int) -> int:
    # 0 if the candidate already sits inside an aligned block (the common
    # case for anything but "lost"); otherwise the gap to the nearest block
    # in the same chain, on whichever side is closer.
    for _ in blocks.overlapping(start, end):
        return 0

    preceding_index = bisect_right(blocks.ref_ends, start) - 1
    following_index = bisect_left(blocks.ref_starts, end)
    distances = []
    if preceding_index >= 0:
        distances.append(start - blocks.blocks[preceding_index].ref_end)
    if following_index < len(blocks.blocks):
        distances.append(blocks.blocks[following_index].ref_start - end)
    return min(distances) if distances else -1


def _nearest_cnee_neighbor_distances(cnee_bed_path: str | Path) -> dict[str, int | None]:
    # Precompute, for every reference CNEE, the distance to its nearest
    # neighboring CNEE annotation on the same chromosome (independent of any
    # chain alignment). None means no other CNEE shares that chromosome.
    by_chrom = group_bed_by_chrom(parse_bed(cnee_bed_path))
    distances: dict[str, int | None] = {}
    for chrom_intervals in by_chrom.values():
        last_index = len(chrom_intervals) - 1
        for index, interval in enumerate(chrom_intervals):
            candidates = []
            if index > 0:
                candidates.append(max(0, interval.start - chrom_intervals[index - 1].end))
            if index < last_index:
                candidates.append(max(0, chrom_intervals[index + 1].start - interval.end))
            distances[interval.cnee_uid] = min(candidates) if candidates else None
    return distances


@dataclass(frozen=True)
class RefNeighborIndex:
    """Directional nearest-neighbor lookup over the reference CNEE annotation.

    Assumes reference CNEEs on a given chromosome do not overlap, same as
    every other reference-BED consumer in this module.
    """

    intervals_by_chrom: dict[str, list[BedInterval]]
    chrom_by_uid: dict[str, str]
    position_by_uid: dict[str, int]

    @classmethod
    def from_bed(cls, cnee_bed_path: str | Path) -> "RefNeighborIndex":
        by_chrom = group_bed_by_chrom(parse_bed(cnee_bed_path))
        chrom_by_uid: dict[str, str] = {}
        position_by_uid: dict[str, int] = {}
        for chrom, intervals in by_chrom.items():
            for position, interval in enumerate(intervals):
                chrom_by_uid[interval.cnee_uid] = chrom
                position_by_uid[interval.cnee_uid] = position
        return cls(by_chrom, chrom_by_uid, position_by_uid)

    def neighbor(self, cnee_uid: str, side: str) -> tuple[str | None, int | None]:
        # Return the immediately preceding ("left") or following ("right")
        # reference CNEE and the gap to it, or (None, None) at a chromosome edge.
        chrom = self.chrom_by_uid.get(cnee_uid)
        if chrom is None:
            return None, None
        intervals = self.intervals_by_chrom[chrom]
        position = self.position_by_uid[cnee_uid]
        current = intervals[position]
        if side == "left":
            if position == 0:
                return None, None
            neighbor = intervals[position - 1]
            return neighbor.cnee_uid, max(0, current.start - neighbor.end)
        if position == len(intervals) - 1:
            return None, None
        neighbor = intervals[position + 1]
        return neighbor.cnee_uid, max(0, neighbor.start - current.end)


@dataclass(frozen=True)
class QueryNeighborIndex:
    """Directional nearest-neighbor lookup over the all-mode query CNEE BED.

    Assumes query CNEEs on a given chromosome do not overlap (this can, in
    principle, be violated by clusters on opposite strands sharing a locus --
    the same assumption RefNeighborIndex makes for the reference BED).
    """

    by_chrom: dict[str, list[QueryBedRecord]]
    starts_by_chrom: dict[str, list[int]]
    ends_by_chrom: dict[str, list[int]]

    @classmethod
    def from_records(
        cls, by_chrom: dict[str, list[QueryBedRecord]]
    ) -> "QueryNeighborIndex":
        starts = {chrom: [r.start for r in records] for chrom, records in by_chrom.items()}
        ends = {chrom: [r.end for r in records] for chrom, records in by_chrom.items()}
        return cls(by_chrom, starts, ends)

    def neighbor(
        self, chrom: str, position: int, direction: str
    ) -> tuple[str | None, int | None]:
        # direction "before": nearest record ending at or before position.
        # direction "after": nearest record starting at or after position.
        records = self.by_chrom.get(chrom)
        if not records:
            return None, None
        if direction == "before":
            index = bisect_right(self.ends_by_chrom[chrom], position) - 1
            if index < 0:
                return None, None
            record = records[index]
            return record.query_id, max(0, position - record.end)
        index = bisect_left(self.starts_by_chrom[chrom], position)
        if index >= len(records):
            return None, None
        record = records[index]
        return record.query_id, max(0, record.start - position)


@dataclass(frozen=True)
class OrthologyIndex:
    """Reverse lookups over the all-mode orthology TSV, by ref and by query CNEE."""

    category_by_ref_uid: dict[str, str]
    category_by_query_id: dict[str, str]
    ref_partners_by_query_id: dict[str, frozenset[str]]

    @classmethod
    def from_path(cls, orthology_path: str | Path) -> "OrthologyIndex":
        category_by_ref_uid: dict[str, str] = {}
        category_by_query_id: dict[str, str] = {}
        ref_partners_by_query_id: dict[str, frozenset[str]] = {}
        path = Path(orthology_path)
        if not path.exists():
            return cls(category_by_ref_uid, category_by_query_id, ref_partners_by_query_id)
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"query_cnee_id", "partner_ref_cnees", "orthology_category"}
            missing = sorted(required.difference(reader.fieldnames or []))
            if missing:
                raise ValueError(f"Orthology file is missing columns: {', '.join(missing)}")
            for row in reader:
                category = row["orthology_category"]
                ref_ids = [ref_id for ref_id in row["partner_ref_cnees"].split(",") if ref_id]
                for ref_id in ref_ids:
                    category_by_ref_uid[ref_id] = category
                query_id = row["query_cnee_id"]
                if query_id != "NA":
                    category_by_query_id[query_id] = category
                    ref_partners_by_query_id[query_id] = frozenset(ref_ids)
        return cls(category_by_ref_uid, category_by_query_id, ref_partners_by_query_id)


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


def _query_flank_anchors(
    row: dict[str, str], record: ChainRecord, blocks: BlockIndex
) -> dict[str, tuple[int | None, str | None]]:
    """Compute, per flank side, the query coordinate + search direction used
    to anchor the nearest-query-CNEE lookup.

    Reference and query chain coordinates are both normalized to plus-strand
    genome coordinates (see chain.py), so walking forward along the reference
    walks *backward* along the query whenever record.query_strand == "-" --
    which side of the reference interval corresponds to which side of the
    query anchor flips accordingly.
    """
    forward = record.query_strand == "+"
    query_start = row["query_aligned_start"]
    query_end = row["query_aligned_end"]
    if query_start != "." and query_end != ".":
        # Candidate has its own aligned query interval (weak/uncertain_loss/
        # one2one): anchor directly on its edges.
        start = int(query_start)
        end = int(query_end)
        if forward:
            return {"left": (start, "before"), "right": (end, "after")}
        return {"left": (end, "after"), "right": (start, "before")}

    # "Lost" candidate: no aligned query bases at all, so interpolate the
    # query anchor from the chain blocks immediately flanking its gap.
    ref_start = int(row["ref_start"])
    ref_end = int(row["ref_end"])
    preceding, following = _flanking_blocks(blocks, ref_start, ref_end)
    if forward:
        left = (preceding.query_end, "before") if preceding else (None, None)
        right = (following.query_start, "after") if following else (None, None)
    else:
        left = (preceding.query_start, "after") if preceding else (None, None)
        right = (following.query_end, "before") if following else (None, None)
    return {"left": left, "right": right}


def _coordinates(chrom: str, start: str | int, end: str | int) -> str:
    # Format a genomic interval as chr:start-end.
    return f"{chrom}:{start}-{end}"


def _candidate_query_projection(target: "SyntenyTarget") -> tuple[str, str]:
    """Return (query_projection_id, query_projection_coordinates) for a
    candidate's own query-side projection -- not a neighbor.

    "lost" candidates have no aligned query bases at all, so both are "NA".
    "uncertain_loss" candidates do have an aligned (non-solid) query
    interval, but it was never graph-assigned a query CNEE id, so the raw
    per-alignment projection name is reported instead of a formal id.
    "weak"/"one2one" candidates report the graph-assigned query CNEE id
    (already resolved onto SyntenyTarget.output_id).
    """
    row = target.row
    query_start = row["query_aligned_start"]
    query_end = row["query_aligned_end"]
    if query_start == "." or query_end == ".":
        return "NA", "NA"
    coordinates = _coordinates(row["query_chrom"], query_start, query_end)
    if target.status in ("weak", "one2one"):
        return target.output_id, coordinates
    return row["query_projection_name"], coordinates


def write_synteny_neighbor_verification(
    index_path: str | Path,
    stats_path: str | Path,
    query_bed_path: str | Path,
    orthology_path: str | Path,
    cnee_bed_path: str | Path,
    out_path: str | Path,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    max_one2one_controls: int = 50000,
) -> SyntenyResult:
    """Write a left/right flanking-CNEE synteny check for weak/lost/uncertain/
    one2one-control candidates.

    For each candidate and each flank side, report the nearest neighboring
    reference CNEE and the nearest neighboring query CNEE -- projected from
    the candidate's own aligned query interval, or (for "lost" candidates
    with no aligned query bases at all) interpolated from the chain blocks
    immediately flanking the candidate's chain gap -- plus a synteny_matches
    flag: true iff the query neighbor is one of the reference neighbor's
    orthology partners.
    """
    rows = _read_stats(stats_path)
    targets = _select_synteny_targets(
        rows,
        query_bed_path,
        orthology_path,
        min_sequence_identity=min_sequence_identity,
        max_fraction_gap=max_fraction_gap,
        max_one2one_controls=max_one2one_controls,
    )
    ref_index = RefNeighborIndex.from_bed(cnee_bed_path)
    query_index = QueryNeighborIndex.from_records(_read_query_bed(query_bed_path))
    orthology_index = OrthologyIndex.from_path(orthology_path)

    payload, records = load_index(index_path)
    chain_path = payload["indexed_chain_path"]
    record_by_id = {record.chain_id: record for record in records}
    block_cache: dict[str, BlockIndex] = {}
    count = 0

    with Path(chain_path).open("rb") as chain_handle:
        with Path(out_path).open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=SYNTENY_NEIGHBOR_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for target in targets:
                cnee_uid = target.row["cnee_uid"]
                chain_id = target.row["chain_id"]
                record = record_by_id[chain_id]
                if chain_id not in block_cache:
                    chain_handle.seek(record.start_byte)
                    chain_text = chain_handle.read(record.byte_length).decode("utf-8")
                    block_cache[chain_id] = BlockIndex.from_blocks(
                        iter_chain_blocks(chain_text)
                    )
                blocks = block_cache[chain_id]
                anchors = _query_flank_anchors(target.row, record, blocks)
                ref_coordinates = _coordinates(
                    target.row["ref_chrom"],
                    target.row["ref_start"],
                    target.row["ref_end"],
                )
                query_projection_id, query_projection_coordinates = (
                    _candidate_query_projection(target)
                )
                for side in ("left", "right"):
                    ref_neighbor_id, ref_distance = ref_index.neighbor(cnee_uid, side)
                    ref_neighbor_status = (
                        orthology_index.category_by_ref_uid.get(ref_neighbor_id, ".")
                        if ref_neighbor_id is not None
                        else "."
                    )
                    anchor, direction = anchors[side]
                    if anchor is None:
                        query_neighbor_id: str | None = None
                        query_distance: int | None = None
                    else:
                        query_neighbor_id, query_distance = query_index.neighbor(
                            record.query_chrom, anchor, direction
                        )
                    query_neighbor_status = (
                        orthology_index.category_by_query_id.get(query_neighbor_id, ".")
                        if query_neighbor_id is not None
                        else "."
                    )
                    if ref_neighbor_id is None or query_neighbor_id is None:
                        synteny_matches = "."
                    else:
                        synteny_matches = str(
                            ref_neighbor_id
                            in orthology_index.ref_partners_by_query_id.get(
                                query_neighbor_id, frozenset()
                            )
                        ).upper()
                    writer.writerow(
                        {
                            "cnee_uid": cnee_uid,
                            "orthology_status": target.status,
                            "query_projection_id": query_projection_id,
                            "ref_coordinates": ref_coordinates,
                            "query_projection_coordinates": query_projection_coordinates,
                            "flank": side,
                            "nearest_ref_cnee_id": ref_neighbor_id or ".",
                            "distance_to_ref_cnee": (
                                "." if ref_distance is None else ref_distance
                            ),
                            "nearest_ref_cnee_status": ref_neighbor_status,
                            "nearest_query_cnee_id": query_neighbor_id or ".",
                            "distance_to_query_cnee": (
                                "." if query_distance is None else query_distance
                            ),
                            "nearest_query_cnee_status": query_neighbor_status,
                            "synteny_matches": synteny_matches,
                        }
                    )
                    count += 1
    return SyntenyResult(rows=count)
