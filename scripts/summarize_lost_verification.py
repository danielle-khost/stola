#!/usr/bin/env python
"""Exploratory summary of "lost" CNEEs, per tests/tunic_lost_verification.png.

For every reference CNEE with candidate chain(s) but zero aligned query bases
(the "lost" category), looks at the chain gap it sits inside -- bounded by
the nearest aligned block on each side, within that CNEE's own best-ranked
candidate chain -- and reports:

  - the reference-side and query-side length of that gap
  - a loss_type classification from their ratio:
      DEL  query gap length < 0.9x reference gap length
      INS  query gap length > 1.1x reference gap length
      RAR  query gap length within +/-10% of the reference gap length
      UNK  the gap is "ambiguous" -- no aligned block on one side at all
           (the CNEE sits at the edge of its chain's aligned span), so
           sizes/loss_type/everything else are unavailable
  - how many ambiguous bases (Ns) the query gap sequence contains
  - for RAR only: whether that same query gap region is covered by an
    aligned block belonging to a DIFFERENT chain (a coarse whole-chain
    bounding-box filter over every chain record, then an exact per-block
    check only for the chains that pass the coarse filter) -- if so, the
    reference coordinates of that alternate alignment are reported,
    otherwise "NA". This only finds a match if the genome aligner already
    chained the query sequence somewhere else; it is not a new sequence
    search.

This is standalone exploration -- it does not modify the main orthology
pipeline or its "lost" category, which v2 verification still leaves alone.

Example:
    PYTHONPATH=src python scripts/summarize_lost_verification.py \\
        --stats tests/test_runs/colLivi/overlap_stats.tsv \\
        --chain-index tests/test_runs/colLivi/work/chain_index.json \\
        --query-fasta tests/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna \\
        --out tests/test_runs/colLivi/lost_verification_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from tunic.chain import iter_chain_blocks
from tunic.index import load_index
from tunic.models import ChainBlock, ChainRecord
from tunic.genome_io import FastaAccessor
from tunic.synteny import BlockIndex, _flanking_blocks

OUT_COLUMNS = [
    "ref_cnee_id",
    "ref_coordinates",
    "cnee_length",
    "ref_gap_length",
    "query_gap_length",
    "loss_type",
    "num_ambiguous_bases",
    "alt_align",
]

RAR_TOLERANCE = 0.1


@dataclass(frozen=True)
class LostCandidate:
    ref_cnee_id: str
    ref_chrom: str
    ref_start: int
    ref_end: int
    chain_id: str


def _select_lost_candidates(stats_path: str | Path) -> list[LostCandidate]:
    # A ref CNEE is "lost" if it has candidate chain(s) but none give any
    # aligned query bases at all -- mirrors classify_orthology's own "lost"
    # definition. Its best-ranked (rank_by_chain_score == 1) candidate row
    # is the chain whose gap gets analyzed.
    best_row_by_uid: dict[str, dict[str, str]] = {}
    any_aligned_by_uid: set[str] = set()
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cnee_uid",
            "ref_chrom",
            "ref_start",
            "ref_end",
            "chain_id",
            "aligned_bases",
            "rank_by_chain_score",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            uid = row["cnee_uid"]
            if int(row["aligned_bases"]) > 0:
                any_aligned_by_uid.add(uid)
                continue
            if int(row["rank_by_chain_score"]) == 1:
                best_row_by_uid[uid] = row

    return [
        LostCandidate(
            ref_cnee_id=uid,
            ref_chrom=row["ref_chrom"],
            ref_start=int(row["ref_start"]),
            ref_end=int(row["ref_end"]),
            chain_id=row["chain_id"],
        )
        for uid, row in best_row_by_uid.items()
        if uid not in any_aligned_by_uid
    ]


class QueryBoundingBoxIndex:
    """Coarse whole-chain query-side bounding-box lookup, for the alt-align
    check -- built once from every ChainRecord's own (already plus-strand
    normalized) query_start/query_end, no per-block parsing required.
    """

    def __init__(self, records: list[ChainRecord]) -> None:
        by_chrom: dict[str, list[ChainRecord]] = defaultdict(list)
        for record in records:
            by_chrom[record.query_chrom].append(record)
        self.by_chrom: dict[str, list[ChainRecord]] = {}
        self.ends_by_chrom: dict[str, list[int]] = {}
        for chrom, chrom_records in by_chrom.items():
            ordered = sorted(chrom_records, key=lambda r: (r.query_start, r.query_end))
            self.by_chrom[chrom] = ordered
            self.ends_by_chrom[chrom] = [r.query_end for r in ordered]

    def overlapping(self, chrom: str, start: int, end: int) -> list[ChainRecord]:
        records = self.by_chrom.get(chrom)
        if not records:
            return []
        ends = self.ends_by_chrom[chrom]
        index = bisect_right(ends, start)
        found = []
        while index < len(records):
            record = records[index]
            if record.query_start >= end:
                break
            if record.query_end > start:
                found.append(record)
            index += 1
        return found


def _chain_text(chain_handle, record: ChainRecord) -> str:
    chain_handle.seek(record.start_byte)
    return chain_handle.read(record.byte_length).decode("utf-8")


def _project_query_overlap_to_ref(
    block: ChainBlock, query_strand: str, query_overlap_start: int, query_overlap_end: int
) -> tuple[int, int]:
    # Inverse of overlap.py's _project_overlap_to_query: given a query-side
    # overlap, find the corresponding reference-side sub-interval.
    if query_strand == "+":
        offset_start = query_overlap_start - block.query_start
        offset_end = query_overlap_end - block.query_start
        ref_start = block.ref_start + offset_start
        ref_end = block.ref_start + offset_end
    else:
        offset_start = block.query_end - query_overlap_end
        offset_end = block.query_end - query_overlap_start
        ref_start = block.ref_start + offset_start
        ref_end = block.ref_start + offset_end
    return min(ref_start, ref_end), max(ref_start, ref_end)


def _find_alt_alignment(
    query_chrom: str,
    query_gap_start: int,
    query_gap_end: int,
    exclude_chain_id: str,
    bbox_index: QueryBoundingBoxIndex,
    record_by_chain_id: dict[str, ChainRecord],
    chain_handle,
    block_cache: dict[str, BlockIndex],
) -> str:
    # Coarse pass: other chains whose whole-chain query bbox overlaps at all.
    coarse_hits = [
        record
        for record in bbox_index.overlapping(query_chrom, query_gap_start, query_gap_end)
        if record.chain_id != exclude_chain_id
    ]
    if not coarse_hits:
        return "NA"
    coarse_hits.sort(key=lambda r: (-r.score, r.chain_id))

    for record in coarse_hits:
        if record.chain_id not in block_cache:
            blocks = iter_chain_blocks(_chain_text(chain_handle, record))
            block_cache[record.chain_id] = BlockIndex.from_blocks(blocks)
        blocks = block_cache[record.chain_id]
        for block in blocks.blocks:
            overlap_start = max(query_gap_start, block.query_start)
            overlap_end = min(query_gap_end, block.query_end)
            if overlap_end <= overlap_start:
                continue
            ref_start, ref_end = _project_query_overlap_to_ref(
                block, record.query_strand, overlap_start, overlap_end
            )
            return f"{record.ref_chrom}:{ref_start}-{ref_end}"
    return "NA"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True, type=Path, help="overlap_stats.tsv (all-mode, raw)")
    parser.add_argument("--chain-index", required=True, type=Path)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument(
        "--work-dir", type=Path, default=None, help="FASTA index cache dir [default: alongside --chain-index]"
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)

    print("Selecting lost CNEEs from overlap stats...", file=sys.stderr)
    lost_candidates = _select_lost_candidates(args.stats)
    print(f"  {len(lost_candidates)} lost CNEEs", file=sys.stderr)

    payload, records = load_index(args.chain_index)
    record_by_chain_id = {record.chain_id: record for record in records}
    chain_path = payload["indexed_chain_path"]
    bbox_index = QueryBoundingBoxIndex(records)

    work_dir = args.work_dir if args.work_dir is not None else args.chain_index.parent
    query_fasta = FastaAccessor(args.query_fasta, work_dir / "query_fasta.idx")

    block_cache: dict[str, BlockIndex] = {}
    rows: list[dict[str, str | int]] = []

    try:
        with Path(chain_path).open("rb") as chain_handle:
            for done, candidate in enumerate(lost_candidates, start=1):
                record = record_by_chain_id[candidate.chain_id]
                if candidate.chain_id not in block_cache:
                    block_cache[candidate.chain_id] = BlockIndex.from_blocks(
                        iter_chain_blocks(_chain_text(chain_handle, record))
                    )
                blocks = block_cache[candidate.chain_id]
                cnee_length = candidate.ref_end - candidate.ref_start
                ref_coordinates = f"{candidate.ref_chrom}:{candidate.ref_start}-{candidate.ref_end}"

                preceding, following = _flanking_blocks(blocks, candidate.ref_start, candidate.ref_end)
                if preceding is None or following is None:
                    rows.append(
                        {
                            "ref_cnee_id": candidate.ref_cnee_id,
                            "ref_coordinates": ref_coordinates,
                            "cnee_length": cnee_length,
                            "ref_gap_length": "NA",
                            "query_gap_length": "NA",
                            "loss_type": "UNK",
                            "num_ambiguous_bases": "NA",
                            "alt_align": "NA",
                        }
                    )
                    continue

                ref_gap_length = following.ref_start - preceding.ref_end
                if record.query_strand == "+":
                    query_gap_start, query_gap_end = preceding.query_end, following.query_start
                else:
                    query_gap_start, query_gap_end = following.query_end, preceding.query_start
                query_gap_length = max(0, query_gap_end - query_gap_start)

                if query_gap_length < (1 - RAR_TOLERANCE) * ref_gap_length:
                    loss_type = "DEL"
                elif query_gap_length > (1 + RAR_TOLERANCE) * ref_gap_length:
                    loss_type = "INS"
                else:
                    loss_type = "RAR"

                num_ambiguous_bases = 0
                if query_gap_length > 0:
                    query_seq = query_fasta.fetch(record.query_chrom, query_gap_start, query_gap_end)
                    num_ambiguous_bases = query_seq.count("N")

                alt_align = "NA"
                if loss_type == "RAR" and query_gap_length > 0:
                    alt_align = _find_alt_alignment(
                        record.query_chrom,
                        query_gap_start,
                        query_gap_end,
                        candidate.chain_id,
                        bbox_index,
                        record_by_chain_id,
                        chain_handle,
                        block_cache,
                    )

                rows.append(
                    {
                        "ref_cnee_id": candidate.ref_cnee_id,
                        "ref_coordinates": ref_coordinates,
                        "cnee_length": cnee_length,
                        "ref_gap_length": ref_gap_length,
                        "query_gap_length": query_gap_length,
                        "loss_type": loss_type,
                        "num_ambiguous_bases": num_ambiguous_bases,
                        "alt_align": alt_align,
                    }
                )
                if done % 2000 == 0:
                    print(f"  processed {done}/{len(lost_candidates)}...", file=sys.stderr)
    finally:
        query_fasta.close()

    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    loss_type_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        loss_type_counts[row["loss_type"]] += 1
    print(f"\nWrote {len(rows)} rows -> {args.out}", file=sys.stderr)
    for loss_type in ("DEL", "INS", "RAR", "UNK"):
        if loss_type in loss_type_counts:
            print(f"  {loss_type}: {loss_type_counts[loss_type]}", file=sys.stderr)


if __name__ == "__main__":
    main()
