"""Candidate chain/CNEE pair generation."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .bed import group_bed_by_chrom, parse_bed
from .index import load_index
from .intervals import overlap_length
from .models import BedInterval, ChainRecord

CANDIDATE_COLUMNS = [
    "cnee_uid",
    "cnee_id",
    "ref_chrom",
    "ref_start",
    "ref_end",
    "chain_id",
    "chain_score",
    "chain_ref_start",
    "chain_ref_end",
    "query_chrom",
    "query_start",
    "query_end",
    "query_strand",
    "span_overlap_bases",
]


def _group_chains_by_ref_chrom(
    records: list[ChainRecord],
) -> dict[str, list[ChainRecord]]:
    grouped: defaultdict[str, list[ChainRecord]] = defaultdict(list)
    for record in records:
        grouped[record.ref_chrom].append(record)
    for chrom in grouped:
        grouped[chrom].sort(key=lambda item: (item.ref_start, item.ref_end, item.chain_id))
    return dict(grouped)


def _without_ncbi_version(chrom: str) -> str:
    prefix, dot, suffix = chrom.rpartition(".")
    if dot and suffix.isdigit():
        return prefix
    return chrom


def _bed_chrom_lookup(
    cnees_by_chrom: dict[str, list[BedInterval]],
) -> dict[str, list[BedInterval]]:
    lookup = dict(cnees_by_chrom)
    stripped_to_originals: defaultdict[str, list[str]] = defaultdict(list)
    for chrom in cnees_by_chrom:
        stripped_to_originals[_without_ncbi_version(chrom)].append(chrom)
    for stripped, originals in stripped_to_originals.items():
        if stripped in lookup or len(originals) != 1:
            continue
        lookup[stripped] = cnees_by_chrom[originals[0]]
    return lookup


def iter_candidate_rows(
    chains: list[ChainRecord], cnees: list[BedInterval]
) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    cnees_by_chrom = group_bed_by_chrom(cnees)
    cnees_by_chain_chrom = _bed_chrom_lookup(cnees_by_chrom)
    chains_by_chrom = _group_chains_by_ref_chrom(chains)

    for chrom, chrom_chains in chains_by_chrom.items():
        chrom_cnees = cnees_by_chain_chrom.get(chrom)
        if not chrom_cnees:
            continue
        first_possible = 0
        for chain in chrom_chains:
            while (
                first_possible < len(chrom_cnees)
                and chrom_cnees[first_possible].end <= chain.ref_start
            ):
                first_possible += 1
            idx = first_possible
            while idx < len(chrom_cnees) and chrom_cnees[idx].start < chain.ref_end:
                cnee = chrom_cnees[idx]
                span_overlap = overlap_length(
                    cnee.start, cnee.end, chain.ref_start, chain.ref_end
                )
                if span_overlap > 0:
                    rows.append(
                        {
                            "cnee_uid": cnee.cnee_uid,
                            "cnee_id": cnee.cnee_id,
                            "ref_chrom": cnee.chrom,
                            "ref_start": cnee.start,
                            "ref_end": cnee.end,
                            "chain_id": chain.chain_id,
                            "chain_score": chain.score,
                            "chain_ref_start": chain.ref_start,
                            "chain_ref_end": chain.ref_end,
                            "query_chrom": chain.query_chrom,
                            "query_start": chain.query_start,
                            "query_end": chain.query_end,
                            "query_strand": chain.query_strand,
                            "span_overlap_bases": span_overlap,
                        }
                    )
                idx += 1
    rows.sort(
        key=lambda row: (
            str(row["ref_chrom"]),
            int(row["ref_start"]),
            int(row["ref_end"]),
            str(row["cnee_uid"]),
            -int(row["chain_score"]),
            str(row["chain_id"]),
        )
    )
    return rows


def write_candidates(index_path: str | Path, bed_path: str | Path, out_path: str | Path) -> int:
    _, records = load_index(index_path)
    cnees = parse_bed(bed_path)
    rows = iter_candidate_rows(records, cnees)
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
