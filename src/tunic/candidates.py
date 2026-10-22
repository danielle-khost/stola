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
    # Group chain spans by target/reference chromosome for efficient scanning.
    grouped: defaultdict[str, list[ChainRecord]] = defaultdict(list)
    for record in records:
        grouped[record.ref_chrom].append(record)
    for chrom in grouped:
        grouped[chrom].sort(key=lambda item: (item.ref_start, item.ref_end, item.chain_id))
    return dict(grouped)


def _without_ncbi_version(chrom: str) -> str:
    # Strip a terminal numeric accession version, for mismatch diagnostics only.
    prefix, dot, suffix = chrom.rpartition(".")
    if dot and suffix.isdigit():
        return prefix
    return chrom


def _check_chromosome_names_match(
    chains_by_chrom: dict[str, list[ChainRecord]],
    cnees_by_chrom: dict[str, list[BedInterval]],
) -> None:
    # Chain and BED chromosome names must match exactly -- no alias fallback.
    # A mismatch that would disappear after stripping an NCBI accession
    # version is a common, confusing failure mode, so it gets a specific,
    # actionable error instead of silently dropping those CNEEs as "missing".
    chain_chroms = set(chains_by_chrom)
    bed_chroms = set(cnees_by_chrom)
    if bed_chroms <= chain_chroms:
        return
    stripped_chain_chroms: dict[str, str] = {}
    for chrom in sorted(chain_chroms):
        stripped_chain_chroms.setdefault(_without_ncbi_version(chrom), chrom)

    mismatches = [
        (bed_chrom, stripped_chain_chroms[_without_ncbi_version(bed_chrom)])
        for bed_chrom in sorted(bed_chroms - chain_chroms)
        if _without_ncbi_version(bed_chrom) in stripped_chain_chroms
    ]
    if not mismatches:
        return
    examples = "; ".join(f"BED {bed!r} vs chain {chain!r}" for bed, chain in mismatches[:5])
    more = f", plus {len(mismatches) - 5} more" if len(mismatches) > 5 else ""
    raise ValueError(
        "Reference BED and chain file use different chromosome naming "
        f"conventions -- e.g. {examples}{more}. toga4cnee requires exact "
        "chromosome ID matches between the BED and chain files; rename one "
        "side to match the other (e.g. add or strip the accession version) "
        "and rerun."
    )


def iter_candidate_rows(
    chains: list[ChainRecord], cnees: list[BedInterval]
) -> list[dict[str, str | int]]:
    # Emit every chain/CNEE pair whose reference-side spans overlap at all.
    rows: list[dict[str, str | int]] = []
    cnees_by_chrom = group_bed_by_chrom(cnees)
    chains_by_chrom = _group_chains_by_ref_chrom(chains)
    _check_chromosome_names_match(chains_by_chrom, cnees_by_chrom)

    for chrom, chrom_chains in chains_by_chrom.items():
        chrom_cnees = cnees_by_chrom.get(chrom)
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
    # Load indexed chains and CNEEs, then write the candidate-pair TSV.
    _, records = load_index(index_path)
    cnees = parse_bed(bed_path)
    rows = iter_candidate_rows(records, cnees)
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
