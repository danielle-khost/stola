#!/usr/bin/env python3
"""Make a small fake UCSC chain file from simuG-style map coordinates."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GapEvent:
    ref_gap_start: int
    target_gap: int
    query_gap: int
    variant_id: str

    @property
    def ref_gap_end(self) -> int:
        # Return the target/reference coordinate immediately after this gap.
        return self.ref_gap_start + self.target_gap


@dataclass(frozen=True)
class ChainSpec:
    chain_id: int
    target_start: int
    target_end: int


@dataclass(frozen=True)
class ChainBuild:
    text: str
    used_events: int
    skipped_events: int


def read_fasta_size(path: str | Path) -> tuple[str, int]:
    # Read the first FASTA record name and total sequence length.
    chrom = ""
    size = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if chrom:
                    break
                chrom = line[1:].split()[0]
                continue
            if chrom:
                size += len(line)
    if not chrom:
        raise ValueError(f"No FASTA record found in {path}")
    return chrom, size


def read_indel_events(map_path: str | Path, ref_chrom: str, query_chrom: str) -> list[GapEvent]:
    # Convert map INDEL rows into chain gap events on reference coordinates.
    events: list[GapEvent] = []
    with Path(map_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "ref_chr",
            "ref_start",
            "ref_end",
            "ref_allele",
            "sim_chr",
            "sim_allele",
            "variant_type",
            "variant_id",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Map file is missing columns: {', '.join(missing)}")
        for row in reader:
            if row["variant_type"] != "INDEL":
                continue
            if row["ref_chr"] != ref_chrom or row["sim_chr"] != query_chrom:
                continue
            event = _event_from_indel_row(row)
            if event.target_gap or event.query_gap:
                events.append(event)
    return sorted(events, key=lambda event: (event.ref_gap_start, event.variant_id))


def _event_from_indel_row(row: dict[str, str]) -> GapEvent:
    # Trim shared allele sequence so dt/dq represent only inserted/deleted bases.
    ref_allele = row["ref_allele"].upper()
    query_allele = row["sim_allele"].upper()
    prefix = _common_prefix_len(ref_allele, query_allele)
    suffix = _common_suffix_len(ref_allele[prefix:], query_allele[prefix:])
    target_gap = len(ref_allele) - prefix - suffix
    query_gap = len(query_allele) - prefix - suffix
    ref_gap_start = int(row["ref_start"]) - 1 + prefix
    return GapEvent(
        ref_gap_start=ref_gap_start,
        target_gap=target_gap,
        query_gap=query_gap,
        variant_id=row["variant_id"],
    )


def _common_prefix_len(left: str, right: str) -> int:
    # Count matching bases from the left side of two alleles.
    prefix = 0
    for left_base, right_base in zip(left, right):
        if left_base != right_base:
            break
        prefix += 1
    return prefix


def _common_suffix_len(left: str, right: str) -> int:
    # Count matching trailing bases after a shared prefix has been removed.
    suffix = 0
    max_suffix = min(len(left), len(right))
    while suffix < max_suffix and left[-suffix - 1] == right[-suffix - 1]:
        suffix += 1
    return suffix


def query_boundary_for_target(target_boundary: int, events: list[GapEvent]) -> int:
    # Translate a reference boundary to query space using completed gap events.
    offset = 0
    for event in events:
        if event.ref_gap_start >= target_boundary:
            break
        offset += event.query_gap - event.target_gap
    return target_boundary + offset


def build_chain(
    spec: ChainSpec,
    target_chrom: str,
    target_size: int,
    query_chrom: str,
    query_size: int,
    events: list[GapEvent],
) -> ChainBuild:
    # Build one plus-strand chain over the requested target interval.
    usable_events_by_start: dict[int, list[GapEvent]] = defaultdict(list)
    skipped_events = 0
    for event in events:
        if event.ref_gap_start < spec.target_start or event.ref_gap_start >= spec.target_end:
            continue
        if event.ref_gap_end >= spec.target_end:
            skipped_events += 1
            continue
        usable_events_by_start[event.ref_gap_start].append(event)

    merged_events = []
    for ref_gap_start in sorted(usable_events_by_start):
        events_at_start = usable_events_by_start[ref_gap_start]
        merged_events.append(
            GapEvent(
                ref_gap_start=ref_gap_start,
                target_gap=sum(event.target_gap for event in events_at_start),
                query_gap=sum(event.query_gap for event in events_at_start),
                variant_id=",".join(event.variant_id for event in events_at_start),
            )
        )

    q_start = query_boundary_for_target(spec.target_start, events)
    rows: list[str] = []
    current_target = spec.target_start
    query_span = 0
    aligned_bases = 0

    for event in merged_events:
        block_size = event.ref_gap_start - current_target
        if block_size < 0:
            skipped_events += 1
            continue
        if block_size == 0 and rows:
            size, target_gap, query_gap = [int(value) for value in rows[-1].split()]
            rows[-1] = f"{size} {target_gap + event.target_gap} {query_gap + event.query_gap}"
        elif block_size == 0:
            skipped_events += 1
            continue
        else:
            rows.append(f"{block_size} {event.target_gap} {event.query_gap}")
            aligned_bases += block_size
        current_target = event.ref_gap_start + event.target_gap
        query_span += block_size + event.query_gap

    final_size = spec.target_end - current_target
    if final_size <= 0:
        raise ValueError(f"Chain {spec.chain_id} has no final aligned block")
    rows.append(str(final_size))
    aligned_bases += final_size
    query_span += final_size
    q_end = q_start + query_span
    score = aligned_bases
    header = (
        f"chain {score} {target_chrom} {target_size} + "
        f"{spec.target_start} {spec.target_end} "
        f"{query_chrom} {query_size} + {q_start} {q_end} {spec.chain_id}"
    )
    return ChainBuild(
        text="\n".join([header, *rows, ""]),
        used_events=len(merged_events),
        skipped_events=skipped_events,
    )


def write_fake_chains(
    ref_fasta: str | Path,
    query_fasta: str | Path,
    map_path: str | Path,
    out_path: str | Path,
) -> None:
    # Write the two-chain fake alignment requested for the dummy genomes.
    ref_chrom, ref_size = read_fasta_size(ref_fasta)
    query_chrom, query_size = read_fasta_size(query_fasta)
    events = read_indel_events(map_path, ref_chrom, query_chrom)
    specs = [
        ChainSpec(chain_id=1, target_start=0, target_end=2000),
        ChainSpec(chain_id=2, target_start=2500, target_end=ref_size),
    ]

    builds = [
        build_chain(spec, ref_chrom, ref_size, query_chrom, query_size, events)
        for spec in specs
    ]
    Path(out_path).write_text("\n".join(build.text for build in builds), encoding="utf-8")

    represented = sum(build.used_events for build in builds)
    skipped_inside_chains = sum(build.skipped_events for build in builds)
    outside_chains = len(events) - represented - skipped_inside_chains
    print(f"reference: {ref_chrom} ({ref_size} bp)")
    print(f"query: {query_chrom} ({query_size} bp)")
    print(f"indel events in map: {len(events)}")
    print(f"indel events represented as chain gaps: {represented}")
    print(f"indel events skipped inside chain intervals: {skipped_inside_chains}")
    print(f"indel events in unaligned interval: {outside_chains}")
    print(f"wrote: {out_path}")


def build_parser() -> argparse.ArgumentParser:
    # Define command-line inputs with defaults for the dummy test data.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-fasta", type=Path, default=Path("tests/dummyRef.fasta"))
    parser.add_argument(
        "--query-fasta", type=Path, default=Path("tests/dummyQuery.simseq.fasta")
    )
    parser.add_argument(
        "--map",
        dest="map_path",
        type=Path,
        default=Path("tests/dummyQuery.refseq2simseq.map.txt"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("tests/dummyRef_to_dummyQuery.fake.chain"),
    )
    return parser


def main() -> None:
    # Generate the fake chain file and print a compact summary.
    args = build_parser().parse_args()
    write_fake_chains(args.ref_fasta, args.query_fasta, args.map_path, args.out)


if __name__ == "__main__":
    main()
