"""Concise per-reference-CNEE summaries from overlap statistics."""

from __future__ import annotations

import csv
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

from .bed import parse_bed

SUMMARY_COLUMNS = [
    "cnee_uid",
    "projection_status",
    "num_query_alignments",
    "overlap_fractions",
    "gap_bases",
    "sequence_identities",
    "sequence_mismatches",
]


@dataclass
class CneeSummary:
    cnee_uid: str
    overlap_fractions: list[str] = field(default_factory=list)
    gap_bases: list[str] = field(default_factory=list)
    sequence_identities: list[str] = field(default_factory=list)
    sequence_mismatches: list[str] = field(default_factory=list)
    candidate_count: int = 0
    aligned_candidate_count: int = 0
    solid_candidate_count: int = 0
    solid_only: bool = False
    require_solid_reference: bool = False

    @property
    def num_query_alignments(self) -> int:
        # Count aligned query projections represented in this summary row.
        return len(self.overlap_fractions)

    @property
    def projection_status(self) -> str:
        # Distinguish aligned, gap-only/lost, and no-candidate/missing CNEEs.
        # With require_solid_reference, reported stats include every aligned
        # row even for weak-only refs, so "aligned" must key off solid
        # candidates specifically rather than num_query_alignments.
        if self.require_solid_reference:
            aligned = self.solid_candidate_count > 0
        else:
            aligned = self.num_query_alignments > 0
        if aligned:
            return "aligned"
        if (
            (self.solid_only or self.require_solid_reference)
            and self.aligned_candidate_count > 0
        ):
            return "uncertain_loss"
        if self.candidate_count > 0:
            return "lost"
        return "missing"

    def row(self) -> dict[str, str | int]:
        # Convert collected per-CNEE values into the output TSV row shape.
        return {
            "cnee_uid": self.cnee_uid,
            "projection_status": self.projection_status,
            "num_query_alignments": self.num_query_alignments,
            "overlap_fractions": ",".join(self.overlap_fractions),
            "gap_bases": ",".join(self.gap_bases),
            "sequence_identities": ",".join(self.sequence_identities),
            "sequence_mismatches": ",".join(self.sequence_mismatches),
        }


def summarize_overlap_stats(
    stats_path: str | Path,
    out_path: str | Path,
    cnee_bed: str | Path | None = None,
    aligned_only: bool = True,
    solid_only: bool = False,
    require_solid_reference: bool = False,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
) -> int:
    """Write a compact per-CNEE summary from an overlap stats TSV.

    By default, only rows with at least one aligned CNEE base are counted as
    query alignments. With solid_only, only rows passing the solid-match
    thresholds are counted as query alignments. With require_solid_reference,
    refs with at least one solid row are marked aligned and report every
    aligned row's stats; weak-only refs are marked uncertain_loss but still
    report their aligned rows' stats (overlap fraction, gap bases, etc.) as
    context rather than being left blank.
    """
    # Aggregate overlap, gap, and sequence metrics by reference CNEE ID.
    summaries: OrderedDict[str, CneeSummary] = OrderedDict()
    if cnee_bed is not None:
        for interval in parse_bed(cnee_bed):
            summaries.setdefault(
                interval.cnee_uid,
                CneeSummary(
                    cnee_uid=interval.cnee_uid,
                    solid_only=solid_only,
                    require_solid_reference=require_solid_reference,
                ),
            )

    with Path(stats_path).open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        required = {
            "cnee_uid",
            "aligned_bases",
            "covered_fraction",
            "gap_bases_within_chain_span",
        }
        if solid_only or require_solid_reference:
            required.add("gap_fraction_within_chain_span")
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Stats file is missing columns: {', '.join(missing)}")
        for row in reader:
            cnee_uid = row["cnee_uid"]
            summary = summaries.setdefault(
                cnee_uid,
                CneeSummary(
                    cnee_uid=cnee_uid,
                    solid_only=solid_only,
                    require_solid_reference=require_solid_reference,
                ),
            )
            summary.solid_only = solid_only
            summary.require_solid_reference = require_solid_reference
            summary.candidate_count += 1
            aligned_bases = int(row["aligned_bases"])
            if aligned_bases > 0:
                summary.aligned_candidate_count += 1
            if (solid_only or require_solid_reference) and stats_row_is_solid(
                row,
                min_sequence_identity=min_sequence_identity,
                max_fraction_gap=max_fraction_gap,
            ):
                summary.solid_candidate_count += 1

    with Path(stats_path).open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        for row in reader:
            summary = summaries[row["cnee_uid"]]
            aligned_bases = int(row["aligned_bases"])
            if aligned_only and aligned_bases <= 0:
                continue
            row_is_solid = (
                stats_row_is_solid(
                    row,
                    min_sequence_identity=min_sequence_identity,
                    max_fraction_gap=max_fraction_gap,
                )
                if solid_only
                else False
            )
            if solid_only and not row_is_solid:
                continue
            # Weak-only refs under require_solid_reference still report their
            # aligned rows' stats as uncertain_loss context (see
            # projection_status above); only solid_only stays strict.
            summary.overlap_fractions.append(row["covered_fraction"])
            summary.gap_bases.append(row["gap_bases_within_chain_span"])
            seq_identity = row.get("seq_identity", ".")
            if seq_identity != ".":
                summary.sequence_identities.append(seq_identity)
            if "seq_mismatches" in row and seq_identity != ".":
                summary.sequence_mismatches.append(row["seq_mismatches"])

    with Path(out_path).open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for summary in summaries.values():
            writer.writerow(summary.row())
    return len(summaries)


def is_solid_match(
    sequence_identity: float | None,
    fraction_gap: float,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    """Apply the shared solid-match criteria to one (identity, gap) pair."""
    # The single source of truth for "solid" used by summary/orthology/plots.
    if fraction_gap >= max_fraction_gap:
        return False
    if sequence_identity is None:
        return True
    return sequence_identity > min_sequence_identity


def stats_row_is_solid(
    row: dict[str, str],
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
) -> bool:
    # Apply solid-match criteria directly to one overlap-stats row.
    if int(row["aligned_bases"]) <= 0:
        return False
    seq_identity_text = row.get("seq_identity", ".")
    seq_identity = None if seq_identity_text in {"", "."} else float(seq_identity_text)
    return is_solid_match(
        seq_identity,
        float(row["gap_fraction_within_chain_span"]),
        min_sequence_identity,
        max_fraction_gap,
    )
