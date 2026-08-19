"""Single-command TOGA4CNEE workflow runner."""

from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .candidates import write_candidates
from .index import build_index
from .orthology import (
    DEFAULT_MIN_ALIGN_SCORE,
    DEFAULT_MIN_GAP_LENGTH,
    DEFAULT_MIN_NEIGHBOR_FRACTION,
    DEFAULT_VERIFICATION_NUM_NEIGHBORS,
    OrthologyResult,
    classify_orthology,
)
from .overlap import OVERLAP_COLUMNS, write_overlap_stats
from .sequence_compare import SequenceComparator
from .synteny import (
    DEFAULT_FLANK_SIZES,
    SyntenyResult,
    write_synteny_neighbor_verification,
    write_synteny_verification,
)
from .summary import stats_row_is_solid, summarize_overlap_stats


@dataclass(frozen=True)
class WorkflowResult:
    chain_count: int
    candidate_count: int
    stats_count: int
    strict_stats_count: int
    query_bed_count: int
    summary_count: int
    strict_summary_count: int
    orthology_result: OrthologyResult
    synteny_result: SyntenyResult
    synteny_neighbor_result: SyntenyResult
    index_path: Path
    candidates_path: Path
    synteny_path: Path
    synteny_neighbors_path: Path
    strict_stats_path: Path
    strict_summary_path: Path
    culled_projections_path: Path
    ref_cnees_path: Path


def run_workflow(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: str | Path,
    out_stats: str | Path,
    out_summary: str | Path | None = None,
    out_orthology: str | Path | None = None,
    out_strict_stats: str | Path | None = None,
    out_strict_summary: str | Path | None = None,
    out_synteny: str | Path | None = None,
    out_synteny_neighbors: str | Path | None = None,
    out_culled_projections: str | Path | None = None,
    out_ref_cnees: str | Path | None = None,
    ref_fasta: str | Path | None = None,
    query_fasta: str | Path | None = None,
    synteny_flank_sizes: Sequence[int] = DEFAULT_FLANK_SIZES,
    synteny_max_one2one_controls: int = 50000,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    synteny_check: bool = True,
    min_align_score: float = DEFAULT_MIN_ALIGN_SCORE,
    min_neighbor_fraction: float = DEFAULT_MIN_NEIGHBOR_FRACTION,
    verification_num_neighbors: int = DEFAULT_VERIFICATION_NUM_NEIGHBORS,
    min_gap_length: int = DEFAULT_MIN_GAP_LENGTH,
    work_dir: str | Path | None = None,
) -> WorkflowResult:
    """Run indexing, candidate discovery, overlap stats, and query BED export."""
    # Normalize output paths, validate paired FASTA inputs, and choose a work dir.
    if not synteny_flank_sizes or any(size < 0 for size in synteny_flank_sizes):
        raise ValueError("synteny_flank_sizes must be one or more non-negative values")
    if synteny_max_one2one_controls < 0:
        raise ValueError("synteny_max_one2one_controls must be non-negative")
    if (ref_fasta is None) != (query_fasta is None):
        raise ValueError("--ref-fasta and --query-fasta must be supplied together")

    out_bed = Path(out_bed)
    out_stats = Path(out_stats)
    out_summary_path = (
        Path(out_summary) if out_summary is not None else out_stats.with_suffix(".summary.tsv")
    )
    out_orthology_path = (
        Path(out_orthology)
        if out_orthology is not None
        else out_stats.with_suffix(".orthology.tsv")
    )
    out_strict_stats_path = (
        Path(out_strict_stats)
        if out_strict_stats is not None
        else _add_suffix_before_extension(out_stats, ".strict")
    )
    out_strict_summary_path = (
        Path(out_strict_summary)
        if out_strict_summary is not None
        else _add_suffix_before_extension(out_stats, ".strict.summary")
    )
    out_synteny_path = (
        Path(out_synteny)
        if out_synteny is not None
        else _add_suffix_before_extension(out_stats, ".synteny")
    )
    out_synteny_neighbors_path = (
        Path(out_synteny_neighbors)
        if out_synteny_neighbors is not None
        else _add_suffix_before_extension(out_stats, ".synteny_neighbors")
    )
    out_culled_projections_path = (
        Path(out_culled_projections)
        if out_culled_projections is not None
        else _add_suffix_before_extension(out_stats, ".culled_projections")
    )
    out_ref_cnees_path = (
        Path(out_ref_cnees)
        if out_ref_cnees is not None
        else _add_suffix_before_extension(out_stats, ".ref_cnees")
    )
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_stats.parent.mkdir(parents=True, exist_ok=True)
    out_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_orthology_path.parent.mkdir(parents=True, exist_ok=True)
    out_strict_stats_path.parent.mkdir(parents=True, exist_ok=True)
    out_strict_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_synteny_path.parent.mkdir(parents=True, exist_ok=True)
    out_synteny_neighbors_path.parent.mkdir(parents=True, exist_ok=True)
    out_culled_projections_path.parent.mkdir(parents=True, exist_ok=True)
    out_ref_cnees_path.parent.mkdir(parents=True, exist_ok=True)

    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="toga4cnee_") as tmpdir:
            return _run_workflow_in_dir(
                chain_file,
                cnee_bed,
                out_bed,
                out_stats,
                out_summary_path,
                out_orthology_path,
                out_strict_stats_path,
                out_strict_summary_path,
                out_synteny_path,
                out_synteny_neighbors_path,
                out_culled_projections_path,
                out_ref_cnees_path,
                ref_fasta,
                query_fasta,
                synteny_flank_sizes,
                synteny_max_one2one_controls,
                min_sequence_identity,
                max_fraction_gap,
                synteny_check,
                min_align_score,
                min_neighbor_fraction,
                verification_num_neighbors,
                min_gap_length,
                Path(tmpdir),
            )

    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    return _run_workflow_in_dir(
        chain_file,
        cnee_bed,
        out_bed,
        out_stats,
        out_summary_path,
        out_orthology_path,
        out_strict_stats_path,
        out_strict_summary_path,
        out_synteny_path,
        out_synteny_neighbors_path,
        out_culled_projections_path,
        out_ref_cnees_path,
        ref_fasta,
        query_fasta,
        synteny_flank_sizes,
        synteny_max_one2one_controls,
        min_sequence_identity,
        max_fraction_gap,
        synteny_check,
        min_align_score,
        min_neighbor_fraction,
        verification_num_neighbors,
        min_gap_length,
        work_path,
    )


def _run_workflow_in_dir(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: Path,
    out_stats: Path,
    out_summary: Path,
    out_orthology: Path,
    out_strict_stats: Path,
    out_strict_summary: Path,
    out_synteny: Path,
    out_synteny_neighbors: Path,
    out_culled_projections: Path,
    out_ref_cnees: Path,
    ref_fasta: str | Path | None,
    query_fasta: str | Path | None,
    synteny_flank_sizes: Sequence[int],
    synteny_max_one2one_controls: int,
    min_sequence_identity: float,
    max_fraction_gap: float,
    synteny_check: bool,
    min_align_score: float,
    min_neighbor_fraction: float,
    verification_num_neighbors: int,
    min_gap_length: int,
    work_dir: Path,
) -> WorkflowResult:
    # Execute the full workflow in order, keeping intermediates under work_dir.
    index_path = work_dir / "chain_index.json"
    candidates_path = work_dir / "candidate_chains.tsv"

    index_payload = build_index(chain_file, index_path)
    candidate_count = write_candidates(index_path, cnee_bed, candidates_path)
    sequence_comparator = (
        SequenceComparator(ref_fasta, query_fasta, work_dir)
        if ref_fasta is not None and query_fasta is not None
        else None
    )
    try:
        stats_count = write_overlap_stats(
            index_path,
            candidates_path,
            out_stats,
            sequence_comparator=sequence_comparator,
        )
        orthology_result = classify_orthology(
            out_stats,
            cnee_bed,
            index_path,
            out_orthology,
            out_bed,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
            synteny_check=synteny_check,
            sequence_comparator=sequence_comparator,
            min_align_score=min_align_score,
            min_neighbor_fraction=min_neighbor_fraction,
            verification_num_neighbors=verification_num_neighbors,
            out_culled_projections=out_culled_projections,
            query_fasta=query_fasta,
            min_gap_length=min_gap_length,
            out_ref_cnees=out_ref_cnees,
        )
        synteny_result = write_synteny_verification(
            index_path,
            out_stats,
            out_bed,
            out_orthology,
            cnee_bed,
            out_synteny,
            flank_sizes=synteny_flank_sizes,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
            max_one2one_controls=synteny_max_one2one_controls,
            sequence_comparator=sequence_comparator,
        )
        synteny_neighbor_result = write_synteny_neighbor_verification(
            index_path,
            out_stats,
            out_bed,
            out_orthology,
            cnee_bed,
            out_synteny_neighbors,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
            max_one2one_controls=synteny_max_one2one_controls,
        )
        summary_count = summarize_overlap_stats(
            out_stats,
            out_summary,
            cnee_bed=cnee_bed,
            require_solid_reference=True,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
        strict_stats_count = _write_strict_stats(
            out_stats,
            out_strict_stats,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
        strict_summary_count = summarize_overlap_stats(
            out_stats,
            out_strict_summary,
            cnee_bed=cnee_bed,
            solid_only=True,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
    finally:
        if sequence_comparator is not None:
            sequence_comparator.close()

    return WorkflowResult(
        chain_count=int(index_payload["chain_count"]),
        candidate_count=candidate_count,
        stats_count=stats_count,
        strict_stats_count=strict_stats_count,
        query_bed_count=orthology_result.query_bed_rows,
        summary_count=summary_count,
        strict_summary_count=strict_summary_count,
        orthology_result=orthology_result,
        synteny_result=synteny_result,
        synteny_neighbor_result=synteny_neighbor_result,
        index_path=index_path,
        candidates_path=candidates_path,
        synteny_path=out_synteny,
        synteny_neighbors_path=out_synteny_neighbors,
        strict_stats_path=out_strict_stats,
        strict_summary_path=out_strict_summary,
        culled_projections_path=out_culled_projections,
        ref_cnees_path=out_ref_cnees,
    )


def _write_strict_stats(
    stats_path: Path,
    out_stats: Path,
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> int:
    # Write a filtered overlap-stats TSV containing only solid candidate rows.
    count = 0
    with stats_path.open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        with out_stats.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=reader.fieldnames or OVERLAP_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                if not stats_row_is_solid(
                    row,
                    min_sequence_identity=min_sequence_identity,
                    max_fraction_gap=max_fraction_gap,
                ):
                    continue
                writer.writerow(row)
                count += 1
    return count


def _add_suffix_before_extension(path: Path, suffix: str) -> Path:
    # Insert a logical suffix before the final file extension.
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")
