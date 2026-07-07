"""Single-command TOGA4CNEE workflow runner."""

from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .candidates import write_candidates
from .index import build_index
from .orthology import OrthologyResult, classify_orthology
from .overlap import OVERLAP_COLUMNS, iter_overlap_rows
from .sequence_compare import SequenceComparator
from .summary import summarize_overlap_stats


@dataclass(frozen=True)
class WorkflowResult:
    chain_count: int
    candidate_count: int
    stats_count: int
    query_bed_count: int
    summary_count: int
    orthology_result: OrthologyResult
    index_path: Path
    candidates_path: Path


def run_workflow(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: str | Path,
    out_stats: str | Path,
    out_summary: str | Path | None = None,
    out_orthology: str | Path | None = None,
    ref_fasta: str | Path | None = None,
    query_fasta: str | Path | None = None,
    flank_size: int = 1000,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    work_dir: str | Path | None = None,
) -> WorkflowResult:
    """Run indexing, candidate discovery, overlap stats, and query BED export."""
    # Normalize output paths, validate paired FASTA inputs, and choose a work dir.
    if flank_size < 0:
        raise ValueError("flank_size must be non-negative")
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
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_stats.parent.mkdir(parents=True, exist_ok=True)
    out_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_orthology_path.parent.mkdir(parents=True, exist_ok=True)

    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="toga4cnee_") as tmpdir:
            return _run_workflow_in_dir(
                chain_file,
                cnee_bed,
                out_bed,
                out_stats,
                out_summary_path,
                out_orthology_path,
                ref_fasta,
                query_fasta,
                flank_size,
                min_sequence_identity,
                max_fraction_gap,
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
        ref_fasta,
        query_fasta,
        flank_size,
        min_sequence_identity,
        max_fraction_gap,
        work_path,
    )


def _run_workflow_in_dir(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: Path,
    out_stats: Path,
    out_summary: Path,
    out_orthology: Path,
    ref_fasta: str | Path | None,
    query_fasta: str | Path | None,
    flank_size: int,
    min_sequence_identity: float,
    max_fraction_gap: float,
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
        stats_count = _write_stats(
            index_path,
            candidates_path,
            out_stats,
            flank_size,
            sequence_comparator=sequence_comparator,
        )
    finally:
        if sequence_comparator is not None:
            sequence_comparator.close()
    orthology_result = classify_orthology(
        out_stats,
        cnee_bed,
        out_orthology,
        out_bed,
        min_sequence_identity=min_sequence_identity,
        max_fraction_gap=max_fraction_gap,
    )
    summary_count = summarize_overlap_stats(out_stats, out_summary, cnee_bed=cnee_bed)

    return WorkflowResult(
        chain_count=int(index_payload["chain_count"]),
        candidate_count=candidate_count,
        stats_count=stats_count,
        query_bed_count=orthology_result.query_bed_rows,
        summary_count=summary_count,
        orthology_result=orthology_result,
        index_path=index_path,
        candidates_path=candidates_path,
    )


def _write_stats(
    index_path: Path,
    candidates_path: Path,
    out_stats: Path,
    flank_size: int,
    sequence_comparator: SequenceComparator | None = None,
) -> int:
    # Stream overlap-stat rows to disk, optionally adding sequence metrics.
    stats_count = 0
    rows = iter_overlap_rows(
        index_path,
        candidates_path,
        flank_size=flank_size,
        sequence_comparator=sequence_comparator,
    )
    with out_stats.open("w", encoding="utf-8", newline="") as stats_handle:
        writer = csv.DictWriter(stats_handle, fieldnames=OVERLAP_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            stats_count += 1
    return stats_count
