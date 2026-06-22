"""Single-command TOGA4CNEE workflow runner."""

from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .candidates import write_candidates
from .index import build_index
from .overlap import OVERLAP_COLUMNS, iter_overlap_rows
from .plots import PlotResult, write_qc_plots
from .query_bed import row_to_bed12
from .sequence_compare import SequenceComparator
from .summary import summarize_overlap_stats


@dataclass(frozen=True)
class WorkflowResult:
    chain_count: int
    candidate_count: int
    stats_count: int
    query_bed_count: int
    summary_count: int
    plot_result: PlotResult
    index_path: Path
    candidates_path: Path


def run_workflow(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: str | Path,
    out_stats: str | Path,
    out_summary: str | Path | None = None,
    out_alignment_histogram: str | Path | None = None,
    out_overlap_density: str | Path | None = None,
    out_sequence_identity_density: str | Path | None = None,
    ref_fasta: str | Path | None = None,
    query_fasta: str | Path | None = None,
    flank_size: int = 1000,
    work_dir: str | Path | None = None,
) -> WorkflowResult:
    """Run indexing, candidate discovery, overlap stats, and query BED export."""
    if flank_size < 0:
        raise ValueError("flank_size must be non-negative")
    if (ref_fasta is None) != (query_fasta is None):
        raise ValueError("--ref-fasta and --query-fasta must be supplied together")

    out_bed = Path(out_bed)
    out_stats = Path(out_stats)
    out_summary_path = (
        Path(out_summary) if out_summary is not None else out_stats.with_suffix(".summary.tsv")
    )
    out_alignment_histogram_path = (
        Path(out_alignment_histogram)
        if out_alignment_histogram is not None
        else out_stats.with_suffix(".num_query_alignments.svg")
    )
    out_overlap_density_path = (
        Path(out_overlap_density)
        if out_overlap_density is not None
        else out_stats.with_suffix(".overlap_fraction_density.svg")
    )
    out_sequence_identity_density_path = (
        Path(out_sequence_identity_density)
        if out_sequence_identity_density is not None
        else out_stats.with_suffix(".sequence_identity_density.svg")
    )
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_stats.parent.mkdir(parents=True, exist_ok=True)
    out_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_alignment_histogram_path.parent.mkdir(parents=True, exist_ok=True)
    out_overlap_density_path.parent.mkdir(parents=True, exist_ok=True)
    out_sequence_identity_density_path.parent.mkdir(parents=True, exist_ok=True)

    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="toga4cnee_") as tmpdir:
            return _run_workflow_in_dir(
                chain_file,
                cnee_bed,
                out_bed,
                out_stats,
                out_summary_path,
                out_alignment_histogram_path,
                out_overlap_density_path,
                out_sequence_identity_density_path,
                ref_fasta,
                query_fasta,
                flank_size,
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
        out_alignment_histogram_path,
        out_overlap_density_path,
        out_sequence_identity_density_path,
        ref_fasta,
        query_fasta,
        flank_size,
        work_path,
    )


def _run_workflow_in_dir(
    chain_file: str | Path,
    cnee_bed: str | Path,
    out_bed: Path,
    out_stats: Path,
    out_summary: Path,
    out_alignment_histogram: Path,
    out_overlap_density: Path,
    out_sequence_identity_density: Path,
    ref_fasta: str | Path | None,
    query_fasta: str | Path | None,
    flank_size: int,
    work_dir: Path,
) -> WorkflowResult:
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
        stats_count, query_bed_count = _write_outputs(
            index_path,
            candidates_path,
            out_stats,
            out_bed,
            flank_size,
            sequence_comparator=sequence_comparator,
        )
    finally:
        if sequence_comparator is not None:
            sequence_comparator.close()
    summary_count = summarize_overlap_stats(out_stats, out_summary, cnee_bed=cnee_bed)
    plot_result = write_qc_plots(
        out_summary,
        out_stats,
        out_alignment_histogram,
        out_overlap_density,
        out_sequence_identity_density,
    )

    return WorkflowResult(
        chain_count=int(index_payload["chain_count"]),
        candidate_count=candidate_count,
        stats_count=stats_count,
        query_bed_count=query_bed_count,
        summary_count=summary_count,
        plot_result=plot_result,
        index_path=index_path,
        candidates_path=candidates_path,
    )


def _write_outputs(
    index_path: Path,
    candidates_path: Path,
    out_stats: Path,
    out_bed: Path,
    flank_size: int,
    sequence_comparator: SequenceComparator | None = None,
) -> tuple[int, int]:
    stats_count = 0
    query_bed_count = 0
    rows = iter_overlap_rows(
        index_path,
        candidates_path,
        flank_size=flank_size,
        sequence_comparator=sequence_comparator,
    )
    with out_stats.open("w", encoding="utf-8", newline="") as stats_handle:
        writer = csv.DictWriter(stats_handle, fieldnames=OVERLAP_COLUMNS, delimiter="\t")
        writer.writeheader()
        with out_bed.open("w", encoding="utf-8") as bed_handle:
            for row in rows:
                writer.writerow(row)
                stats_count += 1
                bed_line = row_to_bed12(row)
                if bed_line is None:
                    continue
                bed_handle.write(bed_line + "\n")
                query_bed_count += 1
    return stats_count, query_bed_count
