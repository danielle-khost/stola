"""Single-command STOLA workflow runner."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

from .bed import bed_missing_name_column, parse_bed, write_named_bed
from .candidates import write_candidates
from .index import build_index
from .orthology import (
    DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    DEFAULT_FLANK_WINDOW,
    DEFAULT_MIN_GAP_LENGTH,
    DEFAULT_QUERY_ID_PREFIX,
    DEFAULT_RELATIVE_SCORE_THRESHOLD,
    DEFAULT_SOFT_DEPTH_LIMIT,
    OrthologyResult,
    classify_orthology,
)
from .overlap import write_overlap_stats
from .plots import OrthologySummaryResult, PlotResult, write_orthology_summary, write_qc_plots
from .sequence_compare import SequenceComparator
from .summary import summarize_overlap_stats

# Fixed output layout under --outdir -- see the README's Output section. Not individually
# overridable; --outdir is the only output-location knob the CLI exposes.
QUERY_BED_FILENAME = "query_cnees.bed"
STATS_FILENAME = "overlap_stats.tsv"
# Summary tables and QC plots, under SUMMARY_DIRNAME (plus the orthology
# summary .svg/.tsv pairs named in plots.py's *_NAME constants).
SUMMARY_DIRNAME = "summary_figures"
SUMMARY_FILENAME = "overlap_stats.summary.tsv"
ALIGNMENT_HISTOGRAM_FILENAME = "num_query_alignments.svg"
OVERLAP_DENSITY_FILENAME = "overlap_fraction_density.svg"
SEQUENCE_IDENTITY_DENSITY_FILENAME = "sequence_identity_density.svg"
QUERY_ORTHOLOGY_FILENAME = "query_ortho_status.tsv"
PRUNED_PROJECTIONS_FILENAME = "pruned_query_projections.tsv"
REF_ORTHOLOGY_FILENAME = "ref_ortho_status.tsv"
WORK_DIRNAME = "work"
NAMED_BED_SUFFIX = "_named.bed"


@dataclass(frozen=True)
class WorkflowResult:
    chain_count: int
    candidate_count: int
    stats_count: int
    query_bed_count: int
    summary_count: int
    summary_plots: PlotResult
    orthology_summary: OrthologySummaryResult
    orthology_result: OrthologyResult
    outdir: Path
    index_path: Path
    candidates_path: Path
    query_bed_path: Path
    stats_path: Path
    summary_dir: Path
    summary_path: Path
    orthology_path: Path
    culled_projections_path: Path
    ref_cnees_path: Path
    named_ref_bed_path: Path | None


def run_workflow(
    chain_file: str | Path,
    cnee_bed: str | Path,
    ref_fasta: str | Path,
    query_fasta: str | Path,
    outdir: str | Path = ".",
    query_id_prefix: str = DEFAULT_QUERY_ID_PREFIX,
    min_sequence_identity: float = 0.7,
    max_fraction_gap: float = 0.1,
    filter_multi_projections: bool = True,
    soft_depth_limit: int = DEFAULT_SOFT_DEPTH_LIMIT,
    relative_score_threshold: float = DEFAULT_RELATIVE_SCORE_THRESHOLD,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    flank_any_coverage_threshold: float = DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    min_gap_length: int = DEFAULT_MIN_GAP_LENGTH,
) -> WorkflowResult:
    """Run the full chain/CNEE workflow, writing every output under outdir."""
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    work_dir = outdir_path / WORK_DIRNAME
    work_dir.mkdir(parents=True, exist_ok=True)

    out_bed = outdir_path / QUERY_BED_FILENAME
    out_stats = outdir_path / STATS_FILENAME
    summary_dir = outdir_path / SUMMARY_DIRNAME
    summary_dir.mkdir(parents=True, exist_ok=True)
    out_summary = summary_dir / SUMMARY_FILENAME
    out_orthology = outdir_path / QUERY_ORTHOLOGY_FILENAME
    out_culled_projections = outdir_path / PRUNED_PROJECTIONS_FILENAME
    out_ref_cnees = outdir_path / REF_ORTHOLOGY_FILENAME

    index_path = work_dir / "chain_index.json"
    candidates_path = work_dir / "candidate_chains.tsv"

    named_ref_bed_path: Path | None = None
    if bed_missing_name_column(cnee_bed):
        # The reference BED lacked a name column, so parse_bed synthesized
        # coordinate-based ids for it. Those ids are already stable across
        # runs on this same file, but write them out so users can see them
        # and compare them directly across separate query runs.
        named_ref_bed_path = outdir_path / f"{Path(cnee_bed).stem}{NAMED_BED_SUFFIX}"
        with warnings.catch_warnings():
            # parse_bed's duplicate/zero-length warnings already fire from
            # the downstream steps below; avoid emitting them twice here.
            warnings.simplefilter("ignore")
            named_intervals = parse_bed(cnee_bed)
        write_named_bed(named_intervals, named_ref_bed_path)

    index_payload = build_index(chain_file, index_path)
    candidate_count = write_candidates(index_path, cnee_bed, candidates_path)
    sequence_comparator = SequenceComparator(ref_fasta, query_fasta, work_dir)
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
            sequence_comparator=sequence_comparator,
            filter_multi_projections=filter_multi_projections,
            soft_depth_limit=soft_depth_limit,
            relative_score_threshold=relative_score_threshold,
            flank_window=flank_window,
            flank_any_coverage_threshold=flank_any_coverage_threshold,
            out_culled_projections=out_culled_projections,
            query_fasta=query_fasta,
            min_gap_length=min_gap_length,
            out_ref_cnees=out_ref_cnees,
            query_id_prefix=query_id_prefix,
        )
        summary_count = summarize_overlap_stats(
            out_stats,
            out_summary,
            cnee_bed=cnee_bed,
            require_solid_reference=True,
            min_sequence_identity=min_sequence_identity,
            max_fraction_gap=max_fraction_gap,
        )
    finally:
        sequence_comparator.close()
    summary_plots = write_qc_plots(
        out_summary,
        out_stats,
        summary_dir / ALIGNMENT_HISTOGRAM_FILENAME,
        summary_dir / OVERLAP_DENSITY_FILENAME,
        summary_dir / SEQUENCE_IDENTITY_DENSITY_FILENAME,
    )
    orthology_summary = write_orthology_summary(out_ref_cnees, out_culled_projections, summary_dir)

    return WorkflowResult(
        chain_count=int(index_payload["chain_count"]),
        candidate_count=candidate_count,
        stats_count=stats_count,
        query_bed_count=orthology_result.query_bed_rows,
        summary_count=summary_count,
        summary_plots=summary_plots,
        orthology_summary=orthology_summary,
        orthology_result=orthology_result,
        outdir=outdir_path,
        index_path=index_path,
        candidates_path=candidates_path,
        query_bed_path=out_bed,
        stats_path=out_stats,
        summary_dir=summary_dir,
        summary_path=out_summary,
        orthology_path=out_orthology,
        culled_projections_path=out_culled_projections,
        ref_cnees_path=out_ref_cnees,
        named_ref_bed_path=named_ref_bed_path,
    )
