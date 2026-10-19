"""Single-command TOGA4CNEE workflow runner."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

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
from .sequence_compare import SequenceComparator
from .synteny import (
    DEFAULT_FLANK_SIZES,
    SyntenyResult,
    write_synteny_neighbor_verification,
    write_synteny_verification,
)
from .summary import summarize_overlap_stats

# Fixed output layout under --outdir -- see output_files.md. Not individually
# overridable; --outdir is the only output-location knob the CLI exposes. The
# synteny reports are exploratory/diagnostic, not main output, so they live
# under work/ alongside the other intermediates rather than in outdir itself.
QUERY_BED_FILENAME = "query_cnees.bed"
STATS_FILENAME = "overlap_stats.tsv"
SUMMARY_FILENAME = "overlap_stats.summary.tsv"
SYNTENY_FILENAME = "overlap_stats.synteny.tsv"
SYNTENY_NEIGHBORS_FILENAME = "overlap_stats.synteny_neighbors.tsv"
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
    orthology_result: OrthologyResult
    synteny_result: SyntenyResult
    synteny_neighbor_result: SyntenyResult
    outdir: Path
    index_path: Path
    candidates_path: Path
    query_bed_path: Path
    stats_path: Path
    summary_path: Path
    synteny_path: Path
    synteny_neighbors_path: Path
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
    synteny_flank_sizes: Sequence[int] = DEFAULT_FLANK_SIZES,
    synteny_max_one2one_controls: int = 50000,
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
    if not synteny_flank_sizes or any(size < 0 for size in synteny_flank_sizes):
        raise ValueError("synteny_flank_sizes must be one or more non-negative values")
    if synteny_max_one2one_controls < 0:
        raise ValueError("synteny_max_one2one_controls must be non-negative")

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    work_dir = outdir_path / WORK_DIRNAME
    work_dir.mkdir(parents=True, exist_ok=True)

    out_bed = outdir_path / QUERY_BED_FILENAME
    out_stats = outdir_path / STATS_FILENAME
    out_summary = outdir_path / SUMMARY_FILENAME
    out_synteny = work_dir / SYNTENY_FILENAME
    out_synteny_neighbors = work_dir / SYNTENY_NEIGHBORS_FILENAME
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
    finally:
        sequence_comparator.close()

    return WorkflowResult(
        chain_count=int(index_payload["chain_count"]),
        candidate_count=candidate_count,
        stats_count=stats_count,
        query_bed_count=orthology_result.query_bed_rows,
        summary_count=summary_count,
        orthology_result=orthology_result,
        synteny_result=synteny_result,
        synteny_neighbor_result=synteny_neighbor_result,
        outdir=outdir_path,
        index_path=index_path,
        candidates_path=candidates_path,
        query_bed_path=out_bed,
        stats_path=out_stats,
        summary_path=out_summary,
        synteny_path=out_synteny,
        synteny_neighbors_path=out_synteny_neighbors,
        orthology_path=out_orthology,
        culled_projections_path=out_culled_projections,
        ref_cnees_path=out_ref_cnees,
        named_ref_bed_path=named_ref_bed_path,
    )
