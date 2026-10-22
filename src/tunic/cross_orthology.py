"""Compare per-reference-CNEE orthology status across multiple pipeline runs.

Takes a list of `toga4cnee run` output directories -- the same reference BED
lifted to different query species -- and offers three run modes (see
dev_docs/compare_cross_orthology.md):

    summary: build one combined TSV where each row is a reference CNEE and
        each column is a query species' final orthology_category for that
        CNEE, read from each directory's ref_ortho_status.tsv, plus two
        bar-chart SVGs of loss/multi-projection frequency across samples.
    sequence: extract FASTA sequence for reference CNEEs that are one2one
        in enough samples (per the `summary`-mode orthology table), one
        <ref_cnee_id>_orthologs.fasta per CNEE.
    cluster: for each ref CNEE lost in enough species (per the `summary`-mode
        orthology table), compute the D statistic of phylogenetic signal
        (Fritz & Purvis 2010) for its present/lost pattern on a given tree.
        Needs toytree/pandas/matplotlib (installed with the conda package;
        with pip, the `cluster` extra); not a dependency of the other two
        modes.

Installed as the `toga4cnee-compare-cross-orthology` command.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from .bed import parse_bed
from .genome_io import open_genome_accessor
from .orthology import ORTHOLOGY_COLUMNS, REF_CNEE_COLUMNS
from .plots import _bar_svg
from .workflow import (
    QUERY_BED_FILENAME,
    QUERY_ORTHOLOGY_FILENAME,
    REF_ORTHOLOGY_FILENAME,
    WORK_DIRNAME,
)

GENERIC_LABEL = "q"  # DEFAULT_QUERY_ID_PREFIX ("q_") stripped of its trailing "_".
FASTA_SUBDIR_NAME = "one2one_ortholog_fasta"
FASTA_BATCH_SIZE = 10000
FASTA_LINE_WIDTH = 60

# Reference-only categories from orthology.py -- no query ortholog was found
# in that sample at all. Any category not recognized as one2one or
# multi-projection (including a table cell that's "NA" or some future
# category) is conservatively treated as lost too, rather than silently
# counted as a usable ortholog.
LOST_CATEGORIES = frozenset(
    {"unaligned", "probable_loss", "assembly_gap_loss", "uncertain_loss", "missing/loss"}
)
MULTI_CATEGORIES = frozenset({"one2many", "many2one", "complex/many2many"})


def _category_kind(category: str | None) -> str:
    if category == "one2one":
        return "one2one"
    if category in MULTI_CATEGORIES:
        return "multi"
    return "lost"


@dataclass
class RunColumn:
    label: str
    category_by_ref: dict[str, str]
    run_dir: Path
    genome_path: Path | None
    # one2one_map (ref_cnee_id -> query_cnee_id) is intentionally NOT a field
    # here -- it's read fresh from QUERY_ORTHOLOGY_FILENAME on demand, one
    # column at a time, inside extract_one2one_orthologs's Pass 2, instead of
    # being pre-loaded for every column up front. That file can have order-
    # 10^5 one2one rows per sample; loading and keeping all of them for
    # hundreds of samples at once (holding all their one2one_maps in memory
    # simultaneously, most of them long before that column's turn in Pass 2)
    # is what actually exhausted memory on real multi-hundred-sample runs --
    # far more than the genome-accessor memory this module already guards.


def _read_run_dirs(list_path: str | Path) -> list[tuple[Path, Path | None, str | None]]:
    # Each line is "outdir", "outdir<TAB>genome", or "outdir<TAB>genome<TAB>label";
    # a blank genome/label column is treated the same as an absent one.
    run_dirs: list[tuple[Path, Path | None, str | None]] = []
    with Path(list_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            run_dir = Path(fields[0].strip())
            genome_path = Path(fields[1].strip()) if len(fields) >= 2 and fields[1].strip() else None
            manual_label = fields[2].strip() if len(fields) >= 3 and fields[2].strip() else None
            run_dirs.append((run_dir, genome_path, manual_label))
    return run_dirs


def _read_query_orthology_file(
    query_orthology_path: Path, build_one2one_map: bool
) -> tuple[str | None, dict[str, str]]:
    # "Peek" at the first real (non-"NA") query_cnee_id to recover the
    # --query-prefix used for this run. When build_one2one_map is set (FASTA
    # extraction is on), also collect every one2one ref->query id pair in
    # the same pass, rather than re-reading this file a second time.
    peeked_label: str | None = None
    seen_first_non_na = False
    one2one_map: dict[str, str] = {}
    with query_orthology_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(ORTHOLOGY_COLUMNS).difference(reader.fieldnames or []))
        if missing:
            raise ValueError(
                f"{query_orthology_path}: missing column(s) {', '.join(missing)}"
            )
        for row in reader:
            query_cnee_id = row["query_cnee_id"]
            has_query_id = bool(query_cnee_id) and query_cnee_id != "NA"
            if not seen_first_non_na and has_query_id:
                peeked_label = query_cnee_id.split("_", 1)[0]
                seen_first_non_na = True
                if not build_one2one_map:
                    break
            if build_one2one_map and has_query_id and row["orthology_category"] == "one2one":
                ref_id = row["partner_ref_cnees"]
                if ref_id and "," not in ref_id:
                    one2one_map[ref_id] = query_cnee_id
    return peeked_label, one2one_map


def _read_ref_orthology(ref_orthology_path: Path) -> dict[str, str]:
    category_by_ref: dict[str, str] = {}
    with ref_orthology_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(REF_CNEE_COLUMNS).difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{ref_orthology_path}: missing column(s) {', '.join(missing)}")
        for row in reader:
            category_by_ref[row["ref_cnee_id"]] = row["orthology_category"]
    return category_by_ref


def _read_query_bed_coords(
    query_bed_path: Path, wanted_ids: set[str]
) -> dict[str, tuple[str, int, int, str]]:
    # Single pass, restricted to the ids we actually need (the cross-run
    # one2one set), so this stays cheap even though query_cnees.bed can be
    # tens of millions of rows -- exits early once every id is found.
    coords: dict[str, tuple[str, int, int, str]] = {}
    if not wanted_ids:
        return coords
    remaining = set(wanted_ids)
    with Path(query_bed_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            name = fields[3]
            if name in remaining:
                coords[name] = (fields[0], int(fields[1]), int(fields[2]), fields[5])
                remaining.discard(name)
                if not remaining:
                    break
    return coords


def _dedupe_label(label: str, used: set[str], run_dir: Path) -> str:
    if label not in used:
        return label
    fallback = run_dir.resolve().name
    if fallback not in used:
        return fallback
    suffix = 2
    while f"{fallback}__{suffix}" in used:
        suffix += 1
    return f"{fallback}__{suffix}"


def _sanitize_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def collect_cross_run_orthology(
    run_dirs: list[tuple[Path, Path | None, str | None]],
    build_fasta_data: bool,
    require_ref_ortho_status: bool = True,
) -> list[RunColumn]:
    """require_ref_ortho_status controls whether each run's
    REF_ORTHOLOGY_FILENAME is required and loaded into the returned column's
    category_by_ref. `sequence` mode gets its categories from an already-built
    `summary`-mode orthology table instead (see _read_combined_orthology_table),
    so it passes False here and fills category_by_ref in from that table
    itself; category_by_ref comes back empty per column when this is False."""
    columns: list[RunColumn] = []
    used_labels: set[str] = set()
    for run_dir, genome_path, manual_label in run_dirs:
        ref_orthology_path = run_dir / REF_ORTHOLOGY_FILENAME
        query_orthology_path = run_dir / QUERY_ORTHOLOGY_FILENAME
        if require_ref_ortho_status and not ref_orthology_path.is_file():
            print(f"warning: skipping {run_dir} (missing {REF_ORTHOLOGY_FILENAME})", file=sys.stderr)
            continue

        if build_fasta_data and not query_orthology_path.is_file():
            # FASTA mode needs this file later (its one2one_map is read on
            # demand per column during extraction, not here), so a manual
            # label can no longer substitute for its absence.
            print(
                f"warning: skipping {run_dir} (missing {QUERY_ORTHOLOGY_FILENAME}, "
                "needed to locate one2one query orthologs for FASTA extraction)",
                file=sys.stderr,
            )
            continue

        if manual_label is not None:
            label = manual_label
        elif query_orthology_path.is_file():
            peeked_label, _ = _read_query_orthology_file(query_orthology_path, False)
            label = (
                peeked_label if peeked_label and peeked_label != GENERIC_LABEL else run_dir.resolve().name
            )
        else:
            print(
                f"warning: {run_dir} has no manual label and no {QUERY_ORTHOLOGY_FILENAME} "
                "to peek at; falling back to its directory name",
                file=sys.stderr,
            )
            label = run_dir.resolve().name

        label = _dedupe_label(label, used_labels, run_dir)
        used_labels.add(label)
        category_by_ref = _read_ref_orthology(ref_orthology_path) if require_ref_ortho_status else {}
        columns.append(
            RunColumn(
                label=label,
                category_by_ref=category_by_ref,
                run_dir=run_dir,
                genome_path=genome_path,
            )
        )

    return columns


def _read_combined_orthology_table(table_path: str | Path) -> dict[str, dict[str, str]]:
    """Read a `summary`-mode combined orthology TSV back in, returning
    {species_label: {ref_cnee_id: orthology_category}} -- i.e. one
    RunColumn.category_by_ref-shaped dict per column, keyed by that column's
    header cell (its species label)."""
    table_path = Path(table_path)
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"{table_path} is empty") from None
        if not header or header[0] != "ref_ID":
            raise ValueError(
                f"{table_path} does not look like a `summary`-mode combined orthology TSV "
                "(expected 'ref_ID' as the first column)"
            )
        labels = header[1:]
        category_by_ref_by_label: dict[str, dict[str, str]] = {label: {} for label in labels}
        for row in reader:
            if not row:
                continue
            ref_id = row[0]
            for label, category in zip(labels, row[1:]):
                category_by_ref_by_label[label][ref_id] = category
    return category_by_ref_by_label


def write_combined_orthology(
    ref_ids: list[str], columns: list[RunColumn], out_path: str | Path
) -> int:
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ref_ID"] + [col.label for col in columns])
        for ref_id in ref_ids:
            writer.writerow([ref_id] + [col.category_by_ref.get(ref_id, "NA") for col in columns])
    return len(ref_ids)


def _compute_frequency_distributions(
    ref_ids: list[str], columns: list[RunColumn]
) -> tuple[list[int], list[int]]:
    # loss_counts[i]/multi_counts[i] = number of ref CNEEs that are lost/
    # multi-projection (see LOST_CATEGORIES/MULTI_CATEGORIES) in exactly i
    # of the n_columns samples. Single pass over ref_ids x columns so both
    # distributions share the same per-ref category lookups.
    n_columns = len(columns)
    loss_counts = [0] * (n_columns + 1)
    multi_counts = [0] * (n_columns + 1)
    for ref_id in ref_ids:
        n_lost = 0
        n_multi = 0
        for col in columns:
            kind = _category_kind(col.category_by_ref.get(ref_id))
            if kind == "lost":
                n_lost += 1
            elif kind == "multi":
                n_multi += 1
        loss_counts[n_lost] += 1
        multi_counts[n_multi] += 1
    return loss_counts, multi_counts


def write_frequency_distribution_plots(
    ref_ids: list[str], columns: list[RunColumn], ref_prefix: str, out_dir: Path
) -> tuple[Path, Path]:
    # Per dev_docs/compare_cross_orthology.md's "Summary figures": how many
    # ref CNEEs were lost in exactly 0, 1, 2, ... of the samples, and same
    # for multi-projection -- one bar chart each.
    n_columns = len(columns)
    labels = [str(i) for i in range(n_columns + 1)]
    loss_counts, multi_counts = _compute_frequency_distributions(ref_ids, columns)

    loss_path = Path(out_dir) / f"{ref_prefix}.loss_frequency_distribution.svg"
    loss_path.write_text(
        _bar_svg(
            title="CNEE loss frequency across samples",
            subtitle=f"Reference CNEEs: {len(ref_ids):,}; samples: {n_columns}",
            x_label="Number of samples lost in",
            y_label="Reference CNEEs",
            labels=labels,
            counts=loss_counts,
        ),
        encoding="utf-8",
    )

    multi_path = Path(out_dir) / f"{ref_prefix}.multiproj_frequency_distribution.svg"
    multi_path.write_text(
        _bar_svg(
            title="CNEE multi-projection frequency across samples",
            subtitle=f"Reference CNEEs: {len(ref_ids):,}; samples: {n_columns}",
            x_label="Number of samples multi-projection in",
            y_label="Reference CNEEs",
            labels=labels,
            counts=multi_counts,
        ),
        encoding="utf-8",
    )
    return loss_path, multi_path


def _write_fasta_records(records: list[tuple[str, str]], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), FASTA_LINE_WIDTH):
                handle.write(seq[i : i + FASTA_LINE_WIDTH] + "\n")


def extract_one2one_orthologs(
    ref_ids: list[str],
    ref_coords_by_id: dict[str, tuple[str, int, int]],
    columns: list[RunColumn],
    ref_fasta: str | Path,
    ref_prefix: str,
    out_dir: Path,
    min_length: int,
    max_seqs: int,
    max_lost_frac: float,
    max_multi_frac: float,
) -> tuple[int, int]:
    """Write one <ref_cnee_id>_orthologs.fasta per ref CNEE whose fraction of
    columns categorized as lost (a reference-only category -- see
    LOST_CATEGORIES) is <= max_lost_frac AND whose fraction categorized as
    multi-projection (one2many/many2one/complex/many2many) is <=
    max_multi_frac, at least min_length bp long (by its reference interval,
    not any query projection), and within the first max_seqs such CNEEs in
    ref_ids order. min_length/max_seqs can be set to 0 to disable that bound.
    A qualifying CNEE's file always has the same
    N+1 records -- the reference sequence plus exactly one per column, in
    column order -- so every file lists the same set of species regardless
    of which ones were actually one2one for that CNEE. A column where it
    isn't one2one gets an all-N placeholder record (same length as the
    reference CNEE) instead of a real sequence, not a failure and not an
    omission.

    Processes one column at a time -- its one2one_map (ref->query id, read
    fresh from QUERY_ORTHOLOGY_FILENAME) and its genome accessor are both
    loaded right before that column's turn and released before the next
    column's -- rather than holding every column's data resident at once.
    This matters at real scale for two independent reasons: a .2bit query
    genome's TwoBitAccessor decodes and caches a full chromosome in memory
    the first time anything on it is fetched (see its docstring), and a
    one2one_map can have order-10^5 entries per sample -- either one, held
    for every column simultaneously instead of one at a time, can add up to
    far more memory than a single genome or map, for hundreds of samples.
    Returns (written_count, hard_failed_count) -- the latter counts CNEEs
    dropped entirely because a column claimed one2one but its query
    sequence couldn't actually be located/fetched, a data-consistency
    problem distinct from "not one2one in that column"."""
    n_columns = len(columns)
    qualifying_refs = []
    for ref_id in ref_ids:
        kinds = [_category_kind(col.category_by_ref.get(ref_id)) for col in columns]
        lost_frac = kinds.count("lost") / n_columns
        multi_frac = kinds.count("multi") / n_columns
        if lost_frac <= max_lost_frac and multi_frac <= max_multi_frac:
            qualifying_refs.append(ref_id)
    if not qualifying_refs:
        return 0, 0

    if min_length > 0:
        too_short = [
            ref_id for ref_id in qualifying_refs
            if (ref_coords_by_id[ref_id][2] - ref_coords_by_id[ref_id][1]) < min_length
        ]
        if too_short:
            print(
                f"Dropping {len(too_short)} ref CNEE(s) shorter than --min-length {min_length}",
                file=sys.stderr,
            )
            too_short_set = set(too_short)
            qualifying_refs = [ref_id for ref_id in qualifying_refs if ref_id not in too_short_set]
        if not qualifying_refs:
            return 0, 0

    if max_seqs > 0 and len(qualifying_refs) > max_seqs:
        print(
            f"{len(qualifying_refs)} ref CNEE(s) pass length filtering; capping at --max-seqs "
            f"{max_seqs} (dropping {len(qualifying_refs) - max_seqs})",
            file=sys.stderr,
        )
        qualifying_refs = qualifying_refs[:max_seqs]

    print(
        f"{len(qualifying_refs)} ref CNEE(s) have lost_frac <= {max_lost_frac:.0%} and "
        f"multi_frac <= {max_multi_frac:.0%} across {n_columns} sample(s); extracting FASTA...",
        file=sys.stderr,
    )

    from Bio.Seq import Seq

    # Pass 1: fetch every reference sequence up front (one genome, opened
    # and closed before any query genome is touched), seeding one record
    # list per surviving ref CNEE.
    ref_work_dir = columns[0].run_dir / WORK_DIRNAME
    ref_work_dir.mkdir(parents=True, exist_ok=True)
    records_by_ref: dict[str, list[tuple[str, str]]] = {}
    failed_refs: set[str] = set()
    ref_accessor = open_genome_accessor(ref_fasta, ref_work_dir / "ref_fasta.idx")
    try:
        for ref_id in qualifying_refs:
            ref_chrom, ref_start, ref_end = ref_coords_by_id[ref_id]
            try:
                ref_seq = ref_accessor.fetch(ref_chrom, ref_start, ref_end).upper()
            except (KeyError, ValueError) as exc:
                print(f"warning: skipping {ref_id} (reference fetch failed: {exc})", file=sys.stderr)
                failed_refs.add(ref_id)
                continue
            records_by_ref[ref_id] = [(f"{ref_prefix} {ref_id} {ref_chrom}:{ref_start}-{ref_end}", ref_seq)]
    finally:
        ref_accessor.close()

    # Pass 2: one query genome at a time. A column where a ref CNEE isn't
    # one2one gets an all-N placeholder record instead of a real sequence
    # (expected, not a failure), so every output file keeps the same
    # species in the same order regardless of which ones matched. A ref
    # CNEE that hard-fails at column K -- one2one per category_by_ref, but
    # its query sequence can't be found/fetched -- is marked failed and
    # skipped (without a further warning) by every later column, matching
    # the original per-ref "stop at the first failure" behavior, just
    # reached via a flag instead of an inner-loop break, since columns are
    # now the outer loop.
    n_placeholder = 0
    for col_index, col in enumerate(columns, start=1):
        # Read on demand rather than upfront in collect_cross_run_orthology
        # -- see RunColumn's comment -- and let it go out of scope (freed)
        # once this column's iteration ends, same as the genome accessor.
        _, one2one_map = _read_query_orthology_file(col.run_dir / QUERY_ORTHOLOGY_FILENAME, True)

        # Restrict this run's query_cnees.bed lookup to only the ids we
        # still need -- that file can be huge, but the wanted set is known.
        wanted_ids = {
            one2one_map[ref_id]
            for ref_id in qualifying_refs
            if ref_id not in failed_refs
            and col.category_by_ref.get(ref_id) == "one2one"
            and ref_id in one2one_map
        }
        coords = _read_query_bed_coords(col.run_dir / QUERY_BED_FILENAME, wanted_ids)

        work_dir = col.run_dir / WORK_DIRNAME
        work_dir.mkdir(parents=True, exist_ok=True)
        accessor = open_genome_accessor(col.genome_path, work_dir / "query_fasta.idx")
        try:
            for ref_id in qualifying_refs:
                if ref_id in failed_refs:
                    continue
                if col.category_by_ref.get(ref_id) != "one2one":
                    n_placeholder += 1
                    ref_start, ref_end = ref_coords_by_id[ref_id][1], ref_coords_by_id[ref_id][2]
                    records_by_ref[ref_id].append((f"{col.label} NA NA", "N" * (ref_end - ref_start)))
                    continue
                query_id = one2one_map.get(ref_id)
                if query_id is None or query_id not in coords:
                    print(
                        f"warning: skipping {ref_id} ({col.label} claims one2one but has no "
                        "query BED coordinates for its ortholog)",
                        file=sys.stderr,
                    )
                    failed_refs.add(ref_id)
                    records_by_ref.pop(ref_id, None)
                    continue
                chrom, start, end, strand = coords[query_id]
                try:
                    seq = accessor.fetch(chrom, start, end)
                except (KeyError, ValueError) as exc:
                    print(f"warning: skipping {ref_id} ({col.label} fetch failed: {exc})", file=sys.stderr)
                    failed_refs.add(ref_id)
                    records_by_ref.pop(ref_id, None)
                    continue
                if strand == "-":
                    seq = str(Seq(seq).reverse_complement())
                records_by_ref[ref_id].append(
                    (f"{col.label} {query_id} {chrom}:{start}-{end}({strand})", seq.upper())
                )
        finally:
            # Closing releases this column's genome from memory (for a
            # .2bit accessor, its decoded-chromosome cache) before the next
            # column's genome is opened.
            accessor.close()

        if col_index % 50 == 0 or col_index == len(columns):
            print(f"  ...{col_index}/{len(columns)} samples processed", file=sys.stderr)

    if n_placeholder:
        print(
            f"Filled {n_placeholder} per-sample sequence(s) with an all-N placeholder "
            "for samples that weren't one2one for that CNEE",
            file=sys.stderr,
        )

    # Pass 3: write out every ref CNEE that didn't hard-fail, batched by its
    # position in qualifying_refs (unchanged by which ones failed).
    fasta_root = out_dir / FASTA_SUBDIR_NAME
    written = 0
    for index, ref_id in enumerate(qualifying_refs):
        if ref_id in failed_refs:
            continue
        batch_dir = fasta_root / f"batch{index // FASTA_BATCH_SIZE + 1}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        _write_fasta_records(records_by_ref[ref_id], batch_dir / f"{_sanitize_filename(ref_id)}_orthologs.fasta")
        written += 1

    return written, len(failed_refs)


def _run_summary_mode(args: argparse.Namespace) -> None:
    run_dirs = _read_run_dirs(args.query_dirs)
    if not run_dirs:
        raise SystemExit(f"{args.query_dirs} lists no output directories")

    columns = collect_cross_run_orthology(run_dirs, build_fasta_data=False)
    if not columns:
        raise SystemExit("No usable run directories found; nothing to write")

    ref_ids = [interval.cnee_uid for interval in parse_bed(args.ref_bed)]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.ref_prefix}_combined_orthology.tsv"
    row_count = write_combined_orthology(ref_ids, columns, out_path)
    print(f"Wrote {row_count} ref CNEE rows x {len(columns)} species -> {out_path}")

    loss_path, multi_path = write_frequency_distribution_plots(ref_ids, columns, args.ref_prefix, args.out_dir)
    print(f"Wrote loss frequency distribution -> {loss_path}")
    print(f"Wrote multi-projection frequency distribution -> {multi_path}")


def _run_sequence_mode(args: argparse.Namespace) -> None:
    if not 0.0 <= args.max_lost_frac <= 1.0:
        raise SystemExit(f"--max-lost-frac must be between 0 and 1 (got {args.max_lost_frac})")
    if not 0.0 <= args.max_multi_frac <= 1.0:
        raise SystemExit(f"--max-multi-frac must be between 0 and 1 (got {args.max_multi_frac})")

    run_dirs = _read_run_dirs(args.query_dirs)
    if not run_dirs:
        raise SystemExit(f"{args.query_dirs} lists no output directories")

    columns = collect_cross_run_orthology(run_dirs, build_fasta_data=True, require_ref_ortho_status=False)
    if not columns:
        raise SystemExit("No usable run directories found; nothing to write")

    missing_genome = [str(col.run_dir) for col in columns if col.genome_path is None]
    if missing_genome:
        raise SystemExit(
            "--query-dirs is missing a genome path (2nd column) for: "
            f"{', '.join(missing_genome)} (required for FASTA extraction)"
        )

    per_label_categories = _read_combined_orthology_table(args.orthology_table)
    for col in columns:
        if col.label not in per_label_categories:
            raise SystemExit(
                f"{args.orthology_table} has no column for run {col.run_dir} (label {col.label!r}); "
                "make sure --query-dirs is the same sample sheet (same run dirs/labels) used to "
                "build this orthology table with `summary` mode"
            )
        col.category_by_ref = per_label_categories[col.label]

    ref_intervals = parse_bed(args.ref_bed)
    ref_ids = [interval.cnee_uid for interval in ref_intervals]
    ref_coords_by_id = {interval.cnee_uid: (interval.chrom, interval.start, interval.end) for interval in ref_intervals}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    written, skipped = extract_one2one_orthologs(
        ref_ids, ref_coords_by_id, columns, args.ref_genome, args.ref_prefix, args.out_dir,
        args.min_length, args.max_seqs, args.max_lost_frac, args.max_multi_frac,
    )
    skip_note = f" ({skipped} skipped)" if skipped else ""
    print(
        f"Wrote {written} ortholog FASTA file(s){skip_note} -> "
        f"{args.out_dir / FASTA_SUBDIR_NAME}"
    )


@dataclass
class TreeScaffold:
    """Everything needed to compute Fritz & Purvis's D statistic on a fixed
    tree topology, precomputed once (not re-derived per CNEE). Node "idx"
    values are toytree's own node ids -- tips are idx 0..n_tips-1, matching
    the row order of toytree.pcm.simulate_continuous_bm's output."""

    tip_idx_by_label: dict[str, int]
    postorder_internal: list[int]  # internal node idx, children before parents
    children: dict[int, list[int]]  # internal node idx -> child idx list
    edges: list[tuple[int, int]]  # (child_idx, parent_idx), every node but the root
    n_tips: int


def _build_tree_scaffold(tree) -> TreeScaffold:
    tip_idx_by_label: dict[str, int] = {}
    postorder_internal: list[int] = []
    children: dict[int, list[int]] = {}
    edges: list[tuple[int, int]] = []
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            tip_idx_by_label[node.name] = node.idx
        else:
            postorder_internal.append(node.idx)
            children[node.idx] = [child.idx for child in node.children]
        if node.up is not None:
            edges.append((node.idx, node.up.idx))
    return TreeScaffold(tip_idx_by_label, postorder_internal, children, edges, len(tip_idx_by_label))


def _compute_d(tip_values: dict[int, float], scaffold: TreeScaffold) -> float:
    # Fritz & Purvis's D: propagate tip values up as each internal node's
    # (unweighted) mean of its children, then sum |child - parent| over
    # every edge in the tree. Used identically for d(obs), each d(random)
    # permutation, and each d(brownian) simulation -- only tip_values differ.
    node_value = dict(tip_values)
    for idx in scaffold.postorder_internal:
        kids = scaffold.children[idx]
        node_value[idx] = sum(node_value[k] for k in kids) / len(kids)
    return sum(abs(node_value[child] - node_value[parent]) for child, parent in scaffold.edges)


def _simulate_random_d_values(n_lost: int, scaffold: TreeScaffold, n_perms: int, rng) -> list[float]:
    # Null distribution 1: no phylogenetic structure -- shuffle which tips
    # are "lost" at random, keeping the same observed count of losses.
    import numpy as np

    tip_idxs = list(scaffold.tip_idx_by_label.values())
    base = np.array([0.0] * n_lost + [1.0] * (scaffold.n_tips - n_lost))
    values = []
    for _ in range(n_perms):
        shuffled = rng.permutation(base)
        values.append(_compute_d(dict(zip(tip_idxs, shuffled)), scaffold))
    return values


def _simulate_brownian_d_values(
    tree, n_lost: int, scaffold: TreeScaffold, n_perms: int, rng
) -> list[float]:
    # Null distribution 2: trait evolves under Brownian motion, then
    # thresholded to binary at the quantile matching the OBSERVED prevalence
    # (n_lost losses) -- not a fixed median -- so the null reflects what a
    # trait this rare/common would look like under pure BM. One call
    # simulates all n_perms replicates at once (one column each).
    import numpy as np

    seed = int(rng.integers(0, 2**32 - 1))
    sim = tree.pcm.simulate_continuous_bm(rates=[1.0] * n_perms, tips_only=True, seed=seed)
    tip_idxs = list(scaffold.tip_idx_by_label.values())
    values = []
    for column in sim.columns:
        tip_continuous = sim[column].to_numpy()
        lost_positions = np.argsort(tip_continuous)[:n_lost]
        binary = np.ones(scaffold.n_tips)
        binary[lost_positions] = 0.0
        values.append(_compute_d(dict(zip(tip_idxs, binary)), scaffold))
    return values


def run_phylo_d_analysis(
    orthology_table: Path,
    tree_path: Path,
    min_lost: int,
    n_perms: int,
    seed: int,
):
    import numpy as np
    import pandas as pd
    import toytree

    tree = toytree.tree(str(tree_path))
    scaffold = _build_tree_scaffold(tree)
    tip_labels = list(scaffold.tip_idx_by_label.keys())

    with Path(orthology_table).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
    if not header or header[0] != "ref_ID":
        raise ValueError(
            f"{orthology_table} does not look like a `summary`-mode combined orthology TSV "
            "(expected 'ref_ID' as the first column)"
        )
    table_species = header[1:]
    missing_tips = sorted(set(tip_labels) - set(table_species))
    if missing_tips:
        raise SystemExit(
            f"Tree tip label(s) not found as columns in {orthology_table}: {', '.join(missing_tips)} "
            "-- tree tip names must exactly match the species labels `summary` mode assigned"
        )
    species_column_index = {label: header.index(label) for label in tip_labels}

    rng = np.random.default_rng(seed)
    null_cache: dict[int, tuple[float, float]] = {}
    rows: list[tuple[str, int, float]] = []
    n_total = n_dropped_na = n_invariant = n_below_min_lost = 0

    with Path(orthology_table).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            n_total += 1
            ref_id = row[0]
            tip_values: dict[int, float] = {}
            n_lost = 0
            has_na = False
            for label in tip_labels:
                category = row[species_column_index[label]]
                if not category or category == "NA":
                    has_na = True
                    break
                is_lost = _category_kind(category) == "lost"
                tip_values[scaffold.tip_idx_by_label[label]] = 0.0 if is_lost else 1.0
                n_lost += int(is_lost)
            if has_na:
                n_dropped_na += 1
                continue
            if n_lost == 0 or n_lost == scaffold.n_tips:
                n_invariant += 1
                continue
            if n_lost < min_lost:
                n_below_min_lost += 1
                continue

            d_obs = _compute_d(tip_values, scaffold)
            if n_lost not in null_cache:
                random_vals = _simulate_random_d_values(n_lost, scaffold, n_perms, rng)
                brownian_vals = _simulate_brownian_d_values(tree, n_lost, scaffold, n_perms, rng)
                null_cache[n_lost] = (float(np.mean(random_vals)), float(np.mean(brownian_vals)))
            mean_random, mean_brownian = null_cache[n_lost]
            denom = mean_random - mean_brownian
            d_stat = (d_obs - mean_brownian) / denom if denom != 0 else float("nan")
            rows.append((ref_id, n_lost, d_stat))

    print(
        f"{n_total} ref CNEEs read; {n_dropped_na} dropped (missing data), "
        f"{n_invariant} skipped (invariant: lost in 0 or all {scaffold.n_tips} species), "
        f"{n_below_min_lost} skipped (lost in fewer than --min-lost {min_lost} species); "
        f"{len(rows)} D values computed across {len(null_cache)} distinct lost-count group(s)",
        file=sys.stderr,
    )
    return pd.DataFrame(rows, columns=["ref_cnee_id", "n_lost", "D"])


def write_phylo_d_density_plot(d_values, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(d_values, bins=30, color="#1f6feb", edgecolor="white")
    ax.axvline(0, color="#0e9594", linestyle="--", linewidth=1, label="D=0 (clumped, BM-like)")
    ax.axvline(1, color="#c2410c", linestyle="--", linewidth=1, label="D=1 (random)")
    ax.set_xlabel("D statistic")
    ax.set_ylabel("Reference CNEEs")
    ax.set_title(f"Phylogenetic signal of CNEE loss (n={len(d_values):,})")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _run_cluster_mode(args: argparse.Namespace) -> None:
    result = run_phylo_d_analysis(
        args.orthology_table, args.tree, args.min_lost, args.n_perms, args.seed,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.ref_prefix}_phylo_d.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(result)} D-statistic rows -> {out_path}")

    valid = result["D"].dropna()
    if len(valid):
        fig_path = args.out_dir / "figures" / f"{args.ref_prefix}.phylo_d_density.png"
        write_phylo_d_density_plot(valid.to_numpy(), fig_path)
        print(f"Wrote D-statistic density plot ({len(valid)} finite values) -> {fig_path}")
    else:
        print("No finite D values to plot (every candidate CNEE had a degenerate null distribution)", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    query_dirs_help = (
        "Text file listing one toga4cnee output directory per line, optionally "
        "with tab-separated columns for that run's query genome (FASTA/2bit, "
        "required for FASTA extraction) and a manual species label"
    )
    ref_prefix_help = (
        "Label for the reference set: used both for the '<ref-prefix>_combined_orthology.tsv' "
        "filename and as the reference sequence's FASTA header"
    )
    out_dir_help = "Directory to write output into [default: current directory]"

    summary_parser = subparsers.add_parser(
        "summary",
        help="Build the combined per-species orthology-category TSV",
        description="Build one combined TSV of per-reference-CNEE orthology category across samples.",
    )
    summary_parser.add_argument("--query-dirs", required=True, type=Path, help=query_dirs_help)
    summary_parser.add_argument(
        "--ref-bed", required=True, type=Path, help="Reference CNEE BED3+ file used for every run"
    )
    summary_parser.add_argument("--ref-prefix", required=True, help=ref_prefix_help)
    summary_parser.add_argument("--out-dir", type=Path, default=Path("."), help=out_dir_help)
    summary_parser.set_defaults(func=_run_summary_mode)

    sequence_parser = subparsers.add_parser(
        "sequence",
        help="Extract FASTA sequence for shared one2one orthologs",
        description="Extract one <ref_cnee_id>_orthologs.fasta per ref CNEE that's one2one in enough samples.",
    )
    sequence_parser.add_argument("--query-dirs", required=True, type=Path, help=query_dirs_help)
    sequence_parser.add_argument(
        "--orthology-table",
        required=True,
        type=Path,
        help=(
            "The '<ref-prefix>_combined_orthology.tsv' built by `summary` mode -- must come "
            "from the same --query-dirs sample sheet (same run dirs/labels), since categories "
            "are looked up by matching each run's label to that table's column header"
        ),
    )
    sequence_parser.add_argument(
        "--ref-bed", required=True, type=Path, help="Reference CNEE BED3+ file used for every run"
    )
    sequence_parser.add_argument(
        "--ref-genome", required=True, type=Path, help="Reference genome (FASTA/2bit)"
    )
    sequence_parser.add_argument("--ref-prefix", required=True, help=ref_prefix_help)
    sequence_parser.add_argument("--out-dir", type=Path, default=Path("."), help=out_dir_help)
    sequence_parser.add_argument(
        "--min-length",
        type=int,
        default=40,
        help=(
            "Minimum reference CNEE length (bp) to include in FASTA extraction, "
            "decided by the reference interval's own length; 0 disables this filter "
            "[default: 40]"
        ),
    )
    sequence_parser.add_argument(
        "--max-seqs",
        type=int,
        default=300000,
        help=(
            "Maximum number of one2one ortholog FASTA files to write; 0 disables this cap "
            "[default: 300000]"
        ),
    )
    sequence_parser.add_argument(
        "--max-lost-frac",
        type=float,
        default=0.3,
        help=(
            "Maximum fraction of samples a ref CNEE can be lost in (a reference-only "
            "category -- unaligned/probable_loss/assembly_gap_loss/uncertain_loss/"
            "missing-loss) and still be included in FASTA extraction [default: 0.3]"
        ),
    )
    sequence_parser.add_argument(
        "--max-multi-frac",
        type=float,
        default=0.1,
        help=(
            "Maximum fraction of samples a ref CNEE can be multi-projection in "
            "(one2many/many2one/complex) and still be included in FASTA extraction "
            "[default: 0.1]"
        ),
    )
    sequence_parser.set_defaults(func=_run_sequence_mode)

    cluster_parser = subparsers.add_parser(
        "cluster",
        help="Phylogenetic signal (Fritz & Purvis D) of CNEE loss across species",
        description=(
            "For each qualifying ref CNEE, compute the D statistic of phylogenetic signal "
            "(Fritz & Purvis 2010) for its present/lost pattern across species on a tree."
        ),
    )
    cluster_parser.add_argument(
        "--orthology-table",
        required=True,
        type=Path,
        help="The '<ref-prefix>_combined_orthology.tsv' built by `summary` mode",
    )
    cluster_parser.add_argument(
        "--tree",
        required=True,
        type=Path,
        help=(
            "Phylogenetic tree of species in Newick format; tip labels must exactly match "
            "the species labels used as columns in --orthology-table"
        ),
    )
    cluster_parser.add_argument("--ref-prefix", required=True, help=ref_prefix_help)
    cluster_parser.add_argument("--out-dir", type=Path, default=Path("."), help=out_dir_help)
    cluster_parser.add_argument(
        "--min-lost",
        type=int,
        default=3,
        help=(
            "Only compute D for a ref CNEE lost in at least this many species (it's always "
            "skipped if lost in 0 or all species, since D is undefined for an invariant trait) "
            "[default: 3]"
        ),
    )
    cluster_parser.add_argument(
        "--n-perms",
        type=int,
        default=10,
        help=(
            "Number of permutations/simulations per distinct lost-species-count group for each "
            "of the two null distributions (random, Brownian motion) [default: 10]"
        ),
    )
    cluster_parser.add_argument(
        "--seed", type=int, default=0, help="Random seed, for reproducible permutations/simulations [default: 0]"
    )
    cluster_parser.set_defaults(func=_run_cluster_mode)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
