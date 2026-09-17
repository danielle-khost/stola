#!/usr/bin/env python
"""Collect per-reference-CNEE orthology status across multiple pipeline runs.

Takes a list of `toga4cnee run` output directories -- the same reference BED
lifted to different query species -- and builds one combined TSV where each
row is a reference CNEE and each column is a query species' final
orthology_category for that CNEE, read from each directory's
ref_ortho_status.tsv 
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from toga4cnee.bed import parse_bed
from toga4cnee.genome_io import open_genome_accessor
from toga4cnee.orthology import ORTHOLOGY_COLUMNS, REF_CNEE_COLUMNS
from toga4cnee.workflow import (
    QUERY_BED_FILENAME,
    QUERY_ORTHOLOGY_FILENAME,
    REF_ORTHOLOGY_FILENAME,
    WORK_DIRNAME,
)

GENERIC_LABEL = "q"  # DEFAULT_QUERY_ID_PREFIX ("q_") stripped of its trailing "_".
FASTA_SUBDIR_NAME = "one2one_ortholog_fasta"
FASTA_BATCH_SIZE = 10000
FASTA_LINE_WIDTH = 60


@dataclass
class RunColumn:
    label: str
    category_by_ref: dict[str, str]
    run_dir: Path
    genome_path: Path | None
    one2one_map: dict[str, str] = field(default_factory=dict)  # ref_cnee_id -> query_cnee_id


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
    run_dirs: list[tuple[Path, Path | None, str | None]], build_fasta_data: bool
) -> list[RunColumn]:
    columns: list[RunColumn] = []
    used_labels: set[str] = set()
    for run_dir, genome_path, manual_label in run_dirs:
        ref_orthology_path = run_dir / REF_ORTHOLOGY_FILENAME
        query_orthology_path = run_dir / QUERY_ORTHOLOGY_FILENAME
        if not ref_orthology_path.is_file():
            print(f"warning: skipping {run_dir} (missing {REF_ORTHOLOGY_FILENAME})", file=sys.stderr)
            continue

        one2one_map: dict[str, str] = {}
        if build_fasta_data:
            # FASTA mode needs this file's actual data (not just a peek), so
            # a manual label can no longer substitute for its absence.
            if not query_orthology_path.is_file():
                print(
                    f"warning: skipping {run_dir} (missing {QUERY_ORTHOLOGY_FILENAME}, "
                    "needed to locate one2one query orthologs for FASTA extraction)",
                    file=sys.stderr,
                )
                continue
            peeked_label, one2one_map = _read_query_orthology_file(query_orthology_path, True)
            label = manual_label or (
                peeked_label if peeked_label and peeked_label != GENERIC_LABEL else run_dir.resolve().name
            )
        elif manual_label is not None:
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
        columns.append(
            RunColumn(
                label=label,
                category_by_ref=_read_ref_orthology(ref_orthology_path),
                run_dir=run_dir,
                genome_path=genome_path,
                one2one_map=one2one_map,
            )
        )

    return columns


def write_combined_orthology(
    ref_ids: list[str], columns: list[RunColumn], out_path: str | Path
) -> int:
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ref_ID"] + [col.label for col in columns])
        for ref_id in ref_ids:
            writer.writerow([ref_id] + [col.category_by_ref.get(ref_id, "NA") for col in columns])
    return len(ref_ids)


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
) -> tuple[int, int]:
    """Write one <ref_cnee_id>_orthologs.fasta per ref CNEE that's one2one
    in every column, at least min_length bp long (by its reference interval,
    not any query projection), and within the first max_seqs such CNEEs in
    ref_ids order. Either bound can be set to 0 to disable it.

    Processes one query genome (column) at a time -- fetch every sequence
    that column contributes, then close its accessor before opening the
    next -- rather than holding every column's genome accessor open at
    once. This matters for .2bit query genomes: TwoBitAccessor decodes and
    caches a full chromosome in memory the first time anything on it is
    fetched (see its docstring), so holding all columns open simultaneously
    means up to len(columns) genomes' worth of decoded sequence resident at
    once -- for hundreds of samples that can far exceed available RAM. One
    genome open at a time bounds that to a single genome regardless of how
    many columns there are. Returns (written_count, skipped_count)."""
    all_one2one = [
        ref_id for ref_id in ref_ids if all(col.category_by_ref.get(ref_id) == "one2one" for col in columns)
    ]
    if not all_one2one:
        return 0, 0

    if min_length > 0:
        too_short = [
            ref_id for ref_id in all_one2one
            if (ref_coords_by_id[ref_id][2] - ref_coords_by_id[ref_id][1]) < min_length
        ]
        if too_short:
            print(
                f"Dropping {len(too_short)} ref CNEE(s) shorter than --min-length {min_length}",
                file=sys.stderr,
            )
            too_short_set = set(too_short)
            all_one2one = [ref_id for ref_id in all_one2one if ref_id not in too_short_set]
        if not all_one2one:
            return 0, 0

    if max_seqs > 0 and len(all_one2one) > max_seqs:
        print(
            f"{len(all_one2one)} ref CNEE(s) pass length filtering; capping at --max-seqs "
            f"{max_seqs} (dropping {len(all_one2one) - max_seqs})",
            file=sys.stderr,
        )
        all_one2one = all_one2one[:max_seqs]

    print(
        f"{len(all_one2one)} ref CNEE(s) are one2one in every sample; extracting FASTA...",
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
        for ref_id in all_one2one:
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

    # Pass 2: one query genome at a time. A ref CNEE that fails at column K
    # is marked failed and skipped (without a warning) by every later
    # column, matching the original per-ref "stop at the first failure"
    # behavior -- just reached via a flag instead of an inner-loop break,
    # since columns are now the outer loop.
    for col_index, col in enumerate(columns, start=1):
        # Restrict this run's query_cnees.bed lookup to only the ids we
        # still need -- that file can be huge, but the wanted set is known.
        wanted_ids = {
            col.one2one_map[ref_id]
            for ref_id in all_one2one
            if ref_id not in failed_refs and ref_id in col.one2one_map
        }
        coords = _read_query_bed_coords(col.run_dir / QUERY_BED_FILENAME, wanted_ids)

        work_dir = col.run_dir / WORK_DIRNAME
        work_dir.mkdir(parents=True, exist_ok=True)
        accessor = open_genome_accessor(col.genome_path, work_dir / "query_fasta.idx")
        try:
            for ref_id in all_one2one:
                if ref_id in failed_refs:
                    continue
                query_id = col.one2one_map.get(ref_id)
                if query_id is None or query_id not in coords:
                    print(
                        f"warning: skipping {ref_id} ({col.label} has no query BED coordinates "
                        "for its one2one ortholog)",
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

    # Pass 3: write out every ref CNEE that survived every column, batched
    # by its position in all_one2one (unchanged by which ones failed).
    fasta_root = out_dir / FASTA_SUBDIR_NAME
    written = 0
    for index, ref_id in enumerate(all_one2one):
        if ref_id in failed_refs:
            continue
        batch_dir = fasta_root / f"batch{index // FASTA_BATCH_SIZE + 1}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        _write_fasta_records(records_by_ref[ref_id], batch_dir / f"{_sanitize_filename(ref_id)}_orthologs.fasta")
        written += 1

    return written, len(failed_refs)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--query-dirs",
        required=True,
        type=Path,
        help=(
            "Text file listing one toga4cnee output directory per line, optionally "
            "with tab-separated columns for that run's query genome (FASTA/2bit, "
            "required for FASTA extraction) and a manual species label"
        ),
    )
    parser.add_argument(
        "--ref-bed", required=True, type=Path, help="Reference CNEE BED3+ file used for every run"
    )
    parser.add_argument(
        "--ref-genome",
        type=Path,
        default=None,
        help="Reference genome (FASTA/2bit); required unless --no-fasta is given",
    )
    parser.add_argument(
        "--ref-prefix",
        required=True,
        help=(
            "Label for the reference set: used both for the '<ref-prefix>_combined_orthology.tsv' "
            "filename and as the reference sequence's FASTA header"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("."),
        help="Directory to write output into [default: current directory]",
    )
    parser.add_argument(
        "--no-fasta",
        action="store_true",
        help="Only write the combined orthology TSV; skip one2one FASTA extraction",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=40,
        help=(
            "Minimum reference CNEE length (bp) to include in FASTA extraction, "
            "decided by the reference interval's own length; 0 disables this filter "
            "[default: 40]"
        ),
    )
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=100000,
        help=(
            "Maximum number of one2one ortholog FASTA files to write; 0 disables this cap "
            "[default: 100000]"
        ),
    )
    args = parser.parse_args(argv)

    fasta_enabled = not args.no_fasta
    if fasta_enabled and args.ref_genome is None:
        raise SystemExit("--ref-genome is required unless --no-fasta is given")

    run_dirs = _read_run_dirs(args.query_dirs)
    if not run_dirs:
        raise SystemExit(f"{args.query_dirs} lists no output directories")

    columns = collect_cross_run_orthology(run_dirs, fasta_enabled)
    if not columns:
        raise SystemExit("No usable run directories found; nothing to write")

    if fasta_enabled:
        missing_genome = [str(col.run_dir) for col in columns if col.genome_path is None]
        if missing_genome:
            raise SystemExit(
                "--query-dirs is missing a genome path (2nd column) for: "
                f"{', '.join(missing_genome)} (required for FASTA extraction; "
                "pass --no-fasta to skip it)"
            )

    ref_intervals = parse_bed(args.ref_bed)
    ref_ids = [interval.cnee_uid for interval in ref_intervals]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.ref_prefix}_combined_orthology.tsv"
    row_count = write_combined_orthology(ref_ids, columns, out_path)
    print(f"Wrote {row_count} ref CNEE rows x {len(columns)} species -> {out_path}")

    if fasta_enabled:
        ref_coords_by_id = {interval.cnee_uid: (interval.chrom, interval.start, interval.end) for interval in ref_intervals}
        written, skipped = extract_one2one_orthologs(
            ref_ids, ref_coords_by_id, columns, args.ref_genome, args.ref_prefix, args.out_dir,
            args.min_length, args.max_seqs,
        )
        skip_note = f" ({skipped} skipped)" if skipped else ""
        print(
            f"Wrote {written} ortholog FASTA file(s){skip_note} -> "
            f"{args.out_dir / FASTA_SUBDIR_NAME}"
        )


if __name__ == "__main__":
    main()
