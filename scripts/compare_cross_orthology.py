#!/usr/bin/env python
"""Collect per-reference-CNEE orthology status across multiple pipeline runs.

Takes a list of `toga4cnee run` output directories -- the same reference BED
lifted to different query species -- and builds one combined TSV where each
row is a reference CNEE and each column is a query species' final
orthology_category for that CNEE, read from each directory's
ref_ortho_status.tsv (see dev_docs/compare_cross_orthology.md).

Each run's column label is chosen in order of preference:
  1. A manual label given as an optional third, tab-separated column in
     --query-dirs (one "outdir<TAB>genome<TAB>label" line per run).
  2. Otherwise, the query_id_prefix baked into that run's
     query_ortho_status.tsv's query_cnee_id values (the text before the
     first "_") -- unless that comes out as "q", the pipeline's generic
     default prefix and uninformative across runs, or the file has no
     usable (non-"NA") query id at all.
  3. Otherwise, that output directory's own folder name.

By default (disable with --no-fasta), also extracts a FASTA of orthologous
sequences for every reference CNEE that is one2one in *every* listed run:
one combined "<ref_cnee_id>_orthologs.fasta" per such CNEE (reference
sequence plus each sample's query ortholog), written under
<out-dir>/one2one_ortholog_fasta/batchN/ in batches of 10000 files, so a
large CNEE set doesn't dump hundreds of thousands of files into one
directory. This needs --ref-genome plus a genome path in --query-dirs'
second column for every run, and reuses each run's own byte-offset FASTA
index (<outdir>/work/query_fasta.idx) built by the original `toga4cnee run`
rather than re-indexing -- so point the genome column at the exact same
FASTA/2bit path that run was originally given.

Example:
    PYTHONPATH=src python scripts/compare_cross_orthology.py \\
        --query-dirs query_dirs.txt \\
        --ref-bed dummyRef.cnees.bed \\
        --ref-genome dummyRef.fasta \\
        --ref-prefix dummyRef \\
        --out-dir .

Where query_dirs.txt looks like:
    out/dummyRef_vs_sppA   genomes/sppA.2bit       sppA
    out/dummyRef_vs_sppB   genomes/sppB.fasta
    out/dummyRef_vs_sppC   genomes/sppC.2bit       sppC
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
) -> tuple[int, int]:
    """Write one <ref_cnee_id>_orthologs.fasta per ref CNEE that's one2one
    in every column. Returns (written_count, skipped_count)."""
    all_one2one = [
        ref_id for ref_id in ref_ids if all(col.category_by_ref.get(ref_id) == "one2one" for col in columns)
    ]
    if not all_one2one:
        return 0, 0

    print(
        f"{len(all_one2one)} ref CNEE(s) are one2one in every sample; extracting FASTA...",
        file=sys.stderr,
    )

    # Restrict each run's query_cnees.bed lookup to only the ids we need --
    # that file can be huge, but the wanted set is known up front.
    coords_by_column = []
    for col in columns:
        wanted_ids = {col.one2one_map[ref_id] for ref_id in all_one2one if ref_id in col.one2one_map}
        coords_by_column.append(_read_query_bed_coords(col.run_dir / QUERY_BED_FILENAME, wanted_ids))

    # Reuse each run's own pipeline-built index instead of re-indexing.
    ref_work_dir = columns[0].run_dir / WORK_DIRNAME
    ref_work_dir.mkdir(parents=True, exist_ok=True)
    ref_accessor = open_genome_accessor(ref_fasta, ref_work_dir / "ref_fasta.idx")
    query_accessors = []
    for col in columns:
        work_dir = col.run_dir / WORK_DIRNAME
        work_dir.mkdir(parents=True, exist_ok=True)
        query_accessors.append(open_genome_accessor(col.genome_path, work_dir / "query_fasta.idx"))

    from Bio.Seq import Seq

    fasta_root = out_dir / FASTA_SUBDIR_NAME
    written = 0
    skipped = 0
    try:
        for index, ref_id in enumerate(all_one2one):
            if index and index % 50000 == 0:
                print(f"  ...{index}/{len(all_one2one)} processed", file=sys.stderr)

            ref_chrom, ref_start, ref_end = ref_coords_by_id[ref_id]
            try:
                ref_seq = ref_accessor.fetch(ref_chrom, ref_start, ref_end).upper()
            except (KeyError, ValueError) as exc:
                print(f"warning: skipping {ref_id} (reference fetch failed: {exc})", file=sys.stderr)
                skipped += 1
                continue

            records = [(f"{ref_prefix} {ref_id} {ref_chrom}:{ref_start}-{ref_end}", ref_seq)]
            ok = True
            for col, coords, accessor in zip(columns, coords_by_column, query_accessors):
                query_id = col.one2one_map.get(ref_id)
                if query_id is None or query_id not in coords:
                    print(
                        f"warning: skipping {ref_id} ({col.label} has no query BED coordinates "
                        "for its one2one ortholog)",
                        file=sys.stderr,
                    )
                    ok = False
                    break
                chrom, start, end, strand = coords[query_id]
                try:
                    seq = accessor.fetch(chrom, start, end)
                except (KeyError, ValueError) as exc:
                    print(f"warning: skipping {ref_id} ({col.label} fetch failed: {exc})", file=sys.stderr)
                    ok = False
                    break
                if strand == "-":
                    seq = str(Seq(seq).reverse_complement())
                records.append((f"{col.label} {query_id} {chrom}:{start}-{end}({strand})", seq.upper()))
            if not ok:
                skipped += 1
                continue

            batch_dir = fasta_root / f"batch{index // FASTA_BATCH_SIZE + 1}"
            batch_dir.mkdir(parents=True, exist_ok=True)
            _write_fasta_records(records, batch_dir / f"{_sanitize_filename(ref_id)}_orthologs.fasta")
            written += 1
    finally:
        ref_accessor.close()
        for accessor in query_accessors:
            accessor.close()

    return written, skipped


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
            ref_ids, ref_coords_by_id, columns, args.ref_genome, args.ref_prefix, args.out_dir
        )
        skip_note = f" ({skipped} skipped)" if skipped else ""
        print(
            f"Wrote {written} ortholog FASTA file(s){skip_note} -> "
            f"{args.out_dir / FASTA_SUBDIR_NAME}"
        )


if __name__ == "__main__":
    main()
