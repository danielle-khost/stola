import contextlib
import csv
import gzip
import io
import random
import struct
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from stola.candidates import write_candidates
from stola.index import build_index, load_index
from stola.overlap import calculate_overlap_rows, write_overlap_stats
from stola.sequence_compare import SequenceComparator
from stola.summary import SUMMARY_COLUMNS, summarize_overlap_stats
from stola.workflow import run_workflow


def _has_networkx() -> bool:
    try:
        import networkx  # noqa: F401
    except ImportError:
        return False
    return True


CHAIN_TEXT = "\n".join(
    [
        "chain 100 chr1 1000 + 100 300 chrQ 1000 + 500 650 1",
        "50 50 0",
        "100",
        "",
        "chain 200 chr1 1000 + 110 240 chrQ2 1000 - 100 230 2",
        "130",
        "",
        "chain 50 chr2 1000 + 100 300 chrQ3 1000 + 400 500 3",
        "50 100 0",
        "50",
        "",
    ]
)


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    chain = tmp_path / "test.chain"
    bed = tmp_path / "cnees.bed"
    chain.write_text(CHAIN_TEXT, encoding="utf-8")
    bed.write_text(
        "\n".join(
            [
                "chr1\t120\t240\tCNEE_PARTIAL",
                "chr2\t175\t225\tCNEE_GAPPED",
                "chr3\t10\t20\tCNEE_NO_CHAIN",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return chain, bed


def _write_fasta(path: Path, records: dict[str, str], width: int = 60) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def _write_fasta_inputs(tmp_path: Path) -> tuple[Path, Path]:
    # Homopolymer ref/query built from a complementary base pair (A/T), so a
    # minus-strand block's reverse-complemented query fetch is still a clean
    # 100% match against the ref regardless of phase/offset -- these tests
    # care about gap fraction (from the chain itself), not sequence content.
    ref_fasta = tmp_path / "ref.fa"
    query_fasta = tmp_path / "query.fa"
    _write_fasta(ref_fasta, {"chr1": "A" * 1000, "chr2": "A" * 1000})
    _write_fasta(query_fasta, {"chrQ": "T" * 1000, "chrQ2": "T" * 1000, "chrQ3": "T" * 1000})
    return ref_fasta, query_fasta


# UCSC .2bit two-bit-per-base encoding, empirically verified by round-tripping
# through Bio.SeqIO.TwoBitIO.TwoBitIterator (Biopython's decoder is a compiled
# extension, so this was confirmed by brute-force trial rather than taken on
# faith from the spec): T=0, C=1, A=2, G=3, packed 4 bases/byte MSB-first.
_TWOBIT_BASE_CODES = {"T": 0, "C": 1, "A": 2, "G": 3}


def _write_twobit(path: Path, records: dict[str, str]) -> None:
    # Hand-build a minimal valid .2bit file: header, name/offset index, then
    # per-sequence dnaSize + nBlocks (N-runs) + maskBlocks (lowercase runs) +
    # packed sequence. No real 2bit-writing tool is used anywhere in this
    # pipeline or in Biopython (which only reads .2bit), so this is the only
    # way to get a small, fast, from-scratch test fixture.
    header = struct.pack("<IIII", 0x1A412743, 0, len(records), 0)
    index_entries = [(name.encode("ascii"), len(name.encode("ascii"))) for name in records]
    index_total_len = sum(1 + n + 4 for _, n in index_entries)
    cursor = len(header) + index_total_len

    def _runs(text: str, predicate) -> list[tuple[int, int]]:
        runs: list[tuple[int, int]] = []
        start = None
        for i, c in enumerate(text):
            if predicate(c) and start is None:
                start = i
            elif not predicate(c) and start is not None:
                runs.append((start, i))
                start = None
        if start is not None:
            runs.append((start, len(text)))
        return runs

    seq_blobs = []
    offsets = []
    for name, seq in records.items():
        offsets.append(cursor)
        upper_seq = seq.upper()
        n_blocks = _runs(upper_seq, lambda c: c == "N")
        mask_blocks = _runs(seq, str.islower)

        codes = [_TWOBIT_BASE_CODES.get(c, 0) for c in upper_seq]
        packed = bytearray()
        for i in range(0, len(codes), 4):
            chunk = codes[i : i + 4]
            chunk = chunk + [0] * (4 - len(chunk))
            packed.append((chunk[0] << 6) | (chunk[1] << 4) | (chunk[2] << 2) | chunk[3])

        rec = struct.pack("<I", len(seq))
        rec += struct.pack("<I", len(n_blocks))
        rec += b"".join(struct.pack("<I", s) for s, _ in n_blocks)
        rec += b"".join(struct.pack("<I", e - s) for s, e in n_blocks)
        rec += struct.pack("<I", len(mask_blocks))
        rec += b"".join(struct.pack("<I", s) for s, _ in mask_blocks)
        rec += b"".join(struct.pack("<I", e - s) for s, e in mask_blocks)
        rec += struct.pack("<I", 0)
        rec += bytes(packed)
        seq_blobs.append(rec)
        cursor += len(rec)

    index_bytes = b"".join(
        struct.pack("<B", n) + name_b + struct.pack("<I", offset)
        for (name_b, n), offset in zip(index_entries, offsets)
    )
    path.write_bytes(header + index_bytes + b"".join(seq_blobs))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _classify_with_filtering(
    tmp_path: Path,
    bed_text: str,
    chain_lines: list[str],
    ref_records: dict[str, str],
    query_records: dict[str, str],
    **classify_kwargs,
):
    # Real index -> candidates -> overlap stats -> classify_orthology chain,
    # with multi-projection filtering on. Returns (result, query orthology
    # rows, pruned-projection rows).
    from stola.orthology import classify_orthology

    tmp_path.mkdir(parents=True, exist_ok=True)
    bed = tmp_path / "cnees.bed"
    chain = tmp_path / "test.chain"
    index = tmp_path / "index.json"
    candidates = tmp_path / "candidates.tsv"
    stats = tmp_path / "stats.tsv"
    orthology = tmp_path / "orthology.tsv"
    query_bed = tmp_path / "query.bed"
    culled = tmp_path / "culled.tsv"
    ref_fasta = tmp_path / "ref.fa"
    query_fasta = tmp_path / "query.fa"

    bed.write_text(bed_text, encoding="utf-8")
    chain.write_text("\n".join(chain_lines), encoding="utf-8")
    _write_fasta(ref_fasta, ref_records)
    _write_fasta(query_fasta, query_records)
    build_index(chain, index)
    write_candidates(index, bed, candidates)
    comparator = SequenceComparator(ref_fasta, query_fasta, tmp_path / "work")
    try:
        write_overlap_stats(index, candidates, stats, sequence_comparator=comparator)
        result = classify_orthology(
            stats,
            bed,
            index,
            orthology,
            query_bed,
            sequence_comparator=comparator,
            filter_multi_projections=True,
            out_culled_projections=culled,
            **classify_kwargs,
        )
    finally:
        comparator.close()
    return result, _read_tsv(orthology), _read_tsv(culled)


class WorkflowTests(unittest.TestCase):
    def test_build_index_records_byte_offsets_and_extractable_plain_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, _ = _write_inputs(tmp_path)
            index = tmp_path / "index.json"

            built = build_index(chain, index)
            records = load_index(index).records

        self.assertEqual(built.records, records)
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].chain_id, "1")
        self.assertEqual(records[0].start_byte, 0)
        self.assertGreater(records[1].start_byte, records[0].start_byte)
        self.assertGreater(records[0].byte_length, 0)

    def test_build_index_auto_decompresses_gzipped_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, _ = _write_inputs(tmp_path)
            gz_chain = tmp_path / "test.chain.gz"
            gz_chain.unlink(missing_ok=True)
            with gzip.open(gz_chain, "wb") as handle:
                handle.write(chain.read_bytes())
            index = tmp_path / "gz-index.json"

            build_index(gz_chain, index)
            chain_index = load_index(index)
            records = chain_index.records

            self.assertEqual(len(records), 3)
            self.assertEqual(chain_index.indexed_chain_path.name, "test.chain")
            self.assertEqual(
                chain_index.indexed_chain_path.read_text(encoding="utf-8"),
                CHAIN_TEXT,
            )
            self.assertEqual(records[2].chain_id, "3")

    def test_candidate_chains_selects_span_overlaps_including_gapped_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            build_index(chain, index)

            count = write_candidates(index, bed, candidates)
            rows = _read_tsv(candidates)

        self.assertEqual(count, 3)
        self.assertEqual(
            [(row["cnee_id"], row["chain_id"]) for row in rows],
            [
                ("CNEE_PARTIAL", "2"),
                ("CNEE_PARTIAL", "1"),
                ("CNEE_GAPPED", "3"),
            ],
        )
        self.assertEqual(rows[2]["span_overlap_bases"], "50")

    def test_candidate_chains_require_exact_chromosome_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, _ = _write_inputs(tmp_path)
            bed = tmp_path / "versioned_cnees.bed"
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            # Chain uses "chr1"; BED uses the same chromosome but with an
            # accession version suffix. There is no alias fallback anymore,
            # so this must raise an informative error instead of silently
            # producing zero candidates for this chromosome.
            bed.write_text("chr1.1\t120\t240\tCNEE_VERSIONED\n", encoding="utf-8")
            build_index(chain, index)

            with self.assertRaises(ValueError) as ctx:
                write_candidates(index, bed, candidates)

        self.assertIn("chr1.1", str(ctx.exception))
        self.assertIn("chr1", str(ctx.exception))

    def test_candidate_chains_allow_unmatched_bed_chromosome_with_no_alias_hint(
        self,
    ) -> None:
        # A BED chromosome absent from the chain entirely (not just a
        # version-suffix mismatch) is a normal "no coverage" case, not an error.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, _ = _write_inputs(tmp_path)
            bed = tmp_path / "other_scaffold.bed"
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            bed.write_text("chrUn_scaffold9\t0\t10\tCNEE_UNPLACED\n", encoding="utf-8")
            build_index(chain, index)

            count = write_candidates(index, bed, candidates)
            rows = _read_tsv(candidates)

        self.assertEqual(count, 0)
        self.assertEqual(rows, [])

    def test_overlap_stats_compute_coverage_gaps_and_ranks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            build_index(chain, index)
            write_candidates(index, bed, candidates)

            rows = calculate_overlap_rows(index, candidates)
            by_pair = {(row["cnee_id"], row["chain_id"]): row for row in rows}

        partial_chain_1 = by_pair[("CNEE_PARTIAL", "1")]
        self.assertEqual(partial_chain_1["cnee_len"], 120)
        self.assertEqual(partial_chain_1["aligned_bases"], 70)
        self.assertEqual(partial_chain_1["gap_bases_within_chain_span"], 50)
        self.assertEqual(partial_chain_1["query_aligned_start"], 520)
        self.assertEqual(partial_chain_1["query_aligned_end"], 590)
        self.assertEqual(partial_chain_1["query_block_count"], 1)
        self.assertEqual(partial_chain_1["query_block_sizes"], "70,")
        self.assertEqual(partial_chain_1["query_block_starts"], "0,")
        self.assertEqual(partial_chain_1["num_candidate_chains_for_cnee"], 2)
        self.assertEqual(partial_chain_1["rank_by_chain_score"], 2)

        partial_chain_2 = by_pair[("CNEE_PARTIAL", "2")]
        self.assertEqual(partial_chain_2["aligned_bases"], 120)
        self.assertEqual(partial_chain_2["gap_bases_within_chain_span"], 0)
        self.assertEqual(partial_chain_2["query_aligned_start"], 770)
        self.assertEqual(partial_chain_2["query_aligned_end"], 890)
        self.assertEqual(partial_chain_2["rank_by_chain_score"], 1)

        gapped = by_pair[("CNEE_GAPPED", "3")]
        self.assertEqual(gapped["aligned_bases"], 0)
        self.assertEqual(gapped["gap_bases_within_chain_span"], 50)
        self.assertEqual(gapped["query_aligned_start"], ".")
        self.assertEqual(gapped["query_aligned_end"], ".")
        self.assertEqual(gapped["query_block_count"], 0)

    def test_end_to_end_writes_overlap_stats_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            count = write_overlap_stats(index, candidates, stats)

            rows = _read_tsv(stats)

        self.assertEqual(count, 3)
        self.assertEqual(rows[0]["cnee_id"], "CNEE_PARTIAL")
        self.assertEqual(rows[0]["chain_id"], "2")
        self.assertEqual(rows[0]["aligned_bases"], "120")
        self.assertEqual(rows[-1]["cnee_id"], "CNEE_GAPPED")

    def test_summary_from_overlap_stats_lists_aligned_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"
            summary = tmp_path / "summary.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            write_overlap_stats(index, candidates, stats)
            count = summarize_overlap_stats(stats, summary, cnee_bed=bed)
            rows = _read_tsv(summary)

        self.assertEqual(count, 3)
        self.assertEqual(list(rows[0]), SUMMARY_COLUMNS)
        self.assertEqual(rows[0]["cnee_uid"], "CNEE_PARTIAL")
        self.assertEqual(rows[0]["num_query_alignments"], "2")
        self.assertEqual(rows[0]["overlap_fractions"], "1.000,0.583")
        self.assertEqual(rows[0]["gap_bases"], "0,50")
        self.assertEqual(rows[1]["cnee_uid"], "CNEE_GAPPED")
        self.assertEqual(rows[1]["num_query_alignments"], "0")
        self.assertEqual(rows[1]["overlap_fractions"], "")
        self.assertEqual(rows[1]["gap_bases"], "")
        self.assertEqual(rows[2]["cnee_uid"], "CNEE_NO_CHAIN")
        self.assertEqual(rows[2]["num_query_alignments"], "0")

    def test_alignment_histogram_reads_summary_rows_past_csv_field_limit(self) -> None:
        from stola.plots import write_alignment_histogram

        n_alignments = 30_000  # ~180k chars per list, past csv's default 131,072
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            summary = tmp_path / "summary.tsv"
            values = ",".join(["0.583"] * n_alignments)
            summary.write_text(
                "\t".join(SUMMARY_COLUMNS) + "\n"
                + f"CNEE_REPEAT\t{n_alignments}\t{values}\t{values}\t{values}\t{values}\n"
                + "CNEE_SINGLE\t1\t1.000\t0\t0.950\t2\n",
                encoding="utf-8",
            )

            total = write_alignment_histogram(summary, tmp_path / "hist.svg")

        self.assertEqual(total, 2)

    def test_cli_rejects_out_of_range_options(self) -> None:
        from stola.cli import build_parser

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
            inputs = [str(chain), str(bed), str(ref_fasta), str(query_fasta)]
            parser = build_parser()

            args = parser.parse_args(
                inputs
                + [
                    "--soft-depth-limit", "1",
                    "--relative-score-threshold", "0",
                    "--flank-window", "1",
                    "--flank-any-coverage-threshold", "0.99",
                    "--min-gap-length", "1",
                ]
            )
            self.assertEqual(args.soft_depth_limit, 1)
            self.assertEqual(args.flank_any_coverage_threshold, 0.99)

            bad_options = [
                ["--soft-depth-limit", "0"],
                ["--relative-score-threshold", "1"],
                ["--relative-score-threshold", "-0.1"],
                ["--relative-score-threshold", "nan"],
                ["--flank-window", "0"],
                ["--flank-any-coverage-threshold", "1.0"],
                ["--min-gap-length", "0"],
                ["--soft-depth-limit", "2.5"],
            ]
            for option in bad_options:
                with self.subTest(option=option), contextlib.redirect_stderr(io.StringIO()) as err:
                    with self.assertRaises(SystemExit):
                        parser.parse_args(inputs + option)
                    self.assertIn(option[0], err.getvalue())

            with contextlib.redirect_stderr(io.StringIO()) as err:
                with self.assertRaises(SystemExit):
                    parser.parse_args([str(tmp_path / "missing.chain")] + inputs[1:])
            self.assertIn("file not found", err.getvalue())

    def test_cli_reports_version(self) -> None:
        from stola import __version__
        from stola.cli import build_parser

        with contextlib.redirect_stdout(io.StringIO()) as out:
            with self.assertRaises(SystemExit):
                build_parser().parse_args(["--version"])
        self.assertEqual(out.getvalue().strip(), f"stola {__version__}")

    def test_cli_passes_every_option_to_the_workflow(self) -> None:
        from unittest import mock

        from stola import cli

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
            argv = [
                str(chain), str(bed), str(ref_fasta), str(query_fasta),
                "--outdir", str(tmp_path / "out"),
                "--query-prefix", "aphWoo",
                "--no-filter-multi-projections",
                "--soft-depth-limit", "7",
                "--relative-score-threshold", "0.3",
                "--flank-window", "500",
                "--flank-any-coverage-threshold", "0.25",
                "--min-gap-length", "20",
                "--keep-intermediates",
            ]
            with mock.patch.object(cli, "run_workflow") as run_workflow:
                with contextlib.redirect_stdout(io.StringIO()):
                    cli.main(argv)

        args, kwargs = run_workflow.call_args
        self.assertEqual(args, (chain, bed, ref_fasta, query_fasta))
        self.assertEqual(
            kwargs,
            {
                "outdir": tmp_path / "out",
                "query_id_prefix": "aphWoo",
                "filter_multi_projections": False,
                "soft_depth_limit": 7,
                "relative_score_threshold": 0.3,
                "flank_window": 500,
                "flank_any_coverage_threshold": 0.25,
                "min_gap_length": 20,
                "keep_intermediates": True,
                "progress": cli._print_progress,
            },
        )

    def test_workflow_reports_each_stage_as_it_finishes(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
            # Record, for each message, whether the stats file and the ref
            # orthology table existed when it was reported.
            stats_path = tmp_path / "overlap_stats.tsv"
            ref_path = tmp_path / "ref_ortho_status.tsv"
            messages: list[tuple[str, bool, bool]] = []
            run_workflow(
                chain,
                bed,
                ref_fasta,
                query_fasta,
                outdir=tmp_path,
                progress=lambda message: messages.append(
                    (message, stats_path.exists(), ref_path.exists())
                ),
            )

        texts = [message for message, _, _ in messages]
        self.assertTrue(texts[0].startswith("Indexed 3 chains"), texts)
        self.assertEqual(
            [text.split(" -> ")[0] for text in texts if text.startswith("Wrote")][:2],
            ["Wrote 3 candidate chain/CNEE pairs", "Wrote 3 all overlap-stat rows"],
        )
        # Messages arrive between stages, not all at the end: the indexing
        # message comes before overlap stats are written, and the overlap
        # message before classification writes the ref orthology table.
        self.assertEqual(messages[0][1:], (False, False))
        overlap = next(m for m in messages if m[0].startswith("Wrote 3 all overlap-stat rows"))
        self.assertEqual(overlap[1:], (True, False))
        self.assertTrue(messages[-2][0].startswith("Wrote final category counts"))
        self.assertTrue(messages[-1][0].startswith("Removed intermediates"))

    def test_workflow_removes_intermediates_after_success(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        for keep in (False, True):
            with self.subTest(keep_intermediates=keep), tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                chain, bed = _write_inputs(tmp_path)
                gz_chain = tmp_path / "test.chain.gz"
                with gzip.open(gz_chain, "wb") as handle:
                    handle.write(chain.read_bytes())
                chain.unlink()
                ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
                outdir = tmp_path / "out"

                result = run_workflow(
                    gz_chain, bed, ref_fasta, query_fasta, outdir=outdir, keep_intermediates=keep
                )

                work_files = sorted(path.name for path in (outdir / "work").iterdir())
                intermediates = ["candidate_chains.tsv", "chain_index.json", "test.chain"]
                genome_indexes = ["query_fasta.idx", "ref_fasta.idx"]
                expected = sorted(intermediates + genome_indexes) if keep else genome_indexes
                self.assertEqual(work_files, expected)
                self.assertTrue(gz_chain.exists())
                # Outputs are unaffected.
                self.assertEqual(len(_read_tsv(result.ref_cnees_path)), 3)
                self.assertTrue(result.query_bed_path.exists())

    def test_workflow_never_deletes_a_plain_chain_input(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            # The plain chain sits inside the run's own work/ folder, where a
            # decompressed copy would go, and must still survive cleanup.
            work_dir = tmp_path / "work"
            work_dir.mkdir()
            chain, bed = _write_inputs(tmp_path)
            plain_chain = chain.rename(work_dir / chain.name)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)

            run_workflow(plain_chain, bed, ref_fasta, query_fasta, outdir=tmp_path)

            self.assertEqual(plain_chain.read_text(encoding="utf-8"), CHAIN_TEXT)
            self.assertFalse((work_dir / "chain_index.json").exists())
            self.assertFalse((work_dir / "candidate_chains.tsv").exists())

    def test_classification_reports_filtering_progress(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # One ref with a candidate on each of two query chromosomes: an
        # initial one2many, so its two projections are filtered at two loci.
        messages: list[str] = []
        with tempfile.TemporaryDirectory() as tmpdir:
            _classify_with_filtering(
                Path(tmpdir),
                "chr1\t100\t150\tMULTI\n",
                [
                    "chain 100 chr1 300 + 50 250 chrQ1 300 + 50 250 1",
                    "200",
                    "",
                    "chain 90 chr1 300 + 50 250 chrQ2 300 + 50 250 2",
                    "200",
                    "",
                ],
                {"chr1": "A" * 300},
                {"chrQ1": "A" * 300, "chrQ2": "A" * 300},
                flank_window=50,
                progress=messages.append,
            )

        self.assertEqual(
            messages,
            [
                "Read 2 aligned projections of 1 reference CNEEs",
                "Filtering 2 projections of 1 multi-projection reference CNEEs, at 2 query loci ...",
                # With only two loci, the 10%..50% milestones all land on the
                # first and collapse into one message.
                "Scored 1/2 query loci (50%)",
                "Scored 2/2 query loci (100%)",
                "Assigning final orthology categories ...",
            ],
        )

    def test_cli_prints_timestamped_progress(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.cli import main

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
            outdir = tmp_path / "out"
            with contextlib.redirect_stdout(io.StringIO()) as out:
                main([str(chain), str(bed), str(ref_fasta), str(query_fasta), "--outdir", str(outdir)])

        lines = out.getvalue().splitlines()
        for line in lines:
            self.assertRegex(line, r"^\[\d\d:\d\d:\d\d\] ")
        self.assertIn("writing outputs under", lines[0])
        self.assertIn("Indexed 3 chains", lines[1])
        self.assertRegex(lines[-1], r"All outputs and intermediates written under .* \(\d+m \d\ds\)$")

    def test_cli_runs_workflow_from_bare_arguments(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.cli import main

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)
            outdir = tmp_path / "out"

            main([str(chain), str(bed), str(ref_fasta), str(query_fasta), "--outdir", str(outdir)])

            ref_rows = _read_tsv(outdir / "ref_ortho_status.tsv")
            plots = {
                name: (outdir / "summary_figures" / name).read_text(encoding="utf-8")
                for name in (
                    "num_query_alignments.svg",
                    "overlap_fraction_density.svg",
                    "sequence_identity_density.svg",
                )
            }

        self.assertEqual(len(ref_rows), 3)
        for name, svg in plots.items():
            self.assertTrue(svg.startswith("<svg"), name)
            self.assertIn("</svg>", svg, name)
        # The fixture's aligned candidate rows give the density plots a curve.
        self.assertIn("<polyline", plots["overlap_fraction_density.svg"])

    def test_single_workflow_writes_stats_and_query_bed(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)

            result = run_workflow(chain, bed, ref_fasta, query_fasta, outdir=tmp_path)
            stats_rows = _read_tsv(result.stats_path)
            bed_lines = result.query_bed_path.read_text(encoding="utf-8").splitlines()
            summary_rows = _read_tsv(result.summary_path)
            orthology_rows = _read_tsv(result.orthology_path)
            summary_figures =sorted(path.name for path in (tmp_path / "summary_figures").iterdir())

        self.assertEqual(result.summary_path, tmp_path / "summary_figures" / "overlap_stats.summary.tsv")
        self.assertEqual(
            summary_figures,
            [
                "final_category_counts.svg",
                "final_category_counts.tsv",
                "initial_to_final_reassignment.svg",
                "initial_to_final_reassignment.tsv",
                "num_query_alignments.svg",
                "overlap_fraction_density.svg",
                "overlap_stats.summary.tsv",
                "rejection_reason_counts.svg",
                "rejection_reason_counts.tsv",
                "sequence_identity_density.svg",
            ],
        )
        self.assertEqual(result.summary_plots.histogram_n, 3)
        summary = result.orthology_summary
        # No query FASTA gap convention can be inferred from this tiny fixture,
        # so the lost ref collapses to missing/loss.
        self.assertEqual(
            summary.final_category_counts,
            {"uncertain_loss": 1, "missing/loss": 1, "unaligned": 1},
        )
        self.assertEqual(summary.reassignment_counts["one2many"], {"uncertain_loss": 1})
        # Both CNEE_PARTIAL projections fail the single-chain flank check;
        # every reason is listed even at zero.
        self.assertEqual(
            summary.rejection_reason_counts, {"depth": 0, "no_flank_support": 2, "low_score": 0}
        )
        self.assertEqual(result.chain_count, 3)
        self.assertEqual(result.candidate_count, 3)
        self.assertEqual(result.stats_count, 3)
        # CNEE_PARTIAL spans two candidate chains (chrQ, chrQ2) -- an initial
        # one2many shape needing filtering. Each chain is the only chain on
        # its own query chromosome in this fixture (no genome-wide
        # competitor at all), so each candidate's own locus has depth==1 --
        # the "isolated" branch, accepted only if its own chain's flanks have
        # real coverage. Neither chain extends meaningfully past its own
        # aligned span here, so both fail that check and get rejected; this
        # ref ends up uncertain_loss with zero query CNEE loci.
        self.assertEqual(result.query_bed_count, 0)
        self.assertEqual(result.summary_count, 3)
        self.assertEqual(result.orthology_result.orthology_rows, 3)
        self.assertEqual(result.orthology_result.culled_projection_rows, 2)
        self.assertEqual(len(stats_rows), 3)
        self.assertEqual(len(summary_rows), 3)
        self.assertEqual(len(orthology_rows), 3)
        self.assertEqual(
            list(orthology_rows[0]),
            ["query_cnee_id", "partner_ref_cnees", "orthology_category", "sequence_identity", "gap_fraction"],
        )
        self.assertEqual(bed_lines, [])
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in orthology_rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_PARTIAL"], "uncertain_loss")
        self.assertEqual(na_categories["CNEE_GAPPED"], "lost")
        self.assertEqual(na_categories["CNEE_NO_CHAIN"], "unaligned")
        summary_alignments = {
            row["cnee_uid"]: row["num_query_alignments"] for row in summary_rows
        }
        self.assertEqual(
            summary_alignments, {"CNEE_PARTIAL": "2", "CNEE_GAPPED": "0", "CNEE_NO_CHAIN": "0"}
        )

    def test_workflow_leaves_named_ref_bed_path_none_when_input_bed_has_names(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)

            result = run_workflow(chain, bed, ref_fasta, query_fasta, outdir=tmp_path)

            self.assertIsNone(result.named_ref_bed_path)
            self.assertEqual(list(tmp_path.glob("*_named.bed")), [])

    def test_workflow_writes_named_ref_bed_when_input_bed_lacks_name_column(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, _ = _write_inputs(tmp_path)
            bed = tmp_path / "cnees.bed"
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t120\t240",
                        "chr2\t175\t225",
                        "chr3\t10\t20",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            ref_fasta, query_fasta = _write_fasta_inputs(tmp_path)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = run_workflow(chain, bed, ref_fasta, query_fasta, outdir=tmp_path)

            self.assertEqual(result.named_ref_bed_path, tmp_path / "cnees_named.bed")
            self.assertEqual(
                result.named_ref_bed_path.read_text(encoding="utf-8").splitlines(),
                [
                    "chr1\t120\t240\tcnee_000001_chr1_120_240",
                    "chr2\t175\t225\tcnee_000002_chr2_175_225",
                    "chr3\t10\t20\tcnee_000003_chr3_10_20",
                ],
            )

    def test_orthology_graph_edges_store_per_projection_metrics(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import (
            _build_orthology_graph,
            _cluster_query_projections,
            _read_projection_evidence,
            _ref_node,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            stats = tmp_path / "stats.tsv"
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                # Good identity but bad gap on its own.
                + "CNEE_A\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\t1\n"
                # Good gap but bad identity on its own.
                + "CNEE_A\tchr1\t0\t50\t1\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\t2\n"
                + "CNEE_B\tchr1\t100\t150\t2\t40\tchrQ\t+\t300\t340\tB#1\t0.25\t.\t1\n",
                encoding="utf-8",
            )

            evidence = _read_projection_evidence(stats)
            clusters = _cluster_query_projections(evidence.projections)
            graph = _build_orthology_graph(["CNEE_A", "CNEE_B"], clusters)

        edge_a = graph.edges[_ref_node("CNEE_A"), "query:000001"]
        self.assertEqual(edge_a["projection_metrics"], [(0.9, 0.5), (0.3, 0.05)])
        self.assertEqual(edge_a["projection_count"], 2)

        edge_b = graph.edges[_ref_node("CNEE_B"), "query:000002"]
        self.assertEqual(edge_b["projection_metrics"], [(None, 0.25)])
        self.assertEqual(edge_b["projection_count"], 1)

    def test_orthology_classification_reports_every_lone_edge_as_one2one(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t0\t50\tCNEE_HIGH_ID",
                        "chr1\t100\t150\tCNEE_SEQ_MISSING",
                        "chr1\t200\t250\tCNEE_LOW_ID",
                        "chr1\t300\t350\tCNEE_GAP_ONLY",
                        "chr1\t400\t450\tCNEE_MISSING",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                + "CNEE_HIGH_ID\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tS#1\t0.05\t0.8\t1\n"
                + "CNEE_SEQ_MISSING\tchr1\t100\t150\t2\t50\tchrQ\t+\t300\t350\tM#1\t0.05\t.\t1\n"
                + "CNEE_LOW_ID\tchr1\t200\t250\t3\t50\tchrQ\t+\t500\t550\tN#1\t0.2\t0.5\t1\n"
                + "CNEE_GAP_ONLY\tchr1\t300\t350\t4\t0\tchrQ\t+\t.\t.\tG#1\t1.0\t.\t1\n",
                encoding="utf-8",
            )

            # No multi-projection refs here, so filtering is off and the
            # index path is never touched.
            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                filter_multi_projections=False,
            )
            rows = _read_tsv(orthology)
            bed_lines = query_bed.read_text(encoding="utf-8").splitlines()

        self.assertEqual(result.orthology_rows, 5)
        # A lone ref/query edge is never filtered, so every aligned ref gets
        # a query locus, whatever its identity or gap fraction.
        self.assertEqual(result.query_bed_rows, 3)
        self.assertEqual(
            bed_lines,
            [
                "chrQ\t100\t150\tq_CNEE_HIGH_ID\t0\t+",
                "chrQ\t300\t350\tq_CNEE_SEQ_MISSING\t0\t+",
                "chrQ\t500\t550\tq_CNEE_LOW_ID\t0\t+",
            ],
        )
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in rows
        }
        self.assertEqual(categories_by_query["q_CNEE_HIGH_ID"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_SEQ_MISSING"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_LOW_ID"], "one2one")
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_GAP_ONLY"], "lost")
        self.assertEqual(na_categories["CNEE_MISSING"], "unaligned")

    def test_classify_orthology_uses_custom_query_id_prefix(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            bed.write_text("chr1\t0\t50\tCNEE_A\n", encoding="utf-8")
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                + "CNEE_A\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.05\t0.8\t1\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                filter_multi_projections=False,
                query_id_prefix="cust_",
            )
            rows = _read_tsv(orthology)
            bed_lines = query_bed.read_text(encoding="utf-8").splitlines()

        self.assertEqual(result.query_bed_rows, 1)
        query_ids = {row["query_cnee_id"] for row in rows}
        self.assertEqual(query_ids, {"cust_CNEE_A"})
        self.assertEqual(bed_lines, ["chrQ\t100\t150\tcust_CNEE_A\t0\t+"])

    def test_classify_orthology_appends_missing_underscore_to_query_id_prefix(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            bed.write_text("chr1\t0\t50\tCNEE_A\n", encoding="utf-8")
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                + "CNEE_A\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.05\t0.8\t1\n",
                encoding="utf-8",
            )

            classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                filter_multi_projections=False,
                query_id_prefix="aphWoo",
            )
            rows = _read_tsv(orthology)

        query_ids = {row["query_cnee_id"] for row in rows}
        self.assertEqual(query_ids, {"aphWoo_CNEE_A"})

    def test_relative_score_filtering_resolves_one2many_to_one2one_by_rejecting_worse_candidate(
        self,
    ) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # A single ref CNEE, MULTI = chr1:100-150, with two candidate
        # chains -- an initial one2many shape:
        #   chain 1: chr1[50,250) <-> chrQ1[50,250), a perfect match
        #     throughout. No other chain ever touches chrQ1, so its locus is
        #     isolated; its own chain extends 50bp past the CNEE core on both
        #     sides, giving full flank coverage from itself -- accepted.
        #   chain 2: chr1[80,220) <-> chrQ2[80,220), a total mismatch
        #     throughout (chr1 is all "A", chrQ2 is all "T"). chain 3
        #     (chr2[500,650) <-> chrQ2[100,250), all "T" on both sides -- a
        #     perfect match) is NOT a candidate for any CNEE (chr2 has no ref
        #     CNEE at all) but it aligns over chain 2's locus on chrQ2, so it
        #     competes there. Over the locus slice chrQ2[50,200), chain 3's
        #     local score is identity 1.0 x aligned fraction 100/150, chain
        #     2's is identity 0 x anything = 0 -- relative score 0, rejected.
        # After filtering, MULTI should have exactly one surviving
        # projection (chain 1) and resolve from one2many down to one2one.
        with tempfile.TemporaryDirectory() as tmpdir:
            result, orthology_rows, culled_rows = _classify_with_filtering(
                Path(tmpdir),
                "chr1\t100\t150\tMULTI\n",
                [
                    "chain 100 chr1 300 + 50 250 chrQ1 300 + 50 250 1",
                    "200",
                    "",
                    "chain 100 chr1 300 + 80 220 chrQ2 300 + 80 220 2",
                    "140",
                    "",
                    "chain 100 chr2 700 + 500 650 chrQ2 300 + 100 250 3",
                    "150",
                    "",
                ],
                {"chr1": "A" * 300, "chr2": "T" * 700},
                {"chrQ1": "A" * 300, "chrQ2": "T" * 300},
                flank_window=50,
            )

        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in orthology_rows
        }
        self.assertEqual(categories_by_query["q_MULTI"], "one2one")
        self.assertEqual(result.query_bed_rows, 1)
        self.assertEqual(result.culled_projection_rows, 1)

        self.assertEqual(len(culled_rows), 1)
        culled_row = culled_rows[0]
        self.assertEqual(culled_row["ref_cnee_id"], "MULTI")
        self.assertEqual(culled_row["query_chrom"], "chrQ2")
        self.assertEqual(culled_row["initial_orthology_category"], "one2many")
        self.assertEqual(culled_row["rejection_reason"], "low_score")
        self.assertEqual(culled_row["locus_depth"], "2")
        self.assertEqual(float(culled_row["local_score"]), 0.0)
        self.assertAlmostEqual(float(culled_row["best_competitor_local_score"]), 100 / 150)
        self.assertEqual(float(culled_row["relative_score"]), 0.0)
        self.assertEqual(culled_row["flank_any_coverage"], ".")
        self.assertEqual(float(culled_row["sequence_identity"]), 0.0)

    def test_relative_score_filtering_keeps_many2one_when_chains_score_alike(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # Two ref CNEEs on different ref chromosomes, each carried by its own
        # chain onto the SAME query locus chrQ[100,150) with identical
        # alignment (perfect match over chrQ[50,250)). Both chains score the
        # same at that locus (relative score 1.0 each), so both projections
        # survive and the locus stays a many2one call -- a single winner per
        # locus is not forced.
        with tempfile.TemporaryDirectory() as tmpdir:
            result, orthology_rows, culled_rows = _classify_with_filtering(
                Path(tmpdir),
                "chr1\t100\t150\tREF_A\nchr2\t100\t150\tREF_B\n",
                [
                    "chain 100 chr1 300 + 50 250 chrQ 300 + 50 250 1",
                    "200",
                    "",
                    "chain 100 chr2 300 + 50 250 chrQ 300 + 50 250 2",
                    "200",
                    "",
                ],
                {"chr1": "A" * 300, "chr2": "A" * 300},
                {"chrQ": "A" * 300},
                flank_window=50,
            )

        self.assertEqual(result.culled_projection_rows, 0)
        self.assertEqual(culled_rows, [])
        self.assertEqual(len(orthology_rows), 1)
        self.assertEqual(orthology_rows[0]["orthology_category"], "many2one")
        self.assertEqual(sorted(orthology_rows[0]["partner_ref_cnees"].split(",")), ["REF_A", "REF_B"])

    def test_soft_depth_limit_rejects_chain_outside_top_chain_scores(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # MULTI = chr1:100-150 has two candidate chains: chain 1 onto an
        # isolated, well-supported chrQ1 locus, and chain 2 (chain score 100)
        # onto chrQ2[100,150), where chain 3 (chain score 500, not a
        # candidate for any CNEE) also aligns. Everything is a perfect
        # match, and chain 2 covers more of the chrQ2 slice than chain 3, so
        # with the default soft_depth_limit chain 2 is the best chain there
        # and MULTI stays one2many. With soft_depth_limit=1 only chain 3 (the
        # higher chain score) is kept at chrQ2, so chain 2's projection is
        # rejected without being scored and MULTI resolves to one2one.
        inputs = (
            "chr1\t100\t150\tMULTI\n",
            [
                "chain 100 chr1 300 + 50 250 chrQ1 300 + 50 250 1",
                "200",
                "",
                "chain 100 chr1 300 + 80 220 chrQ2 300 + 80 220 2",
                "140",
                "",
                "chain 500 chr2 700 + 500 650 chrQ2 300 + 100 250 3",
                "150",
                "",
            ],
            {"chr1": "A" * 300, "chr2": "A" * 700},
            {"chrQ1": "A" * 300, "chrQ2": "A" * 300},
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            _, default_rows, default_culled = _classify_with_filtering(
                Path(tmpdir) / "default", *inputs, flank_window=50
            )
            _, limited_rows, limited_culled = _classify_with_filtering(
                Path(tmpdir) / "limited", *inputs, flank_window=50, soft_depth_limit=1
            )

        self.assertEqual({row["orthology_category"] for row in default_rows}, {"one2many"})
        self.assertEqual(default_culled, [])

        self.assertEqual([row["orthology_category"] for row in limited_rows], ["one2one"])
        self.assertEqual(len(limited_culled), 1)
        culled_row = limited_culled[0]
        self.assertEqual(culled_row["query_chrom"], "chrQ2")
        self.assertEqual(culled_row["rejection_reason"], "depth")
        self.assertEqual(culled_row["locus_depth"], "2")
        self.assertEqual(culled_row["local_score"], ".")

    def test_relative_score_filtering_scores_minus_strand_chain_on_reverse_complement(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # MULTI = chr1:150-200 has two candidate chains (initial one2many):
        #   chain 1: chr1 -> chrQ1, a perfect, isolated, self-supported
        #     match -- accepted.
        #   chain 2: chr1 -> chrQ2's MINUS strand, landing on chrQ2[100,150).
        #     chr1 is built so that, once chain 2's query bases are
        #     reverse-complemented, exactly 40% of that locus matches.
        #   chain 3: chr2 -> chrQ2's plus strand, a perfect match over the
        #     same locus (not a candidate: chr2 has no CNEE).
        # Chain 2's relative score is therefore exactly 0.4 and it is
        # rejected. Wrong minus-strand coordinates or a missing reverse
        # complement would give a different identity.
        rng = random.Random(7)
        query = "".join(rng.choice("ACGT") for _ in range(300))
        complement = str.maketrans("ACGT", "TGCA")
        chr1 = ["A"] * 300
        for offset in range(200):
            # chr1[50 + offset] aligns to minus-strand chrQ2, i.e. the
            # complement of plus-strand position 249 - offset.
            query_pos = 249 - offset
            base = query[query_pos]
            chr1[50 + offset] = base.translate(complement) if query_pos % 5 < 2 else base
        with tempfile.TemporaryDirectory() as tmpdir:
            result, orthology_rows, culled_rows = _classify_with_filtering(
                Path(tmpdir),
                "chr1\t150\t200\tMULTI\n",
                [
                    "chain 100 chr1 300 + 50 250 chrQ1 300 + 50 250 1",
                    "200",
                    "",
                    "chain 100 chr1 300 + 50 250 chrQ2 300 - 50 250 2",
                    "200",
                    "",
                    "chain 100 chr2 300 + 50 250 chrQ2 300 + 50 250 3",
                    "200",
                    "",
                ],
                {"chr1": "".join(chr1), "chr2": query},
                {"chrQ1": "".join(chr1), "chrQ2": query},
                flank_window=50,
            )

        self.assertEqual([row["orthology_category"] for row in orthology_rows], ["one2one"])
        self.assertEqual(result.culled_projection_rows, 1)
        culled_row = culled_rows[0]
        self.assertEqual(culled_row["query_chrom"], "chrQ2")
        self.assertEqual((culled_row["query_start"], culled_row["query_end"]), ("100", "150"))
        self.assertEqual(culled_row["rejection_reason"], "low_score")
        # Both the overlap stage and the filter compare against the reverse complement.
        self.assertAlmostEqual(float(culled_row["sequence_identity"]), 0.4)
        self.assertAlmostEqual(float(culled_row["local_score"]), 0.4)
        self.assertAlmostEqual(float(culled_row["best_competitor_local_score"]), 1.0)
        self.assertAlmostEqual(float(culled_row["relative_score"]), 0.4)

    def test_single_chain_locus_accepts_flank_support_from_another_chain(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")

        # MULTI = chr1:100-150 has two candidate chains (initial one2many).
        # Chain 1 lands on an isolated chrQ1 locus whose flanks it covers
        # itself, so it is accepted. Chain 2 covers only the CNEE,
        # chr1[100,150) -> chrQ2[100,150): it is the only chain with aligned
        # bases at that locus, and it leaves its own flanks unaligned.
        # Chain 3 (not a candidate -- chr2 has no CNEE) aligns both flanks,
        # chrQ2[50,100) and [150,200), but not the locus itself. With chain 3
        # present the flanks are covered by another chain, so chain 2's
        # projection is accepted and MULTI stays one2many; without it, chain
        # 2 is rejected as no_flank_support and MULTI resolves to one2one.
        chains = [
            "chain 100 chr1 300 + 50 250 chrQ1 300 + 50 250 1",
            "200",
            "",
            "chain 100 chr1 300 + 100 150 chrQ2 300 + 100 150 2",
            "50",
            "",
        ]
        flank_chain = ["chain 100 chr2 700 + 500 650 chrQ2 300 + 50 200 3", "50 50 50", "50", ""]
        bed_text = "chr1\t100\t150\tMULTI\n"
        genomes = ({"chr1": "A" * 300, "chr2": "A" * 700}, {"chrQ1": "A" * 300, "chrQ2": "A" * 300})
        with tempfile.TemporaryDirectory() as tmpdir:
            _, supported_rows, supported_culled = _classify_with_filtering(
                Path(tmpdir) / "supported", bed_text, chains + flank_chain, *genomes, flank_window=50
            )
            _, unsupported_rows, unsupported_culled = _classify_with_filtering(
                Path(tmpdir) / "unsupported", bed_text, chains, *genomes, flank_window=50
            )

        self.assertEqual({row["orthology_category"] for row in supported_rows}, {"one2many"})
        self.assertEqual(supported_culled, [])

        self.assertEqual([row["orthology_category"] for row in unsupported_rows], ["one2one"])
        self.assertEqual(len(unsupported_culled), 1)
        culled_row = unsupported_culled[0]
        self.assertEqual(culled_row["query_chrom"], "chrQ2")
        self.assertEqual(culled_row["rejection_reason"], "no_flank_support")
        self.assertEqual(float(culled_row["flank_any_coverage"]), 0.0)

    def test_ref_cnee_tsv_splits_lost_by_query_assembly_gap_overlap(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            ref_cnees = tmp_path / "ref_cnees.tsv"
            chain = tmp_path / "test.chain"
            index = tmp_path / "index.json"
            query_fasta = tmp_path / "query.fa"

            bed.write_text(
                "\n".join(
                    [
                        "chr1\t100\t150\tCNEE_ONE2ONE",
                        "chr1\t300\t350\tCNEE_LOST_GAP",
                        "chr1\t500\t550\tCNEE_LOST_REAL",
                        "chr2\t50\t100\tCNEE_WEAK",
                        "chr3\t10\t60\tCNEE_NO_CHAIN",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            # One chain, 3 aligned blocks: A (ref[100,150) query[1000,1050)),
            # B (ref[400,450) query[1200,1250)), C (ref[700,750) query[1400,1450)).
            # CNEE_LOST_GAP sits in the A/B ref gap, whose query-side gap
            # [1050,1200) will overlap a 50bp N-run in the query FASTA.
            # CNEE_LOST_REAL sits in the B/C ref gap, whose query-side gap
            # [1250,1400) will NOT overlap any N-run.
            chain.write_text(
                "\n".join(
                    [
                        "chain 100 chr1 1000 + 100 750 chrQ 3000 + 1000 1450 1",
                        "50\t250\t150",
                        "50\t250\t150",
                        "50",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            from stola.index import build_index

            build_index(chain, index)

            columns = [
                "cnee_uid",
                "ref_chrom",
                "ref_start",
                "ref_end",
                "chain_id",
                "aligned_bases",
                "query_chrom",
                "query_strand",
                "query_aligned_start",
                "query_aligned_end",
                "query_projection_name",
                "gap_fraction_within_chain_span",
                "seq_identity",
                "rank_by_chain_score",
            ]
            rows = [
                # One clean, high-identity one2one.
                ("CNEE_ONE2ONE", "chr1", 100, 150, "1", 50, "chrQ", "+", "1000", "1050", "ONE2ONE#1", "0.05", "0.8"),
                # Two "lost" candidates (aligned_bases=0) on the same chain.
                ("CNEE_LOST_GAP", "chr1", 300, 350, "1", 0, "chrQ", "+", ".", ".", "LOSTGAP#1", "1.0", "."),
                ("CNEE_LOST_REAL", "chr1", 500, 550, "1", 0, "chrQ", "+", ".", ".", "LOSTREAL#1", "1.0", "."),
                # A weak (low identity, high gap) lone ref/query edge -- with
                # only one candidate at all, it is never subject to
                # filtering and is always kept as plain one2one, with its
                # identity/gap exposed as output columns instead.
                ("CNEE_WEAK", "chr2", 50, 100, "2", 50, "chrQ", "+", "2000", "2050", "WEAK#1", "0.5", "0.9"),
            ]
            stats.write_text(
                "\t".join(columns)
                + "\n"
                + "\n".join(
                    "\t".join([str(field) for field in row] + ["1"]) for row in rows
                )
                + "\n",
                encoding="utf-8",
            )

            # Query FASTA: the deliberate 50bp gap-marker N-run at [1100,1150)
            # (inside CNEE_LOST_GAP's query gap), plus 9 more elsewhere so the
            # convention (>= MIN_RUNS_FOR_INFERENCE) is confidently inferred
            # as "50bp is this assembly's gap marker".
            seq = list("A" * 3000)
            seq[1100:1150] = ["N"] * 50
            padding_start = 1500
            for i in range(9):
                start = padding_start + i * 100
                seq[start : start + 50] = ["N"] * 50
            query_fasta.write_text(
                ">chrQ\n" + "\n".join("".join(seq[i : i + 60]) for i in range(0, len(seq), 60)) + "\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                index,
                orthology,
                query_bed,
                filter_multi_projections=True,
                query_fasta=query_fasta,
                out_ref_cnees=ref_cnees,
            )
            ref_rows = {row["ref_cnee_id"]: row for row in _read_tsv(ref_cnees)}

        self.assertEqual(result.ref_cnee_rows, 5)
        self.assertEqual(
            ref_rows["CNEE_ONE2ONE"],
            {
                "ref_cnee_id": "CNEE_ONE2ONE",
                "orthology_category": "one2one",
                "num_projections_final": "1",
                "num_projections_pruned": "0",
                "pre_pruning_category": "one2one",
            },
        )
        self.assertEqual(ref_rows["CNEE_LOST_GAP"]["orthology_category"], "assembly_gap_loss")
        self.assertEqual(ref_rows["CNEE_LOST_REAL"]["orthology_category"], "probable_loss")
        for ref_id in ("CNEE_LOST_GAP", "CNEE_LOST_REAL"):
            self.assertEqual(ref_rows[ref_id]["num_projections_final"], "0")
            self.assertEqual(ref_rows[ref_id]["pre_pruning_category"], ref_rows[ref_id]["orthology_category"])
        self.assertEqual(ref_rows["CNEE_WEAK"]["orthology_category"], "one2one")
        self.assertEqual(ref_rows["CNEE_WEAK"]["num_projections_final"], "1")
        self.assertEqual(ref_rows["CNEE_WEAK"]["num_projections_pruned"], "0")
        self.assertEqual(ref_rows["CNEE_WEAK"]["pre_pruning_category"], "one2one")
        self.assertEqual(ref_rows["CNEE_NO_CHAIN"]["orthology_category"], "unaligned")
        self.assertEqual(ref_rows["CNEE_NO_CHAIN"]["pre_pruning_category"], "unaligned")

    def test_ref_cnee_tsv_collapses_lost_categories_without_query_fasta(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            ref_cnees = tmp_path / "ref_cnees.tsv"
            chain = tmp_path / "test.chain"
            index = tmp_path / "index.json"

            bed.write_text("chr1\t300\t350\tCNEE_LOST_GAP\n", encoding="utf-8")
            chain.write_text(
                "\n".join(
                    [
                        "chain 100 chr1 1000 + 100 750 chrQ 3000 + 1000 1450 1",
                        "50\t250\t150",
                        "50\t250\t150",
                        "50",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            from stola.index import build_index

            build_index(chain, index)

            columns = [
                "cnee_uid",
                "ref_chrom",
                "ref_start",
                "ref_end",
                "chain_id",
                "aligned_bases",
                "query_chrom",
                "query_strand",
                "query_aligned_start",
                "query_aligned_end",
                "query_projection_name",
                "gap_fraction_within_chain_span",
                "seq_identity",
                "rank_by_chain_score",
            ]
            stats.write_text(
                "\t".join(columns)
                + "\n"
                + "CNEE_LOST_GAP\tchr1\t300\t350\t1\t0\tchrQ\t+\t.\t.\tLOSTGAP#1\t1.0\t.\t1\n",
                encoding="utf-8",
            )

            # No query_fasta at all -- can't even attempt the boundary check.
            result = classify_orthology(
                stats,
                bed,
                index,
                orthology,
                query_bed,
                filter_multi_projections=True,
                out_ref_cnees=ref_cnees,
            )
            ref_rows = {row["ref_cnee_id"]: row for row in _read_tsv(ref_cnees)}

        self.assertEqual(result.ref_cnee_rows, 1)
        self.assertEqual(ref_rows["CNEE_LOST_GAP"]["orthology_category"], "missing/loss")

    def test_orthology_clustering_keeps_opposite_strand_overlaps_separate(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            bed.write_text(
                "chr1\t0\t50\tCNEE_PLUS\nchr1\t100\t150\tCNEE_MINUS\n",
                encoding="utf-8",
            )
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                # Same query coordinates, opposite strands, both high identity.
                + "CNEE_PLUS\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tP#1\t0.05\t0.8\t1\n"
                + "CNEE_MINUS\tchr1\t100\t150\t2\t50\tchrQ\t-\t100\t150\tM#1\t0.05\t0.8\t1\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                filter_multi_projections=False,
            )
            rows = _read_tsv(orthology)

        # Overlapping projections on opposite strands must stay separate
        # loci (two one2one calls), not merge into one many2one locus.
        self.assertEqual(result.query_bed_rows, 2)
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in rows
        }
        self.assertEqual(categories_by_query["q_CNEE_PLUS"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_MINUS"], "one2one")

    def test_duplicate_cnee_ids_are_kept_distinct_via_cnee_uid(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from stola.orthology import classify_orthology

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            summary = tmp_path / "summary.tsv"
            orthology = tmp_path / "orthology.tsv"
            query_bed = tmp_path / "query.bed"
            # Two distinct loci that happen to share a display name.
            bed.write_text(
                "chr1\t0\t50\tCNEE_DUP\nchr1\t100\t150\tCNEE_DUP\n",
                encoding="utf-8",
            )
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "ref_chrom",
                        "ref_start",
                        "ref_end",
                        "chain_id",
                        "aligned_bases",
                        "covered_fraction",
                        "gap_bases_within_chain_span",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                        "rank_by_chain_score",
                    ]
                )
                + "\n"
                + "CNEE_DUP\tchr1\t0\t50\t1\t50\t1.0\t0\tchrQ\t+\t100\t150\tD1#1\t0.05\t0.8\t1\n"
                + "CNEE_DUP__line2\tchr1\t100\t150\t2\t50\t1.0\t0\tchrQ\t+\t500\t550\tD2#1\t0.05\t0.8\t1\n",
                encoding="utf-8",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = classify_orthology(
                    stats,
                    bed,
                    tmp_path / "unused_index.json",
                    orthology,
                    query_bed,
                    filter_multi_projections=False,
                )
                summarize_overlap_stats(stats, summary, cnee_bed=bed)

            orthology_rows = _read_tsv(orthology)
            summary_rows = _read_tsv(summary)

        self.assertTrue(
            any("not unique" in str(w.message) for w in caught),
            "expected a duplicate cnee_id warning",
        )
        self.assertEqual(result.query_bed_rows, 2)
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in orthology_rows
        }
        self.assertEqual(categories_by_query["q_CNEE_DUP"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_DUP__line2"], "one2one")
        summary_uids = {row["cnee_uid"] for row in summary_rows}
        self.assertEqual(summary_uids, {"CNEE_DUP", "CNEE_DUP__line2"})

    def test_fasta_accessor_rebuilds_index_when_fasta_changes(self) -> None:
        from stola.genome_io import FastaAccessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta = tmp_path / "genome.fa"
            index = tmp_path / "genome.idx"
            fasta.write_text(">chr1\nAAAAAAAAAA\nCCCCCCCCCC\n", encoding="utf-8")
            accessor = FastaAccessor(fasta, index)
            accessor.close()

            # Same index path, different genome: the old byte offsets would
            # point into the wrong sequence, so the index must be rebuilt.
            fasta.write_text(">chrX\nGGGG\n>chr1\nTTTTT\nTTTTT\nAAAAA\n", encoding="utf-8")
            accessor = FastaAccessor(fasta, index)
            try:
                self.assertEqual(accessor.fetch("chr1", 8, 12), "TTAA")
            finally:
                accessor.close()

            # An index without a #source line (older stola) is rebuilt too.
            index.write_text(
                "#seq_id\tlength\toffset\tline_bases\tline_width\nchr1\t15\t0\t5\t6\n",
                encoding="utf-8",
            )
            accessor = FastaAccessor(fasta, index)
            try:
                self.assertEqual(accessor.fetch("chr1", 0, 15), "T" * 10 + "A" * 5)
                self.assertIn("chrX", accessor.records)
            finally:
                accessor.close()
            self.assertTrue(index.read_text(encoding="utf-8").startswith("#source\t"))

    def test_reverse_complement_matches_biopython(self) -> None:
        from Bio.Seq import Seq

        from stola.genome_io import reverse_complement

        seq = "ACGTNacgtnRYKMBVDHSWUu"
        self.assertEqual(reverse_complement(seq), str(Seq(seq).reverse_complement()))

    def test_vectorized_block_parser_matches_iter_chain_blocks(self) -> None:
        from stola.chain import BlockIndex, RefBlockIndex, iter_chain_blocks, parse_chain_blocks
        from stola.models import ChainBlock

        chains = [
            # Plus-strand query, with gaps on one or both sides.
            "chain 1000 chr1 1000 + 100 400 chrQ 900 + 200 500 1\n50 10 30\n100 20 0\n120\n",
            # Minus-strand query: reference and query run in opposite directions.
            "chain 900 chr1 1000 + 100 400 chrQ 900 - 200 500 2\n50 10 30\n100 20 0\n120\n",
            # A single block.
            "chain 50 chr1 1000 + 10 60 chrQ 900 - 0 50 3\n50\n",
        ]
        for chain_text in chains:
            with self.subTest(chain=chain_text.splitlines()[0]):
                expected = iter_chain_blocks(chain_text)
                query_blocks = parse_chain_blocks(chain_text)
                parsed = [
                    ChainBlock(int(ref_start), int(ref_end), int(query_start), int(query_end))
                    for ref_start, ref_end, query_start, query_end in zip(
                        query_blocks.ref_starts,
                        query_blocks.ref_ends,
                        query_blocks.query_starts,
                        query_blocks.query_ends,
                    )
                ]
                self.assertEqual(parsed, sorted(expected, key=lambda block: block.query_start))

                # The overlap stage's reference-ordered lookup agrees with BlockIndex.
                ref_index = RefBlockIndex.from_query_blocks(query_blocks)
                block_index = BlockIndex.from_blocks(expected)
                for start, end in [(0, 1000), (95, 105), (140, 175), (155, 165), (260, 300), (399, 450)]:
                    self.assertEqual(
                        ref_index.overlapping(start, end), list(block_index.overlapping(start, end))
                    )

    def test_cross_orthology_label_keeps_whole_query_prefix(self) -> None:
        from stola.cross_orthology import collect_cross_run_orthology
        from stola.orthology import ORTHOLOGY_COLUMNS, REF_CNEE_COLUMNS

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            many_refs = ",".join(f"cnee_{index}" for index in range(20_000))
            run_dirs = []
            for name, prefix in (
                ("run_a", "Anas_platyrhynchos_"),
                ("run_b", "Anas_acuta_"),
                ("run_c", "q_"),
            ):
                run_dir = tmp_path / name
                run_dir.mkdir()
                (run_dir / "query_ortho_status.tsv").write_text(
                    "\t".join(ORTHOLOGY_COLUMNS) + "\n"
                    + "NA\tcnee_1\tunaligned\tNA\tNA\n"
                    # A many2one group big enough to pass csv's default field limit.
                    + f"{prefix}many2one_000001\t{many_refs}\tmany2one\t0.9\t0.0\n"
                    # The query id embeds the ref id with ":" cleaned to "_".
                    + f"{prefix}chr1_10-20\tchr1:10-20\tone2one\t0.9\t0.0\n",
                    encoding="utf-8",
                )
                (run_dir / "ref_ortho_status.tsv").write_text(
                    "\t".join(REF_CNEE_COLUMNS) + "\n", encoding="utf-8"
                )
                run_dirs.append((run_dir, None, None))

            columns = collect_cross_run_orthology(run_dirs, build_fasta_data=False)

        # The default prefix falls back to the run directory's name.
        self.assertEqual(
            [column.label for column in columns], ["Anas_platyrhynchos", "Anas_acuta", "run_c"]
        )

    def test_cross_orthology_summary_mode_combines_runs(self) -> None:
        from stola.cross_orthology import (
            _compute_frequency_distributions,
            _read_run_dirs,
            collect_cross_run_orthology,
            main,
        )
        from stola.orthology import ORTHOLOGY_COLUMNS, REF_CNEE_COLUMNS

        categories = {
            "run_a": {"ce1": "one2one", "ce2": "missing/loss", "ce3": "one2many"},
            "run_b": {"ce1": "one2one", "ce2": "one2one", "ce3": "probable_loss"},
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            ref_bed = tmp_path / "ref.bed"
            ref_bed.write_text(
                "chr1\t0\t50\tce1\nchr1\t100\t150\tce2\nchr1\t200\t250\tce3\n", encoding="utf-8"
            )
            for name, category_by_ref in categories.items():
                run_dir = tmp_path / name
                run_dir.mkdir()
                (run_dir / "ref_ortho_status.tsv").write_text(
                    "\t".join(REF_CNEE_COLUMNS) + "\n"
                    + "".join(
                        f"{ref}\t{category}\t1\t0\t{category}\n"
                        for ref, category in category_by_ref.items()
                    ),
                    encoding="utf-8",
                )
                (run_dir / "query_ortho_status.tsv").write_text(
                    "\t".join(ORTHOLOGY_COLUMNS) + "\n" + f"{name}_ce1\tce1\tone2one\t0.9\t0.0\n",
                    encoding="utf-8",
                )
            sheet = tmp_path / "runs.txt"
            # run_a's label comes from its query ids; run_b's from the sheet's
            # third column, with the genome column left blank.
            sheet.write_text(f"{tmp_path / 'run_a'}\n{tmp_path / 'run_b'}\t\tlabel_b\n", encoding="utf-8")
            out_dir = tmp_path / "out"

            with contextlib.redirect_stdout(io.StringIO()):
                main([
                    "summary",
                    "--query-dirs", str(sheet),
                    "--ref-bed", str(ref_bed),
                    "--ref-prefix", "ref",
                    "--out-dir", str(out_dir),
                ])
            combined = (out_dir / "ref_combined_orthology.tsv").read_text(encoding="utf-8").splitlines()
            plots = sorted(path.name for path in out_dir.glob("*.svg"))
            columns = collect_cross_run_orthology(_read_run_dirs(sheet), build_fasta_data=False)

        self.assertEqual(
            combined,
            [
                "ref_ID\trun_a\tlabel_b",
                "ce1\tone2one\tone2one",
                "ce2\tmissing/loss\tone2one",
                "ce3\tone2many\tprobable_loss",
            ],
        )
        self.assertEqual(
            plots, ["ref.loss_frequency_distribution.svg", "ref.multiproj_frequency_distribution.svg"]
        )
        # Index i = number of reference CNEEs lost (or multi-projection) in
        # exactly i of the two runs; missing/loss and probable_loss both count as lost.
        loss_counts, multi_counts = _compute_frequency_distributions(["ce1", "ce2", "ce3"], columns)
        self.assertEqual(loss_counts, [1, 2, 0])
        self.assertEqual(multi_counts, [2, 1, 0])

    def test_fasta_accessor_fetches_correctly_across_wrapped_lines(self) -> None:
        from stola.genome_io import FastaAccessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta = tmp_path / "genome.fa"
            # Uniform 10-base lines with a legitimately short last line.
            fasta.write_text(
                ">chr1\nAAAAAAAAAA\nCCCCCCCCCC\nGGGGGGG\n",
                encoding="utf-8",
            )
            accessor = FastaAccessor(fasta, tmp_path / "genome.idx")
            try:
                self.assertEqual(
                    accessor.fetch("chr1", 0, 27), "A" * 10 + "C" * 10 + "G" * 7
                )
                self.assertEqual(accessor.fetch("chr1", 8, 12), "AACC")
            finally:
                accessor.close()

    def test_fasta_accessor_rejects_ragged_line_wrapping(self) -> None:
        from stola.genome_io import FastaAccessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta = tmp_path / "genome.fa"
            # Line 2 of 3 has 8 bases instead of the established width of 10.
            fasta.write_text(
                ">chr1\nAAAAAAAAAA\nCCCCCCCC\nGGGGGGGGGG\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError) as ctx:
                FastaAccessor(fasta, tmp_path / "genome.idx")

        self.assertIn("ragged", str(ctx.exception))
        self.assertIn("chr1", str(ctx.exception))

    def test_fasta_accessor_matches_plain_parse_and_reads_each_sequence_once(self) -> None:
        from unittest import mock

        from stola.genome_io import FastaAccessor

        rng = random.Random(3)
        for width, newline in [(7, "\n"), (60, "\n"), (61, "\r\n")]:
            with self.subTest(width=width, newline=newline), tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                records = {
                    name: "".join(rng.choice("ACGTNacgtn") for _ in range(length))
                    for name, length in [("chr1", 500), ("chr2", width * 4), ("empty", 0)]
                }
                fasta = tmp_path / "genome.fa"
                with fasta.open("w", encoding="utf-8", newline="") as handle:
                    for name, seq in records.items():
                        handle.write(f">{name} description{newline}")
                        for i in range(0, len(seq), width):
                            handle.write(seq[i : i + width] + newline)

                accessor = FastaAccessor(fasta, tmp_path / "genome.idx")
                accessor.handle = mock.Mock(wraps=accessor.handle)
                try:
                    for name, seq in records.items():
                        for _ in range(200):
                            start = rng.randint(0, len(seq))
                            end = rng.randint(start, len(seq))
                            self.assertEqual(accessor.fetch(name, start, end), seq[start:end].upper())
                    # Each sequence came from a single read of the file.
                    self.assertEqual(accessor.handle.read.call_count, len(records))
                finally:
                    accessor.close()

    def test_fasta_accessor_rejects_blank_line_inside_a_sequence(self) -> None:
        from stola.genome_io import FastaAccessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta = tmp_path / "genome.fa"
            # A blank line would shift every base after it if the record's
            # byte layout were trusted, so it must be an error instead.
            fasta.write_text(">chr1\nAAAAA\n\nCCCCC\nGG\n", encoding="utf-8")
            accessor = FastaAccessor(fasta, tmp_path / "genome.idx")
            try:
                with self.assertRaises(ValueError) as ctx:
                    accessor.fetch("chr1", 0, 5)
            finally:
                accessor.close()

        self.assertIn("blank lines", str(ctx.exception))

    def test_sniff_genome_format_detects_twobit_vs_fasta(self) -> None:
        from stola.genome_io import sniff_genome_format

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            twobit = tmp_path / "genome.2bit"
            fasta = tmp_path / "genome.fa"
            _write_twobit(twobit, {"chr1": "ACGT"})
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")

            self.assertEqual(sniff_genome_format(twobit), "twobit")
            self.assertEqual(sniff_genome_format(fasta), "fasta")

    def test_twobit_accessor_fetches_bases_n_blocks_and_version_alias(self) -> None:
        from stola.genome_io import TwoBitAccessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            twobit = tmp_path / "genome.2bit"
            # Mixed real bases, an N-block, and a soft-masked (lowercase) run
            # -- fetch() always upper-cases, so the mask shouldn't change the
            # returned bases, only exercise that decoding doesn't choke on it.
            _write_twobit(twobit, {"NC_000001.1": "ACGTNNNNacgtACGT"})
            accessor = TwoBitAccessor(twobit)
            try:
                self.assertEqual(accessor.fetch("NC_000001.1", 0, 16), "ACGTNNNNACGTACGT")
                self.assertEqual(accessor.fetch("NC_000001.1", 4, 8), "NNNN")
                # Unversioned alias resolves unambiguously to the one record.
                self.assertEqual(accessor.fetch("NC_000001", 0, 4), "ACGT")
                with self.assertRaises(KeyError):
                    accessor.fetch("no_such_chrom", 0, 1)
                with self.assertRaises(ValueError):
                    accessor.fetch("NC_000001.1", 0, 100)
            finally:
                accessor.close()

    def test_truncated_twobit_is_rejected_when_opened(self) -> None:
        from stola.genome_io import TwoBitAccessor, check_twobit_complete, iter_named_sequences

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            twobit = tmp_path / "genome.2bit"
            # Mimic an interrupted download: the name index at the front of
            # the file survives, chr1 is whole, chr2 is cut mid-sequence, and
            # chr3 starts past end-of-file (which Biopython alone silently
            # reads back as a zero-length sequence). Each record here is 16
            # header bytes + 10 packed bytes, so dropping the last 31 bytes
            # removes all of chr3 and the final 5 bytes of chr2.
            _write_twobit(twobit, {name: "ACGT" * 10 for name in ("chr1", "chr2", "chr3")})
            full = twobit.read_bytes()
            check_twobit_complete(twobit)
            twobit.write_bytes(full[:-31])

            with self.assertRaises(ValueError) as ctx:
                TwoBitAccessor(twobit)
            self.assertIn("truncated", str(ctx.exception))
            self.assertIn("2 of 3 sequences", str(ctx.exception))
            self.assertIn("'chr2'", str(ctx.exception))
            with self.assertRaises(ValueError):
                list(iter_named_sequences(twobit))

            # Cut inside the name index itself (16-byte header + partial index).
            twobit.write_bytes(full[:20])
            with self.assertRaises(ValueError) as ctx:
                check_twobit_complete(twobit)
            self.assertIn("sequence-name index", str(ctx.exception))

    def test_open_genome_accessor_dispatches_on_sniffed_format(self) -> None:
        from stola.genome_io import FastaAccessor, TwoBitAccessor, open_genome_accessor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            twobit = tmp_path / "genome.2bit"
            fasta = tmp_path / "genome.fa"
            _write_twobit(twobit, {"chr1": "ACGT"})
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")

            twobit_accessor = open_genome_accessor(twobit, tmp_path / "unused.idx")
            fasta_accessor = open_genome_accessor(fasta, tmp_path / "genome.idx")
            try:
                self.assertIsInstance(twobit_accessor, TwoBitAccessor)
                self.assertIsInstance(fasta_accessor, FastaAccessor)
                self.assertEqual(twobit_accessor.fetch("chr1", 0, 4), "ACGT")
                self.assertEqual(fasta_accessor.fetch("chr1", 0, 4), "ACGT")
            finally:
                twobit_accessor.close()
                fasta_accessor.close()

    def test_iter_n_runs_matches_between_fasta_and_twobit_for_same_content(self) -> None:
        from stola.assembly_gaps import iter_n_runs

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            twobit = tmp_path / "genome.2bit"
            fasta = tmp_path / "genome.fa"
            records = {
                "chr1": "ACGT" + "N" * 100 + "ACGT" + "N" * 25 + "A",
                "chr2": "N" * 5 + "ACGTACGT",
            }
            _write_twobit(twobit, records)
            _write_fasta(fasta, records)

            fasta_runs = {(r.chrom, r.start, r.end) for r in iter_n_runs(fasta)}
            twobit_runs = {(r.chrom, r.start, r.end) for r in iter_n_runs(twobit)}
            self.assertEqual(fasta_runs, twobit_runs)
            self.assertEqual(
                fasta_runs,
                {("chr1", 4, 104), ("chr1", 108, 133), ("chr2", 0, 5)},
            )

    def test_sequence_comparator_supports_mixed_fasta_and_twobit_genomes(self) -> None:
        try:
            import Bio  # noqa: F401
        except ImportError:
            self.skipTest("Biopython is not installed")
        from stola.chain import iter_chain_blocks
        from stola.sequence_compare import aligned_segments_for_interval

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            ref_fasta = tmp_path / "ref.fa"
            ref_twobit = tmp_path / "ref.2bit"
            query_fasta = tmp_path / "query.fa"
            _write_fasta(ref_fasta, {"chr1": "ACGTACGTAC"})
            _write_twobit(ref_twobit, {"chr1": "ACGTACGTAC"})
            _write_fasta(query_fasta, {"chrQ": "ACGTACGAAC"})  # one mismatch at position 7

            chain_text = "chain 100 chr1 10 + 0 10 chrQ 10 + 0 10 1\n10\n"
            blocks = iter_chain_blocks(chain_text)
            segments = aligned_segments_for_interval(blocks, 0, 10, "+")

            baseline = SequenceComparator(ref_fasta, query_fasta, tmp_path / "work_fasta")
            mixed = SequenceComparator(ref_twobit, query_fasta, tmp_path / "work_mixed")
            try:
                baseline_metrics = baseline.compare("chr1", 0, 10, "chrQ", "+", segments)
                mixed_metrics = mixed.compare("chr1", 0, 10, "chrQ", "+", segments)
            finally:
                baseline.close()
                mixed.close()

            self.assertEqual(mixed_metrics, baseline_metrics)
            self.assertEqual(baseline_metrics["seq_mismatches"], 1)

    def test_sequence_comparison_metrics_for_plus_strand_projection(self) -> None:
        try:
            import Bio  # noqa: F401
        except ImportError:
            self.skipTest("Biopython is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            ref_fasta = tmp_path / "ref.fa"
            query_fasta = tmp_path / "query.fa"
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            ref_fasta.write_text(
                ">chr1.1\n" + ("N" * 120) + ("A" * 120) + ("N" * 100) + "\n"
                ">chr2\n" + ("C" * 400) + "\n",
                encoding="utf-8",
            )
            query_fasta.write_text(
                ">chrQ.1\n"
                + ("N" * 520)
                + ("A" * 30)
                + ("G" * 40)
                + ("N" * 500)
                + "\n"
                ">chrQ2\n"
                + ("N" * 770)
                + ("T" * 120)
                + ("N" * 200)
                + "\n"
                ">chrQ3\n"
                + ("N" * 1000)
                + "\n",
                encoding="utf-8",
            )
            build_index(chain, index)
            write_candidates(index, bed, candidates)
            comparator = SequenceComparator(ref_fasta, query_fasta, tmp_path)
            try:
                from stola.overlap import iter_overlap_rows

                seq_rows = list(
                    iter_overlap_rows(
                        index,
                        candidates,
                        sequence_comparator=comparator,
                    )
                )
            finally:
                comparator.close()

        by_pair = {(row["cnee_id"], row["chain_id"]): row for row in seq_rows}
        plus = by_pair[("CNEE_PARTIAL", "1")]
        self.assertEqual(plus["seq_ref_len"], 120)
        self.assertEqual(plus["seq_query_len"], 70)
        self.assertEqual(plus["seq_compared_bases"], 70)
        self.assertEqual(plus["seq_matches"], 30)
        self.assertEqual(plus["seq_mismatches"], 40)
        self.assertAlmostEqual(float(plus["seq_identity"]), 30 / 70)

        minus = by_pair[("CNEE_PARTIAL", "2")]
        self.assertEqual(minus["seq_compared_bases"], 120)
        self.assertEqual(minus["seq_matches"], 120)
        self.assertEqual(minus["seq_mismatches"], 0)
        self.assertEqual(minus["seq_identity"], 1.0)

        gapped = by_pair[("CNEE_GAPPED", "3")]
        self.assertEqual(gapped["seq_compared_bases"], 0)
        self.assertEqual(gapped["seq_identity"], ".")


if __name__ == "__main__":
    unittest.main()
