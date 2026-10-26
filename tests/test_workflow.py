import csv
import gzip
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
from stola.summary import summarize_overlap_stats
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

            payload = build_index(chain, index)
            _, records = load_index(index)

        self.assertEqual(payload["chain_count"], 3)
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

            payload = build_index(gz_chain, index)
            _, records = load_index(index)

            self.assertEqual(payload["chain_count"], 3)
            self.assertEqual(Path(payload["indexed_chain_path"]).name, "test.chain")
            self.assertEqual(
                Path(payload["indexed_chain_path"]).read_text(encoding="utf-8"),
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

    def test_summary_from_overlap_stats_uses_aligned_rows_by_default(self) -> None:
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
        self.assertEqual(rows[0]["cnee_uid"], "CNEE_PARTIAL")
        self.assertEqual(rows[0]["projection_status"], "aligned")
        self.assertEqual(rows[0]["num_query_alignments"], "2")
        self.assertEqual(rows[0]["overlap_fractions"], "1.0,0.5833333333333334")
        self.assertEqual(rows[0]["gap_bases"], "0,50")
        self.assertEqual(rows[1]["cnee_uid"], "CNEE_GAPPED")
        self.assertEqual(rows[1]["projection_status"], "lost")
        self.assertEqual(rows[1]["num_query_alignments"], "0")
        self.assertEqual(rows[1]["overlap_fractions"], "")
        self.assertEqual(rows[1]["gap_bases"], "")
        self.assertEqual(rows[2]["cnee_uid"], "CNEE_NO_CHAIN")
        self.assertEqual(rows[2]["projection_status"], "missing")
        self.assertEqual(rows[2]["num_query_alignments"], "0")

    def test_solid_aware_summary_marks_weak_only_refs_uncertain_loss(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bed = tmp_path / "cnees.bed"
            stats = tmp_path / "stats.tsv"
            summary = tmp_path / "summary.tsv"
            bed.write_text(
                "chr1\t0\t50\tCNEE_SOLID_PLUS_WEAK\n"
                "chr1\t100\t150\tCNEE_WEAK_ONLY\n",
                encoding="utf-8",
            )
            stats.write_text(
                "cnee_uid\taligned_bases\tcovered_fraction\t"
                "gap_bases_within_chain_span\tgap_fraction_within_chain_span\t"
                "seq_identity\n"
                "CNEE_SOLID_PLUS_WEAK\t50\t1.0\t0\t0.05\t0.8\n"
                "CNEE_SOLID_PLUS_WEAK\t25\t0.5\t25\t0.5\t0.8\n"
                "CNEE_WEAK_ONLY\t25\t0.5\t25\t0.5\t0.8\n",
                encoding="utf-8",
            )

            count = summarize_overlap_stats(
                stats,
                summary,
                cnee_bed=bed,
                require_solid_reference=True,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
            )
            rows = {row["cnee_uid"]: row for row in _read_tsv(summary)}

        self.assertEqual(count, 2)
        self.assertEqual(rows["CNEE_SOLID_PLUS_WEAK"]["projection_status"], "aligned")
        self.assertEqual(rows["CNEE_SOLID_PLUS_WEAK"]["num_query_alignments"], "2")
        self.assertEqual(rows["CNEE_SOLID_PLUS_WEAK"]["overlap_fractions"], "1.0,0.5")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["projection_status"], "uncertain_loss")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["num_query_alignments"], "1")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["overlap_fractions"], "0.5")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["gap_bases"], "25")

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
            work_files = sorted(path.name for path in (tmp_path / "work").iterdir())
            summary_figures = sorted(path.name for path in (tmp_path / "summary_figures").iterdir())

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
        summary_statuses = {
            row["cnee_uid"]: row["projection_status"] for row in summary_rows
        }
        self.assertEqual(summary_statuses["CNEE_GAPPED"], "lost")
        self.assertEqual(summary_statuses["CNEE_NO_CHAIN"], "missing")
        # The synteny reports are no longer part of the pipeline.
        self.assertFalse([name for name in work_files if "synteny" in name])

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
        import networkx as nx

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
            graph = _build_orthology_graph(nx, ["CNEE_A", "CNEE_B"], clusters)

        edge_a = graph.edges[_ref_node("CNEE_A"), "query:000001"]
        self.assertEqual(edge_a["projection_metrics"], [(0.9, 0.5), (0.3, 0.05)])
        self.assertEqual(edge_a["projection_count"], 2)

        edge_b = graph.edges[_ref_node("CNEE_B"), "query:000002"]
        self.assertEqual(edge_b["projection_metrics"], [(None, 0.25)])
        self.assertEqual(edge_b["projection_count"], 1)

    def test_orthology_classification_uses_solid_edges_only(self) -> None:
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
                        "chr1\t0\t50\tCNEE_SOLID",
                        "chr1\t100\t150\tCNEE_SEQ_MISSING",
                        "chr1\t200\t250\tCNEE_NONSOLID",
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
                + "CNEE_SOLID\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tS#1\t0.05\t0.8\t1\n"
                + "CNEE_SEQ_MISSING\tchr1\t100\t150\t2\t50\tchrQ\t+\t300\t350\tM#1\t0.05\t.\t1\n"
                + "CNEE_NONSOLID\tchr1\t200\t250\t3\t50\tchrQ\t+\t500\t550\tN#1\t0.2\t0.5\t1\n"
                + "CNEE_GAP_ONLY\tchr1\t300\t350\t4\t0\tchrQ\t+\t.\t.\tG#1\t1.0\t.\t1\n",
                encoding="utf-8",
            )

            # synteny_check=False: this test is about the solid/non-solid
            # distinction driving categorization, not the verification feature
            # itself (covered separately), so pruning is disabled and the
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
        # CNEE_NONSOLID is now a real (weak) one2one query locus rather than
        # a standalone uncertain_loss row -- it still gets a query CNEE id.
        self.assertEqual(result.query_bed_rows, 3)
        self.assertEqual(
            bed_lines,
            [
                "chrQ\t100\t150\tq_CNEE_SOLID\t0\t+",
                "chrQ\t300\t350\tq_CNEE_SEQ_MISSING\t0\t+",
                "chrQ\t500\t550\tq_CNEE_NONSOLID\t0\t+",
            ],
        )
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in rows
        }
        self.assertEqual(categories_by_query["q_CNEE_SOLID"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_SEQ_MISSING"], "one2one")
        # A lone ref/query edge is always reported as "one2one" now, even if
        # non-solid -- the one2one(weak) label only ever applied to the raw,
        # pre-verification graph shape (see out_culled_projections).
        self.assertEqual(categories_by_query["q_CNEE_NONSOLID"], "one2one")
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_GAP_ONLY"], "lost")
        self.assertEqual(na_categories["CNEE_MISSING"], "unaligned")

    def test_orthology_keeps_lone_edge_as_one2one_when_no_single_projection_is_solid(
        self,
    ) -> None:
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
                # Good identity, bad gap on its own.
                + "CNEE_A\tchr1\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\t1\n"
                # Good gap, bad identity on its own. Same query cluster as A#1
                # (overlapping query coordinates) so both land on one edge.
                + "CNEE_A\tchr1\t0\t50\t1\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\t2\n",
                encoding="utf-8",
            )

            # synteny_check=False: no other CNEEs exist to verify against, and
            # this test is specifically about the solid/non-solid distinction.
            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                filter_multi_projections=False,
            )
            rows = _read_tsv(orthology)

        # Blending would give (max identity 0.9, min gap 0.05), which passes
        # both thresholds. Neither projection alone does, so the edge isn't
        # solid -- but it's still the only ref/query pair here, so it's kept
        # as a one2one call (the weak/solid distinction only matters for the
        # raw pre-verification category, not final output) rather than
        # dropped to uncertain_loss.
        self.assertEqual(result.query_bed_rows, 1)
        self.assertEqual(result.orthology_rows, 1)
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in rows
        }
        self.assertEqual(categories_by_query["q_CNEE_A"], "one2one")

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
                # One clean solid one2one.
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
                # Same query coordinates, opposite strands, both solid.
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

    def test_orthology_diagnostic_plots_use_category_groups(self) -> None:
        from stola.plots import write_orthology_diagnostic_plots

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            stats = tmp_path / "stats.tsv"
            orthology = tmp_path / "orthology.tsv"
            category_plot = tmp_path / "categories.svg"
            seq_plot = tmp_path / "seq.svg"
            gap_plot = tmp_path / "gap.svg"
            stats.write_text(
                "\t".join(
                    [
                        "cnee_uid",
                        "aligned_bases",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                    ]
                )
                + "\n"
                + "CNEE_SOLID\t50\t0.05\t0.8\n"
                + "CNEE_NONSOLID\t50\t0.2\t0.5\n"
                + "CNEE_LOST\t0\t1.0\t.\n",
                encoding="utf-8",
            )
            orthology.write_text(
                "\t".join(
                    [
                        "query_cnee_id",
                        "partner_ref_cnees",
                        "orthology_category",
                    ]
                )
                + "\n"
                + "q_CNEE_SOLID\tCNEE_SOLID\tone2one\n"
                + "NA\tCNEE_NONSOLID\tuncertain_loss\n"
                + "NA\tCNEE_LOST\tlost\n",
                encoding="utf-8",
            )

            result = write_orthology_diagnostic_plots(
                stats,
                orthology,
                category_plot,
                seq_plot,
                gap_plot,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
            )
            category_text = category_plot.read_text(encoding="utf-8")
            seq_text = seq_plot.read_text(encoding="utf-8")
            gap_text = gap_plot.read_text(encoding="utf-8")

        self.assertEqual(result.category_counts["one2one"], 1)
        self.assertEqual(result.category_counts["uncertain_loss"], 1)
        self.assertEqual(result.category_counts["lost"], 1)
        self.assertEqual(result.sequence_identity_counts["one2one"], 1)
        self.assertEqual(result.sequence_identity_counts["uncertain_loss"], 1)
        self.assertEqual(result.gap_fraction_counts["lost"], 1)
        self.assertIn("Orthology category counts", category_text)
        self.assertIn("Sequence identity by orthology category", seq_text)
        self.assertIn("Gap fraction by orthology category", gap_text)

    def test_multi_alignment_metric_plots_filter_by_summary_counts(self) -> None:
        from stola.plots import write_multi_alignment_metric_densities

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            summary = tmp_path / "summary.tsv"
            stats = tmp_path / "stats.tsv"
            seq_plot = tmp_path / "multi_seq.svg"
            gap_plot = tmp_path / "multi_gap.svg"
            summary.write_text(
                "cnee_uid\tnum_query_alignments\n"
                "CNEE_SINGLE\t1\n"
                "CNEE_MULTI\t2\n"
                "CNEE_ZERO\t0\n",
                encoding="utf-8",
            )
            stats.write_text(
                "cnee_uid\tgap_fraction_within_chain_span\tseq_identity\n"
                "CNEE_SINGLE\t0.05\t0.9\n"
                "CNEE_MULTI\t0.1\t0.8\n"
                "CNEE_MULTI\t0.4\t.\n"
                "CNEE_ZERO\t1.0\t.\n",
                encoding="utf-8",
            )

            result = write_multi_alignment_metric_densities(
                summary,
                stats,
                seq_plot,
                gap_plot,
            )
            seq_text = seq_plot.read_text(encoding="utf-8")
            gap_text = gap_plot.read_text(encoding="utf-8")

        self.assertEqual(result.multi_alignment_ref_count, 1)
        self.assertEqual(result.sequence_identity_density_n, 1)
        self.assertEqual(result.gap_fraction_density_n, 2)
        self.assertIn("Sequence identity density for multi-query CNEEs", seq_text)
        self.assertIn("Gap fraction density for multi-query CNEEs", gap_text)

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
