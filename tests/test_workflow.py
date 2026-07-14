import csv
import gzip
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from toga4cnee.candidates import write_candidates
from toga4cnee.index import build_index, load_index
from toga4cnee.overlap import calculate_overlap_rows, write_overlap_stats
from toga4cnee.sequence_compare import SequenceComparator
from toga4cnee.summary import summarize_overlap_stats
from toga4cnee.workflow import run_workflow


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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["num_query_alignments"], "0")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["overlap_fractions"], "")

    def test_single_workflow_writes_stats_and_query_bed(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            stats = tmp_path / "overlap_stats.tsv"
            query_bed = tmp_path / "query_cnees.bed"
            summary = tmp_path / "summary.tsv"
            orthology = tmp_path / "orthology.tsv"
            work_dir = tmp_path / "work"

            result = run_workflow(
                chain,
                bed,
                query_bed,
                stats,
                out_summary=summary,
                out_orthology=orthology,
                work_dir=work_dir,
            )
            stats_rows = _read_tsv(stats)
            bed_lines = query_bed.read_text(encoding="utf-8").splitlines()
            summary_rows = _read_tsv(summary)
            orthology_rows = _read_tsv(orthology)
            strict_stats_rows = _read_tsv(tmp_path / "overlap_stats.strict.tsv")
            strict_bed_lines = (tmp_path / "query_cnees.strict.bed").read_text(
                encoding="utf-8"
            ).splitlines()
            strict_summary_rows = _read_tsv(
                tmp_path / "overlap_stats.strict.summary.tsv"
            )
            strict_orthology_rows = _read_tsv(
                tmp_path / "overlap_stats.strict.orthology.tsv"
            )
            synteny_rows = _read_tsv(tmp_path / "overlap_stats.synteny.tsv")

        self.assertEqual(result.chain_count, 3)
        self.assertEqual(result.candidate_count, 3)
        self.assertEqual(result.stats_count, 3)
        self.assertEqual(result.strict_stats_count, 1)
        self.assertEqual(result.query_bed_count, 2)
        self.assertEqual(result.strict_query_bed_count, 1)
        self.assertEqual(result.summary_count, 3)
        self.assertEqual(result.strict_summary_count, 3)
        self.assertEqual(result.orthology_result.orthology_rows, 4)
        self.assertEqual(result.strict_orthology_result.orthology_rows, 3)
        # 2 synteny candidates x 3 default flank sizes (10000, 1000, 100).
        self.assertEqual(result.synteny_result.rows, 6)
        self.assertEqual(len(stats_rows), 3)
        self.assertEqual(len(strict_stats_rows), 1)
        self.assertEqual(len(summary_rows), 3)
        self.assertEqual(len(strict_summary_rows), 3)
        self.assertEqual(len(orthology_rows), 4)
        self.assertEqual(len(strict_orthology_rows), 3)
        self.assertEqual(len(synteny_rows), 6)
        self.assertEqual(
            list(orthology_rows[0]),
            ["query_cnee_id", "partner_ref_cnees", "orthology_category"],
        )
        self.assertEqual(
            bed_lines,
            [
                "chrQ\t520\t590\tq_CNEE_PARTIAL_1\t0\t+",
                "chrQ2\t770\t890\tq_CNEE_PARTIAL_2\t0\t-",
            ],
        )
        self.assertEqual(
            strict_bed_lines,
            [
                "chrQ2\t770\t890\tq_CNEE_PARTIAL\t0\t-",
            ],
        )
        categories = {
            row["query_cnee_id"]: row["orthology_category"] for row in orthology_rows
        }
        self.assertEqual(categories["q_CNEE_PARTIAL_1"], "one2many")
        self.assertEqual(categories["q_CNEE_PARTIAL_2"], "one2many")
        strict_categories = {
            row["query_cnee_id"]: row["orthology_category"]
            for row in strict_orthology_rows
        }
        self.assertEqual(strict_categories["q_CNEE_PARTIAL"], "one2one")
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in orthology_rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_GAPPED"], "lost")
        self.assertEqual(na_categories["CNEE_NO_CHAIN"], "missing")
        summary_statuses = {
            row["cnee_uid"]: row["projection_status"] for row in summary_rows
        }
        self.assertEqual(summary_statuses["CNEE_GAPPED"], "lost")
        self.assertEqual(summary_statuses["CNEE_NO_CHAIN"], "missing")
        strict_summary_by_id = {
            row["cnee_uid"]: row for row in strict_summary_rows
        }
        self.assertEqual(
            strict_summary_by_id["CNEE_PARTIAL"]["num_query_alignments"], "1"
        )
        self.assertEqual(strict_summary_by_id["CNEE_PARTIAL"]["projection_status"], "aligned")
        synteny_statuses = {row["orthology_status"] for row in synteny_rows}
        self.assertEqual(synteny_statuses, {"weak", "lost"})
        self.assertEqual(
            {int(row["flank_size"]) for row in synteny_rows}, {10000, 1000, 100}
        )
        self.assertIn(
            "q_CNEE_PARTIAL_1",
            {row["cnee_uid"] for row in synteny_rows if row["orthology_status"] == "weak"},
        )

    def test_orthology_graph_edges_store_per_projection_metrics_and_require_single_projection_solid(
        self,
    ) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        import networkx as nx

        from toga4cnee.orthology import (
            _build_orthology_graph,
            _cluster_query_projections,
            _edge_is_solid,
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
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                    ]
                )
                + "\n"
                # Good identity but bad gap on its own.
                + "CNEE_A\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\n"
                # Good gap but bad identity on its own.
                + "CNEE_A\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\n"
                + "CNEE_B\t40\tchrQ\t+\t300\t340\tB#1\t0.25\t.\n",
                encoding="utf-8",
            )

            evidence = _read_projection_evidence(stats)
            clusters = _cluster_query_projections(evidence.projections)
            graph = _build_orthology_graph(nx, ["CNEE_A", "CNEE_B"], clusters)

        edge_a = graph.edges[_ref_node("CNEE_A"), "query:000001"]
        self.assertEqual(edge_a["projection_metrics"], [(0.9, 0.5), (0.3, 0.05)])
        self.assertEqual(edge_a["projection_count"], 2)
        # Blending would give (max identity 0.9, min gap 0.05), which passes
        # both thresholds -- but no single projection does, so the edge must
        # not be solid.
        self.assertFalse(
            _edge_is_solid(edge_a, min_sequence_identity=0.7, max_fraction_gap=0.1)
        )

        edge_b = graph.edges[_ref_node("CNEE_B"), "query:000002"]
        self.assertEqual(edge_b["projection_metrics"], [(None, 0.25)])
        self.assertEqual(edge_b["projection_count"], 1)

    def test_orthology_classification_uses_solid_edges_only(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from toga4cnee.orthology import classify_orthology

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
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                    ]
                )
                + "\n"
                + "CNEE_SOLID\t50\tchrQ\t+\t100\t150\tS#1\t0.05\t0.8\n"
                + "CNEE_SEQ_MISSING\t50\tchrQ\t+\t300\t350\tM#1\t0.05\t.\n"
                + "CNEE_NONSOLID\t50\tchrQ\t+\t500\t550\tN#1\t0.2\t0.5\n"
                + "CNEE_GAP_ONLY\t0\tchrQ\t+\t.\t.\tG#1\t1.0\t.\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                orthology,
                query_bed,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
            )
            rows = _read_tsv(orthology)
            bed_lines = query_bed.read_text(encoding="utf-8").splitlines()

        self.assertEqual(result.orthology_rows, 5)
        self.assertEqual(result.query_bed_rows, 2)
        self.assertEqual(
            bed_lines,
            [
                "chrQ\t100\t150\tq_CNEE_SOLID\t0\t+",
                "chrQ\t300\t350\tq_CNEE_SEQ_MISSING\t0\t+",
            ],
        )
        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in rows
        }
        self.assertEqual(categories_by_query["q_CNEE_SOLID"], "one2one")
        self.assertEqual(categories_by_query["q_CNEE_SEQ_MISSING"], "one2one")
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_NONSOLID"], "uncertain_loss")
        self.assertEqual(na_categories["CNEE_GAP_ONLY"], "lost")
        self.assertEqual(na_categories["CNEE_MISSING"], "missing")

    def test_orthology_falls_back_to_uncertain_loss_when_no_single_projection_is_solid(
        self,
    ) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from toga4cnee.orthology import classify_orthology

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
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                    ]
                )
                + "\n"
                # Good identity, bad gap on its own.
                + "CNEE_A\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\n"
                # Good gap, bad identity on its own. Same query cluster as A#1
                # (overlapping query coordinates) so both land on one edge.
                + "CNEE_A\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                orthology,
                query_bed,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
            )
            rows = _read_tsv(orthology)

        # Blending would give (max identity 0.9, min gap 0.05), which passes
        # both thresholds. Neither projection alone does, so CNEE_A must not
        # be called orthologous -- it falls back to uncertain_loss.
        self.assertEqual(result.query_bed_rows, 0)
        self.assertEqual(result.orthology_rows, 1)
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_A"], "uncertain_loss")

    def test_orthology_clustering_keeps_opposite_strand_overlaps_separate(self) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from toga4cnee.orthology import classify_orthology

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
                        "aligned_bases",
                        "query_chrom",
                        "query_strand",
                        "query_aligned_start",
                        "query_aligned_end",
                        "query_projection_name",
                        "gap_fraction_within_chain_span",
                        "seq_identity",
                    ]
                )
                + "\n"
                # Same query coordinates, opposite strands, both solid.
                + "CNEE_PLUS\t50\tchrQ\t+\t100\t150\tP#1\t0.05\t0.8\n"
                + "CNEE_MINUS\t50\tchrQ\t-\t100\t150\tM#1\t0.05\t0.8\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                orthology,
                query_bed,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
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
        from toga4cnee.orthology import classify_orthology

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
                    ]
                )
                + "\n"
                + "CNEE_DUP\t50\t1.0\t0\tchrQ\t+\t100\t150\tD1#1\t0.05\t0.8\n"
                + "CNEE_DUP__line2\t50\t1.0\t0\tchrQ\t+\t500\t550\tD2#1\t0.05\t0.8\n",
                encoding="utf-8",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = classify_orthology(
                    stats,
                    bed,
                    orthology,
                    query_bed,
                    min_sequence_identity=0.7,
                    max_fraction_gap=0.1,
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

    def test_synteny_verification_writes_long_format_rows_per_custom_flank_size(
        self,
    ) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from toga4cnee.synteny import write_synteny_verification

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = _write_inputs(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"
            query_bed = tmp_path / "query_cnees.bed"
            orthology = tmp_path / "orthology.tsv"
            synteny = tmp_path / "synteny.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            write_overlap_stats(index, candidates, stats)
            from toga4cnee.orthology import classify_orthology

            classify_orthology(stats, bed, orthology, query_bed, output_mode="all")

            result = write_synteny_verification(
                index,
                stats,
                query_bed,
                orthology,
                synteny,
                flank_sizes=(500, 50),
            )
            rows = _read_tsv(synteny)

        self.assertEqual(result.rows, len(rows))
        self.assertGreater(result.rows, 0)
        # Every candidate gets a row for each requested flank size, and only
        # the requested sizes -- not the module default.
        self.assertEqual(result.rows % 2, 0)
        self.assertEqual({int(row["flank_size"]) for row in rows}, {500, 50})
        cnee_uids = {row["cnee_uid"] for row in rows}
        for cnee_uid in cnee_uids:
            sizes_for_uid = {
                int(row["flank_size"]) for row in rows if row["cnee_uid"] == cnee_uid
            }
            self.assertEqual(sizes_for_uid, {500, 50})

    def test_synteny_selection_caps_one2one_controls(self) -> None:
        from toga4cnee.synteny import _select_synteny_targets

        rows = [
            {
                "cnee_id": "CNEE_A",
                "cnee_uid": "CNEE_A",
                "aligned_bases": "10",
                "gap_fraction_within_chain_span": "0.0",
                "seq_identity": "1.0",
                "query_chrom": "chrQ",
                "query_aligned_start": "10",
                "query_aligned_end": "20",
            },
            {
                "cnee_id": "CNEE_B",
                "cnee_uid": "CNEE_B",
                "aligned_bases": "10",
                "gap_fraction_within_chain_span": "0.0",
                "seq_identity": "1.0",
                "query_chrom": "chrQ",
                "query_aligned_start": "30",
                "query_aligned_end": "40",
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            query_bed = tmp_path / "query.bed"
            orthology = tmp_path / "orthology.tsv"
            query_bed.write_text(
                "chrQ\t10\t20\tq_CNEE_A\t0\t+\n"
                "chrQ\t30\t40\tq_CNEE_B\t0\t+\n",
                encoding="utf-8",
            )
            orthology.write_text(
                "query_cnee_id\tpartner_ref_cnees\torthology_category\n"
                "q_CNEE_A\tCNEE_A\tone2one\n"
                "q_CNEE_B\tCNEE_B\tone2one\n",
                encoding="utf-8",
            )

            targets = _select_synteny_targets(
                rows,
                query_bed,
                orthology,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
                max_one2one_controls=1,
            )

        self.assertEqual(len(targets), 1)
        self.assertEqual(targets[0].status, "one2one")
        self.assertEqual(targets[0].output_id, "q_CNEE_A")

    def test_orthology_diagnostic_plots_use_category_groups(self) -> None:
        from toga4cnee.plots import write_orthology_diagnostic_plots

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
        from toga4cnee.plots import write_multi_alignment_metric_densities

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
        from toga4cnee.sequence_compare import FastaAccessor

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
        from toga4cnee.sequence_compare import FastaAccessor

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
                from toga4cnee.overlap import iter_overlap_rows

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
