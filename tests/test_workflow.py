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
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["num_query_alignments"], "1")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["overlap_fractions"], "0.5")
        self.assertEqual(rows["CNEE_WEAK_ONLY"]["gap_bases"], "25")

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
            strict_summary_rows = _read_tsv(
                tmp_path / "overlap_stats.strict.summary.tsv"
            )
            synteny_rows = _read_tsv(tmp_path / "overlap_stats.synteny.tsv")

        self.assertEqual(result.chain_count, 3)
        self.assertEqual(result.candidate_count, 3)
        self.assertEqual(result.stats_count, 3)
        self.assertEqual(result.strict_stats_count, 1)
        # CNEE_PARTIAL spans two candidate chains (chrQ, chrQ2) -- an initial
        # one2many shape needing verification. It's also alone on chr1, so it
        # has no reference-CNEE neighbor on either flank at all -- under v2's
        # neighbor_fraction criterion, no trusted neighbor means no synteny
        # evidence either way, which fails outright regardless of align_score.
        # BOTH projections get pruned (not just the weak chrQ one), so this
        # ref ends up uncertain_loss with zero query CNEE loci, rather than
        # keeping the solid chrQ2 projection as a lone one2one.
        self.assertEqual(result.query_bed_count, 0)
        self.assertEqual(result.summary_count, 3)
        self.assertEqual(result.strict_summary_count, 3)
        self.assertEqual(result.orthology_result.orthology_rows, 3)
        self.assertEqual(result.orthology_result.culled_projection_rows, 2)
        # 2 synteny candidates (weak + lost; CNEE_PARTIAL no longer has a
        # one2one control row since it isn't one2one anymore) x 3 default
        # flank sizes (10000, 1000, 100).
        self.assertEqual(result.synteny_result.rows, 6)
        self.assertEqual(len(stats_rows), 3)
        self.assertEqual(len(strict_stats_rows), 1)
        self.assertEqual(len(summary_rows), 3)
        self.assertEqual(len(strict_summary_rows), 3)
        self.assertEqual(len(orthology_rows), 3)
        self.assertEqual(len(synteny_rows), 6)
        self.assertEqual(
            list(orthology_rows[0]),
            ["query_cnee_id", "partner_ref_cnees", "orthology_category"],
        )
        self.assertEqual(bed_lines, [])
        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in orthology_rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CNEE_PARTIAL"], "uncertain_loss")
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
        # The pruned (non-solid, unverifiable) chrQ projection never got a
        # graph-assigned query id, so it's reported under its raw projection name.
        self.assertIn(
            "CNEE_PARTIAL#1",
            {row["cnee_uid"] for row in synteny_rows if row["orthology_status"] == "weak"},
        )
        # Every CNEE in this fixture is on its own chromosome.
        for row in synteny_rows:
            self.assertEqual(row["distance_to_nearest_cnee_neighbor"], ".")
        for row in synteny_rows:
            if row["orthology_status"] == "lost":
                self.assertEqual(row["distance_to_nearest_aligned_block"], "25")

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
                    ]
                )
                + "\n"
                # Good identity but bad gap on its own.
                + "CNEE_A\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\n"
                # Good gap but bad identity on its own.
                + "CNEE_A\t0\t50\t1\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\n"
                + "CNEE_B\t100\t150\t2\t40\tchrQ\t+\t300\t340\tB#1\t0.25\t.\n",
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
                    ]
                )
                + "\n"
                + "CNEE_SOLID\t0\t50\t1\t50\tchrQ\t+\t100\t150\tS#1\t0.05\t0.8\n"
                + "CNEE_SEQ_MISSING\t100\t150\t2\t50\tchrQ\t+\t300\t350\tM#1\t0.05\t.\n"
                + "CNEE_NONSOLID\t200\t250\t3\t50\tchrQ\t+\t500\t550\tN#1\t0.2\t0.5\n"
                + "CNEE_GAP_ONLY\t300\t350\t4\t0\tchrQ\t+\t.\t.\tG#1\t1.0\t.\n",
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
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
                synteny_check=False,
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
        self.assertEqual(na_categories["CNEE_MISSING"], "missing")

    def test_orthology_keeps_lone_edge_as_one2one_when_no_single_projection_is_solid(
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
                    ]
                )
                + "\n"
                # Good identity, bad gap on its own.
                + "CNEE_A\t0\t50\t1\t50\tchrQ\t+\t100\t150\tA#1\t0.5\t0.9\n"
                # Good gap, bad identity on its own. Same query cluster as A#1
                # (overlapping query coordinates) so both land on one edge.
                + "CNEE_A\t0\t50\t1\t45\tchrQ\t+\t120\t180\tA#2\t0.05\t0.3\n",
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
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
                synteny_check=False,
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

    def test_synteny_verification_v2_promotes_syntenic_weak_call_and_culls_mismatched_one(
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
            culled = tmp_path / "culled.tsv"
            # A/B/D/E/F/G are solid one2one anchors, colinear on both genomes.
            # CAND_PASS is weak but sits exactly where a syntenic candidate
            # should (between B and D on both ref and query) -- it should
            # survive verification and merge into plain "one2one" output.
            # CAND_FAIL is also weak, bracketed by trusted anchors F/G, but
            # its own query projection sits nowhere near F or G's query loci
            # (a rearranged/spurious hit) -- it should get culled.
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t0\t50\tA",
                        "chr1\t100\t150\tB",
                        "chr1\t200\t250\tCAND_PASS",
                        "chr1\t300\t350\tD",
                        "chr1\t400\t450\tE",
                        "chr1\t500\t550\tF",
                        "chr1\t600\t650\tCAND_FAIL",
                        "chr1\t700\t750\tG",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            columns = [
                "cnee_uid",
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
            ]
            rows = [
                ("A", 0, 50, 1000, 1050, "A#1", 0.05, 0.9),
                ("B", 100, 150, 1100, 1150, "B#1", 0.05, 0.9),
                ("CAND_PASS", 200, 250, 1200, 1250, "CANDP#1", 0.5, 0.9),
                ("D", 300, 350, 1300, 1350, "D#1", 0.05, 0.9),
                ("E", 400, 450, 1400, 1450, "E#1", 0.05, 0.9),
                ("F", 500, 550, 1500, 1550, "F#1", 0.05, 0.9),
                # Rearranged: lands far from F/G's query neighborhood instead
                # of between them.
                ("CAND_FAIL", 600, 650, 50, 100, "CANDF#1", 0.5, 0.9),
                ("G", 700, 750, 1700, 1750, "G#1", 0.05, 0.9),
            ]
            stats.write_text(
                "\t".join(columns)
                + "\n"
                + "\n".join(
                    "\t".join(
                        [
                            ref_id,
                            str(ref_start),
                            str(ref_end),
                            "1",
                            "50",
                            "chrQ",
                            "+",
                            str(q_start),
                            str(q_end),
                            proj_name,
                            str(gap),
                            str(identity),
                        ]
                    )
                    for ref_id, ref_start, ref_end, q_start, q_end, proj_name, gap, identity in rows
                )
                + "\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
                synteny_check=True,
                verification_num_neighbors=1,
                out_culled_projections=culled,
            )
            orthology_rows = _read_tsv(orthology)
            culled_rows = _read_tsv(culled)

        categories_by_query = {
            row["query_cnee_id"]: row["orthology_category"] for row in orthology_rows
        }
        # CAND_PASS survived and is reported as plain "one2one", not
        # "one2one(weak)".
        self.assertEqual(categories_by_query["q_CAND_PASS"], "one2one")
        for ref_id in ("A", "B", "D", "E", "F", "G"):
            self.assertEqual(categories_by_query[f"q_{ref_id}"], "one2one")

        na_categories = {
            row["partner_ref_cnees"]: row["orthology_category"]
            for row in orthology_rows
            if row["query_cnee_id"] == "NA"
        }
        self.assertEqual(na_categories["CAND_FAIL"], "uncertain_loss")
        self.assertEqual(result.query_bed_rows, 7)
        self.assertEqual(result.culled_projection_rows, 1)

        self.assertEqual(len(culled_rows), 1)
        self.assertEqual(
            culled_rows[0],
            {
                "projection_name": "CANDF#1",
                "ref_cnee_id": "CAND_FAIL",
                "initial_orthology_category": "one2one(weak)",
                "query_chrom": "chrQ",
                "query_start": "50",
                "query_end": "100",
                "sequence_identity": "0.9",
                "gap_fraction": "0.5",
            },
        )

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
                    ]
                )
                + "\n"
                # Same query coordinates, opposite strands, both solid.
                + "CNEE_PLUS\t0\t50\t1\t50\tchrQ\t+\t100\t150\tP#1\t0.05\t0.8\n"
                + "CNEE_MINUS\t100\t150\t2\t50\tchrQ\t-\t100\t150\tM#1\t0.05\t0.8\n",
                encoding="utf-8",
            )

            result = classify_orthology(
                stats,
                bed,
                tmp_path / "unused_index.json",
                orthology,
                query_bed,
                min_sequence_identity=0.7,
                max_fraction_gap=0.1,
                synteny_check=False,
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
                    ]
                )
                + "\n"
                + "CNEE_DUP\t0\t50\t1\t50\t1.0\t0\tchrQ\t+\t100\t150\tD1#1\t0.05\t0.8\n"
                + "CNEE_DUP__line2\t100\t150\t2\t50\t1.0\t0\tchrQ\t+\t500\t550\tD2#1\t0.05\t0.8\n",
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
                    min_sequence_identity=0.7,
                    max_fraction_gap=0.1,
                    synteny_check=False,
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

            classify_orthology(stats, bed, index, orthology, query_bed)

            result = write_synteny_verification(
                index,
                stats,
                query_bed,
                orthology,
                bed,
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
        # _write_inputs puts every CNEE on its own chromosome, so none has a
        # neighboring CNEE.
        for row in rows:
            self.assertEqual(row["distance_to_nearest_cnee_neighbor"], ".")
        # CNEE_GAPPED sits entirely in chain 3's gap between blocks
        # (100-150) and (250-300); it is 25 bases from each block edge.
        lost_rows = [row for row in rows if row["orthology_status"] == "lost"]
        self.assertTrue(lost_rows)
        for row in lost_rows:
            self.assertEqual(row["distance_to_nearest_aligned_block"], "25")
        cnee_uids = {row["cnee_uid"] for row in rows}
        for cnee_uid in cnee_uids:
            sizes_for_uid = {
                int(row["flank_size"]) for row in rows if row["cnee_uid"] == cnee_uid
            }
            self.assertEqual(sizes_for_uid, {500, 50})

    def test_distance_to_nearest_block_handles_overlap_and_gaps(self) -> None:
        from toga4cnee.models import ChainBlock
        from toga4cnee.synteny import BlockIndex, _distance_to_nearest_block

        blocks = BlockIndex.from_blocks(
            [
                ChainBlock(ref_start=100, ref_end=150, query_start=0, query_end=50),
                ChainBlock(ref_start=250, ref_end=300, query_start=50, query_end=100),
            ]
        )

        self.assertEqual(_distance_to_nearest_block(blocks, 120, 140), 0)
        self.assertEqual(_distance_to_nearest_block(blocks, 175, 225), 25)
        self.assertEqual(_distance_to_nearest_block(blocks, 180, 200), 30)
        self.assertEqual(_distance_to_nearest_block(blocks, 0, 50), 50)
        self.assertEqual(_distance_to_nearest_block(blocks, 350, 400), 50)

    def test_nearest_cnee_neighbor_distances_computes_gaps_on_same_chromosome(
        self,
    ) -> None:
        from toga4cnee.synteny import _nearest_cnee_neighbor_distances

        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "cnees.bed"
            bed.write_text(
                "chr1\t0\t50\tA\n"
                "chr1\t100\t150\tB\n"
                "chr1\t500\t550\tC\n"
                "chr2\t0\t10\tD\n",
                encoding="utf-8",
            )

            distances = _nearest_cnee_neighbor_distances(bed)

        self.assertEqual(distances["A"], 50)
        self.assertEqual(distances["B"], 50)
        self.assertEqual(distances["C"], 350)
        self.assertIsNone(distances["D"])

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

    def test_ref_neighbor_index_directional_lookup_and_chromosome_edges(self) -> None:
        from toga4cnee.synteny import RefNeighborIndex

        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "cnees.bed"
            bed.write_text(
                "chr1\t0\t50\tA\n"
                "chr1\t100\t150\tB\n"
                "chr1\t500\t550\tC\n"
                "chr2\t0\t10\tD\n",
                encoding="utf-8",
            )
            index = RefNeighborIndex.from_bed(bed)

        self.assertEqual(index.neighbor("B", "left"), ("A", 50))
        self.assertEqual(index.neighbor("B", "right"), ("C", 350))
        self.assertEqual(index.neighbor("A", "left"), (None, None))
        self.assertEqual(index.neighbor("C", "right"), (None, None))
        self.assertEqual(index.neighbor("D", "left"), (None, None))
        self.assertEqual(index.neighbor("D", "right"), (None, None))

    def test_query_neighbor_index_directional_lookup(self) -> None:
        from toga4cnee.synteny import QueryBedRecord, QueryNeighborIndex

        by_chrom = {
            "chrQ": [
                QueryBedRecord(chrom="chrQ", start=0, end=50, query_id="qA"),
                QueryBedRecord(chrom="chrQ", start=200, end=250, query_id="qB"),
            ]
        }
        index = QueryNeighborIndex.from_records(by_chrom)

        self.assertEqual(index.neighbor("chrQ", 100, "before"), ("qA", 50))
        self.assertEqual(index.neighbor("chrQ", 100, "after"), ("qB", 100))
        self.assertEqual(index.neighbor("chrQ", 10, "before"), (None, None))
        self.assertEqual(index.neighbor("chrQ2", 10, "before"), (None, None))

    def test_orthology_index_assigns_category_to_every_partner_ref(self) -> None:
        from toga4cnee.synteny import OrthologyIndex

        with tempfile.TemporaryDirectory() as tmpdir:
            orthology = Path(tmpdir) / "orthology.tsv"
            orthology.write_text(
                "query_cnee_id\tpartner_ref_cnees\torthology_category\n"
                "q_many\tCNEE_X,CNEE_Y\tmany2one\n"
                "NA\tCNEE_LOST\tlost\n",
                encoding="utf-8",
            )
            index = OrthologyIndex.from_path(orthology)

        self.assertEqual(index.category_by_ref_uid["CNEE_X"], "many2one")
        self.assertEqual(index.category_by_ref_uid["CNEE_Y"], "many2one")
        self.assertEqual(index.category_by_ref_uid["CNEE_LOST"], "lost")
        self.assertEqual(index.category_by_query_id["q_many"], "many2one")
        self.assertEqual(
            index.ref_partners_by_query_id["q_many"], frozenset({"CNEE_X", "CNEE_Y"})
        )
        self.assertNotIn("NA", index.category_by_query_id)

    def test_query_flank_anchors_own_interval_forward_and_reverse_strand(self) -> None:
        from toga4cnee.models import ChainRecord
        from toga4cnee.synteny import BlockIndex, _query_flank_anchors

        blocks = BlockIndex.from_blocks([])
        row = {
            "query_aligned_start": "100",
            "query_aligned_end": "200",
            "ref_start": "0",
            "ref_end": "50",
        }

        def _record(query_strand: str) -> ChainRecord:
            return ChainRecord(
                chain_id="1",
                score=0,
                ref_chrom="chr1",
                ref_size=1000,
                ref_strand="+",
                ref_start=0,
                ref_end=50,
                query_chrom="chrQ",
                query_size=1000,
                query_strand=query_strand,
                query_start=100,
                query_end=200,
                start_byte=0,
                byte_length=0,
            )

        forward = _query_flank_anchors(row, _record("+"), blocks)
        self.assertEqual(forward["left"], (100, "before"))
        self.assertEqual(forward["right"], (200, "after"))

        reverse = _query_flank_anchors(row, _record("-"), blocks)
        self.assertEqual(reverse["left"], (200, "after"))
        self.assertEqual(reverse["right"], (100, "before"))

    def test_query_flank_anchors_interpolates_from_flanking_blocks_for_lost_candidates(
        self,
    ) -> None:
        from toga4cnee.models import ChainBlock, ChainRecord
        from toga4cnee.synteny import BlockIndex, _query_flank_anchors

        blocks = BlockIndex.from_blocks(
            [
                ChainBlock(ref_start=0, ref_end=100, query_start=1000, query_end=1100),
                ChainBlock(ref_start=300, ref_end=400, query_start=1200, query_end=1300),
            ]
        )
        row = {
            "query_aligned_start": ".",
            "query_aligned_end": ".",
            "ref_start": "200",
            "ref_end": "250",
        }

        def _record(query_strand: str) -> ChainRecord:
            return ChainRecord(
                chain_id="1",
                score=0,
                ref_chrom="chr1",
                ref_size=1000,
                ref_strand="+",
                ref_start=0,
                ref_end=400,
                query_chrom="chrQ",
                query_size=1000,
                query_strand=query_strand,
                query_start=1000,
                query_end=1300,
                start_byte=0,
                byte_length=0,
            )

        forward = _query_flank_anchors(row, _record("+"), blocks)
        self.assertEqual(forward["left"], (1100, "before"))
        self.assertEqual(forward["right"], (1200, "after"))

        reverse = _query_flank_anchors(row, _record("-"), blocks)
        self.assertEqual(reverse["left"], (1000, "after"))
        self.assertEqual(reverse["right"], (1300, "before"))

        # A gap at the very edge of the chain's own aligned span: no
        # preceding block on the left, so that side is undeterminable
        # regardless of strand.
        edge_blocks = BlockIndex.from_blocks(
            [ChainBlock(ref_start=300, ref_end=400, query_start=1200, query_end=1300)]
        )
        edge_row = {
            "query_aligned_start": ".",
            "query_aligned_end": ".",
            "ref_start": "0",
            "ref_end": "50",
        }
        edge_forward = _query_flank_anchors(edge_row, _record("+"), edge_blocks)
        self.assertEqual(edge_forward["left"], (None, None))
        self.assertEqual(edge_forward["right"], (1200, "after"))
        edge_reverse = _query_flank_anchors(edge_row, _record("-"), edge_blocks)
        self.assertEqual(edge_reverse["left"], (None, None))
        self.assertEqual(edge_reverse["right"], (1300, "before"))

    def _write_lost_candidate_fixture(self, tmp_path: Path) -> tuple[Path, Path]:
        # A chain with a single gap (ref 100-300 <-> query 1100-1200) and three
        # reference CNEEs: CNEE_A and CNEE_C flank the gap and align solidly to
        # either side of it; CNEE_B sits entirely inside the gap ("lost").
        chain_text = "\n".join(
            [
                "chain 100 chr1 1000 + 0 400 chrQ 1000 + 1000 1300 10",
                "100 200 100",
                "100",
                "",
            ]
        )
        chain = tmp_path / "gap.chain"
        bed = tmp_path / "gap_cnees.bed"
        chain.write_text(chain_text, encoding="utf-8")
        bed.write_text(
            "chr1\t0\t50\tCNEE_A\n"
            "chr1\t200\t250\tCNEE_B\n"
            "chr1\t350\t400\tCNEE_C\n",
            encoding="utf-8",
        )
        return chain, bed

    def test_synteny_neighbor_verification_matches_true_when_query_neighbors_are_ref_neighbors_orthologs(
        self,
    ) -> None:
        if not _has_networkx():
            self.skipTest("networkx is not installed")
        from toga4cnee.orthology import classify_orthology
        from toga4cnee.synteny import write_synteny_neighbor_verification

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = self._write_lost_candidate_fixture(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"
            query_bed = tmp_path / "query_cnees.bed"
            orthology = tmp_path / "orthology.tsv"
            synteny_neighbors = tmp_path / "synteny_neighbors.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            write_overlap_stats(index, candidates, stats)
            classify_orthology(stats, bed, index, orthology, query_bed)

            result = write_synteny_neighbor_verification(
                index, stats, query_bed, orthology, bed, synteny_neighbors
            )
            rows = _read_tsv(synteny_neighbors)

        self.assertEqual(result.rows, len(rows))
        by_flank = {
            (row["cnee_uid"], row["flank"]): row
            for row in rows
            if row["orthology_status"] == "lost"
        }
        left = by_flank[("CNEE_B", "left")]
        right = by_flank[("CNEE_B", "right")]

        self.assertEqual(left["nearest_ref_cnee_id"], "CNEE_A")
        self.assertEqual(left["distance_to_ref_cnee"], "150")
        self.assertEqual(left["nearest_ref_cnee_status"], "one2one")
        self.assertEqual(left["nearest_query_cnee_id"], "q_CNEE_A")
        self.assertEqual(left["distance_to_query_cnee"], "50")
        self.assertEqual(left["nearest_query_cnee_status"], "one2one")
        self.assertEqual(left["synteny_matches"], "TRUE")

        self.assertEqual(right["nearest_ref_cnee_id"], "CNEE_C")
        self.assertEqual(right["distance_to_ref_cnee"], "100")
        self.assertEqual(right["nearest_query_cnee_id"], "q_CNEE_C")
        self.assertEqual(right["distance_to_query_cnee"], "50")
        self.assertEqual(right["synteny_matches"], "TRUE")

        # CNEE_B is "lost" -- no query projection exists at all.
        self.assertEqual(left["ref_coordinates"], "chr1:200-250")
        self.assertEqual(left["query_projection_id"], "NA")
        self.assertEqual(left["query_projection_coordinates"], "NA")

        # CNEE_A and CNEE_C are themselves solid one2one candidates (positive
        # controls), each with a real graph-assigned query projection.
        by_status = {row["cnee_uid"]: row for row in rows if row["orthology_status"] == "one2one"}
        self.assertEqual(by_status["CNEE_A"]["ref_coordinates"], "chr1:0-50")
        self.assertEqual(by_status["CNEE_A"]["query_projection_id"], "q_CNEE_A")
        self.assertEqual(
            by_status["CNEE_A"]["query_projection_coordinates"], "chrQ:1000-1050"
        )
        self.assertEqual(by_status["CNEE_C"]["query_projection_id"], "q_CNEE_C")
        self.assertEqual(
            by_status["CNEE_C"]["query_projection_coordinates"], "chrQ:1250-1300"
        )

    def test_candidate_query_projection_by_status(self) -> None:
        from toga4cnee.synteny import SyntenyTarget, _candidate_query_projection

        lost_target = SyntenyTarget(
            row={
                "query_aligned_start": ".",
                "query_aligned_end": ".",
                "query_chrom": "chrQ",
                "query_projection_name": "CNEE_LOST#1",
            },
            status="lost",
            output_id="CNEE_LOST",
        )
        self.assertEqual(_candidate_query_projection(lost_target), ("NA", "NA"))

        uncertain_target = SyntenyTarget(
            row={
                "query_aligned_start": "100",
                "query_aligned_end": "150",
                "query_chrom": "chrQ",
                "query_projection_name": "CNEE_UNC#2",
            },
            status="uncertain_loss",
            output_id="CNEE_UNC",
        )
        self.assertEqual(
            _candidate_query_projection(uncertain_target),
            ("CNEE_UNC#2", "chrQ:100-150"),
        )

        weak_target = SyntenyTarget(
            row={
                "query_aligned_start": "200",
                "query_aligned_end": "250",
                "query_chrom": "chrQ",
                "query_projection_name": "CNEE_WEAK#3",
            },
            status="weak",
            output_id="q_CNEE_WEAK_1",
        )
        self.assertEqual(
            _candidate_query_projection(weak_target),
            ("q_CNEE_WEAK_1", "chrQ:200-250"),
        )

        one2one_target = SyntenyTarget(
            row={
                "query_aligned_start": "300",
                "query_aligned_end": "350",
                "query_chrom": "chrQ",
                "query_projection_name": "CNEE_ONE#4",
            },
            status="one2one",
            output_id="q_CNEE_ONE",
        )
        self.assertEqual(
            _candidate_query_projection(one2one_target),
            ("q_CNEE_ONE", "chrQ:300-350"),
        )

    def test_synteny_neighbor_verification_flags_mismatch_for_rearranged_query_neighbor(
        self,
    ) -> None:
        from toga4cnee.synteny import write_synteny_neighbor_verification

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = self._write_lost_candidate_fixture(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"
            query_bed = tmp_path / "query_cnees.bed"
            orthology = tmp_path / "orthology.tsv"
            synteny_neighbors = tmp_path / "synteny_neighbors.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            write_overlap_stats(index, candidates, stats)

            # Hand-author query BED/orthology so the query CNEE physically
            # nearest the left side of CNEE_B's projected gap (q_WRONG, right
            # where CNEE_A's ortholog would sit) is *not* actually derived
            # from CNEE_A -- simulating a rearrangement. The right side stays
            # a genuine match.
            query_bed.write_text(
                "chrQ\t1000\t1050\tq_WRONG\t0\t+\n"
                "chrQ\t1250\t1300\tq_CNEE_C\t0\t+\n",
                encoding="utf-8",
            )
            orthology.write_text(
                "query_cnee_id\tpartner_ref_cnees\torthology_category\n"
                "q_WRONG\tCNEE_OTHER\tone2one\n"
                "q_CNEE_C\tCNEE_C\tone2one\n",
                encoding="utf-8",
            )

            result = write_synteny_neighbor_verification(
                index, stats, query_bed, orthology, bed, synteny_neighbors
            )
            rows = _read_tsv(synteny_neighbors)

        self.assertEqual(result.rows, len(rows))
        by_flank = {
            (row["cnee_uid"], row["flank"]): row
            for row in rows
            if row["orthology_status"] == "lost"
        }
        left = by_flank[("CNEE_B", "left")]
        right = by_flank[("CNEE_B", "right")]

        self.assertEqual(left["nearest_ref_cnee_id"], "CNEE_A")
        self.assertEqual(left["nearest_query_cnee_id"], "q_WRONG")
        self.assertEqual(left["synteny_matches"], "FALSE")

        self.assertEqual(right["nearest_ref_cnee_id"], "CNEE_C")
        self.assertEqual(right["nearest_query_cnee_id"], "q_CNEE_C")
        self.assertEqual(right["synteny_matches"], "TRUE")

    def test_synteny_neighbor_verification_reports_dot_when_undeterminable(self) -> None:
        from toga4cnee.synteny import write_synteny_neighbor_verification

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            chain, bed = self._write_lost_candidate_fixture(tmp_path)
            index = tmp_path / "index.json"
            candidates = tmp_path / "candidates.tsv"
            stats = tmp_path / "overlap_stats.tsv"
            query_bed = tmp_path / "query_cnees.bed"
            orthology = tmp_path / "orthology.tsv"
            synteny_neighbors = tmp_path / "synteny_neighbors.tsv"

            build_index(chain, index)
            write_candidates(index, bed, candidates)
            write_overlap_stats(index, candidates, stats)
            # No query BED / orthology rows at all -- every neighbor lookup
            # comes up empty, so nothing can be resolved.
            query_bed.write_text("", encoding="utf-8")
            orthology.write_text(
                "query_cnee_id\tpartner_ref_cnees\torthology_category\n",
                encoding="utf-8",
            )

            result = write_synteny_neighbor_verification(
                index, stats, query_bed, orthology, bed, synteny_neighbors
            )
            rows = _read_tsv(synteny_neighbors)

        self.assertEqual(result.rows, len(rows))
        lost_rows = [row for row in rows if row["orthology_status"] == "lost"]
        self.assertTrue(lost_rows)
        for row in lost_rows:
            self.assertEqual(row["nearest_query_cnee_id"], ".")
            self.assertEqual(row["distance_to_query_cnee"], ".")
            self.assertEqual(row["nearest_query_cnee_status"], ".")
            self.assertEqual(row["synteny_matches"], ".")
            # Ref-side neighbors are independent of query data and still resolve.
            self.assertIn(row["nearest_ref_cnee_id"], {"CNEE_A", "CNEE_C"})

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
