import sys
import tempfile
import unittest
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tunic.bed import bed_missing_name_column, parse_bed, write_named_bed
from tunic.chain import iter_chain_blocks, parse_chain_header
from tunic.intervals import overlap_length


class ParserTests(unittest.TestCase):
    def test_parse_chain_header_normalizes_query_minus_strand(self) -> None:
        header = "chain 200 chr1 1000 + 110 240 chrQ 1000 - 100 230 2"

        record = parse_chain_header(header, start_byte=17, byte_length=42)

        self.assertEqual(record.chain_id, "2")
        self.assertEqual(record.score, 200)
        self.assertEqual(record.ref_chrom, "chr1")
        self.assertEqual(record.ref_start, 110)
        self.assertEqual(record.ref_end, 240)
        self.assertEqual(record.query_strand, "-")
        self.assertEqual(record.query_start, 770)
        self.assertEqual(record.query_end, 900)
        self.assertEqual(record.start_byte, 17)
        self.assertEqual(record.byte_length, 42)

    def test_iter_chain_blocks_handles_plus_and_minus_query_strands(self) -> None:
        plus_chain = "\n".join(
            [
                "chain 100 chr1 1000 + 100 300 chrQ 1000 + 500 650 1",
                "50 50 0",
                "100",
                "",
            ]
        )
        minus_chain = "\n".join(
            [
                "chain 200 chr1 1000 + 110 240 chrQ 1000 - 100 230 2",
                "130",
                "",
            ]
        )

        plus_blocks = iter_chain_blocks(plus_chain)
        minus_blocks = iter_chain_blocks(minus_chain)

        self.assertEqual(plus_blocks[0].ref_start, 100)
        self.assertEqual(plus_blocks[0].ref_end, 150)
        self.assertEqual(plus_blocks[0].query_start, 500)
        self.assertEqual(plus_blocks[0].query_end, 550)
        self.assertEqual(plus_blocks[1].ref_start, 200)
        self.assertEqual(plus_blocks[1].ref_end, 300)
        self.assertEqual(plus_blocks[1].query_start, 550)
        self.assertEqual(plus_blocks[1].query_end, 650)
        self.assertEqual(minus_blocks[0].ref_start, 110)
        self.assertEqual(minus_blocks[0].ref_end, 240)
        self.assertEqual(minus_blocks[0].query_start, 770)
        self.assertEqual(minus_blocks[0].query_end, 900)

    def test_parse_bed_accepts_bed3_and_bed4_with_stable_unique_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "cnees.bed"
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t10\t20",
                        "chr1\t30\t40\tCNEE_A",
                        "chr1\t50\t60\tCNEE_A",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            intervals = parse_bed(bed)

        self.assertEqual(intervals[0].cnee_id, "cnee_000001_chr1_10_20")
        self.assertEqual(intervals[0].cnee_uid, "cnee_000001_chr1_10_20")
        self.assertEqual(intervals[1].cnee_id, "CNEE_A")
        self.assertEqual(intervals[1].cnee_uid, "CNEE_A")
        self.assertEqual(intervals[2].cnee_id, "CNEE_A")
        self.assertEqual(intervals[2].cnee_uid, "CNEE_A__line3")

    def test_parse_bed_warns_on_zero_length_intervals(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "cnees.bed"
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t10\t20\tCNEE_NORMAL",
                        "chr1\t30\t30\tCNEE_ZERO_LEN",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                intervals = parse_bed(bed)

        self.assertEqual(len(intervals), 2)
        self.assertTrue(
            any("zero-length" in str(w.message) for w in caught),
            "expected a zero-length CNEE warning",
        )

    def test_bed_missing_name_column_detects_bed3_and_blank_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            named_bed = Path(tmpdir) / "named.bed"
            named_bed.write_text("chr1\t10\t20\tCNEE_A\n", encoding="utf-8")
            unnamed_bed = Path(tmpdir) / "unnamed.bed"
            unnamed_bed.write_text("chr1\t10\t20\n", encoding="utf-8")
            blank_name_bed = Path(tmpdir) / "blank.bed"
            blank_name_bed.write_text("chr1\t10\t20\t\n", encoding="utf-8")

            self.assertFalse(bed_missing_name_column(named_bed))
            self.assertTrue(bed_missing_name_column(unnamed_bed))
            self.assertTrue(bed_missing_name_column(blank_name_bed))

    def test_write_named_bed_uses_cnee_uid_in_file_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "cnees.bed"
            bed.write_text(
                "\n".join(
                    [
                        "chr1\t10\t20",
                        "chr1\t30\t40\tCNEE_A",
                        "chr1\t50\t60\tCNEE_A",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            intervals = parse_bed(bed)
            out_path = Path(tmpdir) / "cnees_named.bed"

            count = write_named_bed(intervals, out_path)

            self.assertEqual(count, 3)
            self.assertEqual(
                out_path.read_text(encoding="utf-8").splitlines(),
                [
                    "chr1\t10\t20\tcnee_000001_chr1_10_20",
                    "chr1\t30\t40\tCNEE_A",
                    "chr1\t50\t60\tCNEE_A__line3",
                ],
            )

    def test_overlap_length_uses_half_open_intervals(self) -> None:
        self.assertEqual(overlap_length(10, 20, 20, 30), 0)
        self.assertEqual(overlap_length(10, 25, 20, 30), 5)
        self.assertEqual(overlap_length(10, 30, 12, 18), 6)


if __name__ == "__main__":
    unittest.main()
