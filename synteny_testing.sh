source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biopython

PYTHONPATH=src python - <<'PY'
from toga4cnee.sequence_compare import SequenceComparator
from toga4cnee.synteny import write_synteny_verification

work_dir = "tests/test_runs/pigeonChickenStrict/work"

comparator = SequenceComparator(
    "tests/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna",
    "tests/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna",
    work_dir,
)

try:
    result = write_synteny_verification(
        index_path="tests/test_runs/pigeonChickenStrict/work/chain_index.json",
        stats_path="tests/test_runs/pigeonChickenStrict/overlap_stats.tsv",
        query_bed_path="tests/test_runs/pigeonChickenStrict/query_cnees.bed",
        orthology_path="tests/test_runs/pigeonChickenStrict/overlap_stats.orthology.tsv",
        out_path="tests/test_runs/pigeonChickenStrict/overlap_stats.synteny.flank10000.tsv",
        flank_sizes=[10000],
        min_sequence_identity=0.7,
        max_fraction_gap=0.1,
        sequence_comparator=comparator,
    )
    print(f"Wrote {result.rows} synteny rows")
finally:
    comparator.close()
PY
