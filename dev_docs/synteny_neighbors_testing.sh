source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biopython

PYTHONPATH=src python - <<'PY'
import time
from toga4cnee.synteny import write_synteny_neighbor_verification

t0 = time.time()
result = write_synteny_neighbor_verification(
    index_path="tests/test_runs/pigeonChickenStrict/work/chain_index.json",
    stats_path="tests/test_runs/pigeonChickenStrict/overlap_stats.tsv",
    query_bed_path="tests/test_runs/pigeonChickenStrict/query_cnees.bed",
    orthology_path="tests/test_runs/pigeonChickenStrict/overlap_stats.orthology.tsv",
    cnee_bed_path="tests/gGal7_cnee_final_data_filtered.fixed.bed",
    out_path="tests/test_runs/pigeonChickenStrict/overlap_stats.synteny_neighbors.tsv",
    min_sequence_identity=0.7,
    max_fraction_gap=0.1,
)
print(f"Wrote {result.rows} synteny neighbor-check rows in {time.time() - t0:.1f}s")
PY
