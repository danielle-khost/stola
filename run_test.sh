source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biopython

CHAIN="tests/HLgalGal7.HLmorBass1B.allfilled.chain"
REFBED="tests/gGal7_cnee_final_data_filtered.fixed.bed"
OUTPATH="tests/test_runs/morBass"
mkdir -p $OUTPATH
REFFASTA="tests/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"
QUERYFASTA="tests/GCA_031468815.1_bMorBas2.hap2_genomic.fna"
PREFIX="morBass"

PYTHONPATH=src python -m stola.cli \
  $CHAIN \
  $REFBED \
  $REFFASTA \
  $QUERYFASTA \
  --outdir $OUTPATH \
  --max-fraction-gap 0.1 \
  --min-sequence-identity 0.7


python tests/plot_stola_outputs.py \
  --summary ${OUTPATH}/summary_figures/overlap_stats.summary.tsv \
  --stats ${OUTPATH}/overlap_stats.tsv \
  --orthology ${OUTPATH}/query_ortho_status.tsv \
  --out-dir ${OUTPATH}/figures \
  --prefix $PREFIX \
  --max-fraction-gap 0.1 \
  --min-sequence-identity 0.7
