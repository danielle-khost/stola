source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biopython

CHAIN="tests/HLgalGal7.HLmorBass1B.allfilled.chain"
REFBED="tests/gGal7_cnee_final_data_filtered.fixed.bed"
OUTPATH="tests/test_runs/morBass"
mkdir -p $OUTPATH
REFFASTA="tests/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"
QUERYFASTA="tests/GCA_031468815.1_bMorBas2.hap2_genomic.fna"
PREFIX="morBass"

PYTHONPATH=src python -m toga4cnee.cli run \
  $CHAIN \
  $REFBED \
  --out-bed ${OUTPATH}/query_cnees.bed \
  --out-stats ${OUTPATH}/overlap_stats.tsv \
  --out-summary ${OUTPATH}/overlap_stats.summary.tsv \
  --out-orthology ${OUTPATH}/overlap_stats.orthology.tsv \
  --ref-fasta $REFFASTA \
  --query-fasta $QUERYFASTA \
  --max-fraction-gap 0.1 \
  --min-sequence-identity 0.7 \
  --work-dir ${OUTPATH}/work


python tests/plot_toga4cnee_outputs.py \
  --summary ${OUTPATH}/overlap_stats.summary.tsv \
  --stats ${OUTPATH}/overlap_stats.tsv \
  --orthology ${OUTPATH}/overlap_stats.orthology.tsv \
  --out-dir ${OUTPATH}/figures \
  --prefix $PREFIX \
  --max-fraction-gap 0.1 \
  --min-sequence-identity 0.7
