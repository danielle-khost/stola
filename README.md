# toga4cnee
Development repo for toga-like annotations of CNEEs

## Current CLI

Run the first workflow pass with one command:

```bash
toga4cnee alignment.chain.gz reference_cnees.bed \
  --out-bed query_cnees.bed \
  --out-stats overlap_stats.tsv \
  --ref-fasta reference.fa \
  --query-fasta query.fa
```

The command indexes the chain file, finds candidate chain/CNEE pairs, computes
overlap statistics for every candidate, clusters overlapping query-side
projections into graph-assigned query CNEEs, writes those query CNEEs to BED,
writes graph-based orthology calls at `overlap_stats.orthology.tsv`, writes a
concise summary at `overlap_stats.summary.tsv`, and writes SVG QC plots by
default: `overlap_stats.num_query_alignments.svg`,
`overlap_stats.overlap_fraction_density.svg`, and
`overlap_stats.sequence_identity_density.svg`. Use `--work-dir` to keep
intermediate index and candidate TSV files.

The current orthology pass draws an edge for every aligned query projection,
then classifies connected components as `one2one`, `one2many`, `many2one`,
`lost`, `missing`, or the catch-all `complex/many2many`. Filtering logic is
intentionally deferred so test runs can expose the raw duplication/loss
patterns. The orthology TSV has three columns: `query_cnee_id`,
`partner_ref_cnees`, and `orthology_category`; lost and missing reference CNEEs
are written with `NA` as the query ID.

When both `--ref-fasta` and `--query-fasta` are supplied, the stats and summary
also include direct chain-aligned sequence identity and mismatch counts. These
metrics require Biopython.

Create the summary from an existing stats file:

```bash
toga4cnee summarize overlap_stats.tsv --cnee-bed reference_cnees.bed --out cnee_summary.tsv
```

Create the QC plots from existing outputs:

```bash
toga4cnee plot cnee_summary.tsv overlap_stats.tsv \
  --out-alignment-histogram num_query_alignments.svg \
  --out-overlap-density overlap_fraction_density.svg \
  --out-sequence-identity-density sequence_identity_density.svg
```
