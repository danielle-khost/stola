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
concise summary at `overlap_stats.summary.tsv`. Use `--work-dir` to keep
intermediate index and candidate TSV files. Diagnostic figures are generated
separately from completed outputs.

The current orthology pass draws an edge for every aligned query projection,
stores `sequence_identity` and `fraction_gap` on that edge, then classifies the
connected components formed by solid edges. By default, a solid edge has
`fraction_gap < 0.9`; when FASTAs are supplied it must also have
`sequence_identity > 0.7`. Tune these with `--max-fraction-gap` and
`--min-sequence-identity`. Categories are `one2one`, `one2many`, `many2one`,
`lost`, `missing`, `uncertain_loss`, or the catch-all `complex/many2many`. The
orthology TSV has three columns: `query_cnee_id`, `partner_ref_cnees`, and
`orthology_category`; non-query reference outcomes are written with `NA` as the
query ID.

When both `--ref-fasta` and `--query-fasta` are supplied, the stats and summary
also include direct chain-aligned sequence identity and mismatch counts. These
metrics require Biopython.

Create the summary from an existing stats file:

```bash
toga4cnee summarize overlap_stats.tsv --cnee-bed reference_cnees.bed --out cnee_summary.tsv
```

Create the original QC plots from existing outputs:

```bash
toga4cnee plot cnee_summary.tsv overlap_stats.tsv \
  --out-alignment-histogram num_query_alignments.svg \
  --out-overlap-density overlap_fraction_density.svg \
  --out-sequence-identity-density sequence_identity_density.svg
```

For test/evaluation runs, generate the full diagnostic figure set, including
orthology category counts and category-grouped sequence identity/gap fraction
densities:

```bash
python tests/plot_toga4cnee_outputs.py \
  --summary overlap_stats.summary.tsv \
  --stats overlap_stats.tsv \
  --orthology overlap_stats.orthology.tsv \
  --out-dir figures
```
