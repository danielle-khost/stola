# toga4cnee
Development repo for toga-like annotations of CNEEs

## Current CLI

Run the first workflow pass with one command:

```bash
toga4cnee alignment.chain.gz reference_cnees.bed \
  --out-bed query_cnees.bed \
  --out-stats overlap_stats.tsv
```

The command indexes the chain file, finds candidate chain/CNEE pairs, computes
overlap statistics for every candidate, writes BED12 query projections for
candidate pairs with aligned CNEE bases, writes a concise summary at
`overlap_stats.summary.tsv`, and writes two SVG QC plots by default:
`overlap_stats.num_query_alignments.svg` and
`overlap_stats.overlap_fraction_density.svg`. Use `--work-dir` to keep
intermediate index and candidate TSV files.

Create the summary from an existing stats file:

```bash
toga4cnee summarize overlap_stats.tsv --cnee-bed reference_cnees.bed --out cnee_summary.tsv
```

Create the QC plots from existing outputs:

```bash
toga4cnee plot cnee_summary.tsv overlap_stats.tsv \
  --out-alignment-histogram num_query_alignments.svg \
  --out-overlap-density overlap_fraction_density.svg
```
