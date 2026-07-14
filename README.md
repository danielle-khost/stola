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
overlap statistics for every candidate, writes an all-output BED/orthology TSV,
writes strict companion outputs, and writes concise summaries. The requested
`--out-bed`, `--out-stats`, `--out-summary`, and `--out-orthology` paths are
the all outputs. Strict outputs are written by default as `.strict` companions,
for example `query_cnees.strict.bed`, `overlap_stats.strict.tsv`,
`overlap_stats.strict.summary.tsv`, and `overlap_stats.strict.orthology.tsv`.
Use `--work-dir` to keep intermediate index and candidate TSV files.
Diagnostic figures are generated separately from completed outputs.

The reference BED and chain file must use exactly the same chromosome/contig
naming convention -- there is no accession-version alias fallback. A BED
chromosome that only fails to match because of a stripped or added accession
version (e.g. `chr1.1` vs `chr1`) raises an informative error rather than
silently being treated as having no chain coverage; rename one side to match
the other and rerun.

The current orthology pass draws an edge for every aligned query projection and
classifies graph components. Query-side projections are only merged into the
same query locus when they overlap on the same strand; overlapping projections
on opposite strands become separate loci. By default, a projection is solid
when it has
`fraction_gap < 0.1`; when FASTAs are supplied it must also have
`sequence_identity > 0.7`. Both criteria must be met by the same projection --
metrics are never blended across the multiple projections that can share one
ref/query edge, so an edge is solid only if at least one of its projections
passes both checks on its own. Tune these with `--max-fraction-gap` and
`--min-sequence-identity`. All output keeps every aligned query projection for
a reference CNEE that has at least one solid projection. Strict output keeps
only solid projections; references with weak aligned projections but no solid
projection are retained as `uncertain_loss` in both modes. Categories are
`one2one`, `one2many`, `many2one`, `lost`, `missing`, `uncertain_loss`, or the
catch-all `complex/many2many`. The orthology TSV has three columns:
`query_cnee_id`, `partner_ref_cnees`, and `orthology_category`; non-query
reference outcomes are written with `NA` as the query ID.

Reference CNEE identity is tracked internally by the unique `cnee_uid`, not the
display `cnee_id` (BED column 4), so two BED intervals that happen to share a
name are never merged into one reference locus; `partner_ref_cnees` and the
summary TSV report `cnee_uid` values, which match `cnee_id` unless the input
BED has duplicate names. Parsing a BED with duplicate names emits a warning,
as does parsing a BED with zero-length intervals (which can never overlap a
chain and always end up `missing`).

The workflow also writes an exploratory flanking-synteny verification table at
`overlap_stats.synteny.tsv`. This table reports weak all-mode projections,
`lost` and `uncertain_loss` candidates, and a capped sample of `one2one`
positive controls. It is long-format: each candidate gets one row per
requested flank size, with columns `cnee_uid`, `orthology_status`,
`flank_size`, `flank_seq_identity`, and `flank_fraction_gap`. Use
`--synteny-flank-sizes` (one or more values, default `10000 1000 100`) to set
the reference bases upstream/downstream checked at each scale, and
`--synteny-max-one2one-controls` to control the maximum number of one2one
candidates, default `50000`.

When both `--ref-fasta` and `--query-fasta` are supplied, the stats and summary
also include direct chain-aligned sequence identity and mismatch counts. These
metrics require Biopython. FASTA indexing requires uniform line width within
each record (every line except the last must be the same length); a ragged
FASTA raises a clear error naming the offending sequence rather than silently
computing wrong offsets -- re-wrap it to a fixed width (e.g.
`seqkit seq -w 60 in.fa > out.fa`) and rerun.

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
