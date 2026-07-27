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
`flank_size`, `flank_seq_identity`, `flank_fraction_gap`,
`distance_to_nearest_aligned_block`, and `distance_to_nearest_cnee_neighbor`.
The last two are properties of the candidate itself (not the flank check), so
they repeat identically across a candidate's rows: distance to the nearest
aligned block in the same chain (`0` if the candidate already overlaps one --
mainly informative for `lost` candidates, which by definition don't), and
distance to the nearest neighboring reference CNEE annotation on the same
chromosome (`.` if it's the only CNEE on that chromosome). Use
`--synteny-flank-sizes` (one or more values, default `10000 1000 100`) to set
the reference bases upstream/downstream checked at each scale, and
`--synteny-max-one2one-controls` to control the maximum number of one2one
candidates, default `50000`.

A second, complementary synteny check is written to
`overlap_stats.synteny_neighbors.tsv`. Where the table above compares raw
flanking sequence, this one is purely positional: for the same candidates
(weak, lost, uncertain_loss, one2one controls), it asks whether the query
region corresponding to the candidate is still flanked by the *same*
reference CNEEs that flank it in the reference -- the kind of check that
catches a rearrangement a sequence-identity scan alone wouldn't. It is
long-format with one row per candidate per flank side (`left`/`right`),
with columns `cnee_uid`, `orthology_status`, `query_projection_id`,
`ref_coordinates`, `query_projection_coordinates`, `flank`,
`nearest_ref_cnee_id`, `distance_to_ref_cnee`, `nearest_ref_cnee_status`,
`nearest_query_cnee_id`, `distance_to_query_cnee`,
`nearest_query_cnee_status`, and `synteny_matches`.

`query_projection_id`, `ref_coordinates`, and `query_projection_coordinates`
are properties of the candidate itself (not the flank check), so they repeat
identically across a candidate's two rows. `ref_coordinates` is always the
candidate's reference CNEE interval, as `chr:start-end`. `query_projection_id`
and `query_projection_coordinates` describe the candidate's own query-side
projection: for `weak`/`one2one` candidates, the graph-assigned query CNEE id
and its `chr:start-end` coordinates; for `uncertain_loss` candidates (aligned
but never solid enough to be graph-assigned a query CNEE), the raw
per-alignment projection name and coordinates instead; for `lost` candidates
(no aligned query bases at all), both are `NA`.

For each side, `nearest_ref_cnee_id`/`distance_to_ref_cnee` report the
nearest neighboring reference CNEE annotation (independent of any chain),
while `nearest_query_cnee_id`/`distance_to_query_cnee` report the nearest
query CNEE to the corresponding projected position -- the candidate's own
aligned query interval for weak/uncertain_loss/one2one rows, or (for `lost`
candidates) a position interpolated from the chain blocks immediately
flanking the candidate's gap. `synteny_matches` is `TRUE`/`FALSE` when the
query neighbor is/isn't one of the reference neighbor's orthology partners,
or `.` when either neighbor can't be resolved (chromosome/chain edge, or no
query CNEE nearby). Use `--out-synteny-neighbors` to change its path (default
derived from `--out-stats`); it reuses `--synteny-max-one2one-controls`,
`--min-sequence-identity`, and `--max-fraction-gap` and needs no ref/query
FASTA, since it does no sequence comparison.

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
