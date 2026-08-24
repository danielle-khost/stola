# toga4cnee
Development repo for toga-like annotations of CNEEs

## Installation

### Prerequisites 
*Dedicated conda environment TBD*
For now, the pipeline and its dependencies are installed via `pip`:

```
git clone https://github.com/danielle-khost/toga4cnee.git

cd toga4cnee
pip install .
```

## Current CLI

Run the first workflow pass with one command:

```bash
toga4cnee alignment.chain.gz reference_cnees.bed reference.fa/2bit query.fa/2bit \
  --outdir results/ --query_prefix q
```

The reference genome, BED and chain file must use exactly the same chromosome/contig
naming convention.

### Required input
- `.chain` alignment file of query species vs reference. Can be uncompressed or `.gz`
- BED annotation file of reference CNEEs. Format must be either:
  - `BED3`: standard 3 column BED with NO name column, or
  - `BED4+`: expanded BED with names in the *4th column*
- Reference + query genomes: can be in FASTA or UCSC 2bit format (can mix and match)
  - *Note: due to a limitation of Biopython, 2bit files for large (>4Gbase) genomes cannot be read; convert them to FASTA first instead*

Be sure that the naming conventions match between the chain, BED and genome files!

### Optional arguments
`--outdir`: specify where output files will be generated (default: current directory `./`)

`--query_prefix`: the prefix (e.g. species ID) to add to each query CNEE ID. Defaults to `q`, e.g. `q_ce1234`.


## Pipeline overview

*TBD*

### Orthology categories
- `one2one`: single reference CNEE with single query projection
- `one2many`: single reference CNEE duplicated in query (>= 2 projections)
- `many2one`: multiple reference CNEEs match single query CNEE
- `complex`: multiple reference CNEEs match multiple query CNEEs

*Query only:*
- `pruned`: query projection that failed filtering step

*Reference only:*
- `unaligned`: ref CNEE does not overlap **any** chain aligned blocks OR gaps
- `probable_loss`: ref CNEE occurs in chain gap that does NOT look like contig boundary
- `assembly_gap_loss`: ref CNEE occurs in chain gap that is likely a contig boundary
- `uncertain_loss`: ref CNEE with zero projections that passed the filtering steps

## Output
`query_cnees.bed`: annotation file of final CNEE projections to the query genome

`query_ortho_status.tsv`: key linking query projection IDs to reference IDs, with their orthology relationship

`ref_ortho_stats.tsv`: key listing each reference CNEE, its inital (pre-pruning) orthology category, number of projections before and after pruning, and final orthology category

`pruned_query_projections.tsv`: key of all query projections that failed filtering step (i.e. were removed)

`overlap_stats.tsv`: raw data table with overlap stats for each reference CNEE

`overlap_stats_summary.tsv`: summary of `overlap_stats.tsv` with lists for each ref CNEE of number of projections (i.e. overlaps), fraction overlapping aligned bases (and number of gaps), sequence identity vs query (and number of mismatches)


## Summary figures

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
  --orthology query_ortho_status.tsv \
  --out-dir figures
```
