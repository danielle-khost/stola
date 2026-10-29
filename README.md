# `stola`
[![DOI](https://zenodo.org/badge/1267684507.svg)](https://doi.org/10.5281/zenodo.23779867)

Inspired by [TOGA2](https://github.com/hillerlab/TOGA2), `stola` uses a pairwise chain alignment between a query and reference to annotate conserved non-exonic elements (CNEEs) in the query as well as infer orthological relationships between query and reference for each CNEE.

## Installation

### Conda
stola is packaged for conda (a bioconda release is planned). Until then, build the package from a checkout and install it into its own environment:

```
git clone https://github.com/harvardinformatics/stola.git
cd stola

conda create -n conda-build -c conda-forge conda-build
conda run -n conda-build conda build recipe -c conda-forge -c bioconda --output-folder dist/conda

conda create -n stola -c ./dist/conda -c conda-forge -c bioconda stola
conda activate stola
```

This installs the `stola` pipeline and the `stola-compare-cross-orthology` utility, with all of their dependencies.

### pip
You can also install from a checkout with `pip` instead:

```
cd stola/
pip install .
```

## Pipeline overview
CNEE annotations in the reference are overlapped with chain alignments to identify candidate CNEE projections in the query. If a reference CNEE overlaps only a single query locus, it is accepted as a one-to-one ortholog (one2one); for multi-projection cases, candidates are filtered to remove spurious/repeat-driven alignments as follows:

### Multi-projection filtering
Projections of ref CNEEs with more than one candidate (i.e. putative one2many, many2one, complex) are evaluated based on their alignment over the CNEE + flanking regions vs competing alignment at that locus. For each query locus, every chain aligned gets a local score (sequence identity over the locus x aligned fraction over the locus plus flanks). The local score of each projection is normalized per-locus by dividing by the max-scoring projection (R = S / S(max)) and a projection is kept if its relative score is above a fraction of the best chain at that locus (R > 0.5).

### Orthology categories
- `one2one`: single reference CNEE with single query projection
- `one2many`: single reference CNEE duplicated in query (>= 2 projections)
- `many2one`: single query CNEE duplicated in reference
- `complex/many2many`: multiple reference CNEEs match multiple query CNEEs

*Reference only:*
- `unaligned`: ref CNEE does not overlap **any** chain aligned blocks OR gaps
- `probable_loss`: ref CNEE occurs in chain gap that does NOT look like contig boundary
- `assembly_gap_loss`: ref CNEE occurs in chain gap that is likely a contig boundary
- `uncertain_loss`: ref CNEE with zero projections that passed the filtering steps


## Useage

Run the workflow with one command:

```bash
stola alignment.chain.gz reference_cnees.bed reference.fa/2bit query.fa/2bit \
  --outdir [path] --query-prefix [string]
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

`--query-prefix`: the prefix (e.g. species ID) to add to each query CNEE ID. Defaults to `q`, e.g. `q_ce1234`.

`--relative-score-threshold`: minimum fraction of the locus's best local score needed to keep a projection (default `0.5`)

`--soft-depth-limit`: number of chains scored per locus, taking the best by chain score; projections on any other chain are pruned as repeat-driven pileups (default `20`)

`--flank-window`: flank size (bp, each side) included in the aligned fraction (default `10000`)

`--flank-any-coverage-threshold`: at loci where only one chain aligns, the minimum fraction of the projection's flanks that must be aligned by any chain (default `0.5`)

`--no-filter-multi-projections`: skip this filtering step entirely


## Output
`query_cnees.bed`: annotation file of final CNEE projections to the query genome

`query_ortho_status.tsv`: key linking query projection IDs to reference IDs, with their orthology relationship

`ref_ortho_stats.tsv`: key listing each reference CNEE, its inital (pre-pruning) orthology category, number of projections before and after pruning, and final orthology category

`pruned_query_projections.tsv`: key of all query projections that failed filtering step (i.e. were removed)

`overlap_stats.tsv`: raw data table with overlap stats for each reference CNEE

`summary_figures/`: summary tables and plots:
- `overlap_stats_summary.tsv`: summary of `overlap_stats.tsv` with lists for each ref CNEE of number of projections (i.e. overlaps), fraction overlapping aligned bases (and number of gaps), sequence identity vs query (and number of mismatches)
- `final_category_counts.svg/.tsv`: reference CNEEs per final orthology category
- `initial_to_final_reassignment.svg/.tsv`: how filtering recategorized multi-projection reference CNEEs
- `rejection_reason_counts.svg/.tsv`: pruned projections by reason (`depth`, `no_flank_support`, `low_score`)
- `num_query_alignments.svg`, `overlap_fraction_density.svg`, `sequence_identity_density.svg`

## Analysis
### Utility scripts
`stola-compare-cross-orthology` (installed with the package): collects orthology info across multiple runs of the same reference against different queries. Organized into run modes (see `--help` on each subcommand):
- `summary`: build the combined per-species orthology-category TSV
- `sequence`: extract FASTA sequence of shared one2one orthologs
- `cluster`: compute the D statistic of phylogenetic signal (Fritz & Purvis 2010) for CNEE loss on a tree

#### Input
The `run_dirs.txt` file should be a list of full paths to output directories from multiple runs, one directory per line. 

You may optionally include a third column (tab separated) with a species prefix to use for that sample; otherwise, the script will attempt to use the query prefix present in the output files from that run. Example:

```
out/ref_vs_queryA  path/to/queryA.fasta|2bit  sppA
out/ref_vs_queryB  path/to/queryB.fasta|2bit
out/ref_vs_queryC  path/to/queryC.fasta|2bit sppC
```

#### summary
Useage:
```
stola-compare-cross-orthology summary \\
    --query-dirs run_dirs.txt \\
    --ref-bed reference_cnees.bed \\
    --ref-prefix ref \\
    --out-dir .
```

`summary` mode also writes two bar-chart SVGs alongside the combined TSV: how many reference CNEEs were lost in exactly 0, 1, 2, ... of the samples (`<ref-prefix>.loss_frequency_distribution.svg`), and the same for multi-projection (`<ref-prefix>.multiproj_frequency_distribution.svg`).

#### sequence
Useage:
```
stola-compare-cross-orthology sequence \\
    --query-dirs run_dirs.txt \\
    --orthology-table ref_combined_orthology.tsv \\
    --ref-bed reference_cnees.bed \\
    --ref-prefix ref \\
    --ref-genome path/to/ref.fasta|2bit \\
    --min-length 40 \\
    --max-seqs 300000 \\
    --max-lost-frac 0.3 \\
    --max-multi-frac 0.1 \\
    --out-dir .
```

`sequence` mode always outputs one2one orthologs, and reads its per-sample orthology categories from the TSV that `summary` mode already built for the *same* `run_dirs.txt` (pass it via `--orthology-table`) rather than recomputing them. Include the argument `--ref-genome` with the path to the reference 2bit/FASTA genome, and in the `run_dirs.txt` file add a second column (tab separated) for each sample with a path to the query 2bit/FASTA genome. This will generate one FASTA file for each ortholog set, which are output to subdirectories contained in `one2one_ortholog_fasta/`. A CNEE qualifies for extraction if it's lost in at most `--max-lost-frac` of samples AND multi-projection in at most `--max-multi-frac` of samples (a sample where it isn't one2one gets an all-N placeholder sequence in its output file rather than being excluded outright).
By default the script outputs only CNEEs >40bp (`--min-length`) and a maximum of 300000 individual CNEE FASTA files (`--max-seqs`). Set either of these two parameters to `0` to disable filtering and/or output all shared orthologs. 

#### cluster
Useage:
```
stola-compare-cross-orthology cluster \\
    --orthology-table ref_combined_orthology.tsv \\
    --tree species.nwk \\
    --ref-prefix dummyRef \\
    --min-lost 3 \\
    --n-perms 10 \\
    --out-dir .
```

`cluster` mode computes, per ref CNEE lost in at least `--min-lost` species (default 3; always skipped if lost in 0 or all species, since D is undefined for an invariant trait), the Fritz & Purvis D statistic of phylogenetic signal for its present/lost pattern on a given tree.

The tree's tip labels must exactly match the species labels used as columns in `--orthology-table`. Writes `<ref-prefix>_phylo_d.tsv` (`ref_cnee_id`, `n_lost`, `D`) plus a density plot of D values at `figures/<ref-prefix>.phylo_d_density.png`. D≈0 indicates a clumped/Brownian-motion-like loss pattern (phylogenetic signal); D≈1 indicates no phylogenetic signal (random with respect to the tree). `--n-perms` controls how many permutations/simulations build each of the two null distributions (random, Brownian motion) — these are computed once per distinct count of lost species (not per CNEE), since both null distributions depend only on that count, not on which particular species are lost.
