# Utility: compare_cross_orthology (`toga4cnee-compare-cross-orthology`, src/toga4cnee/cross_orthology.py)
Goal: take output across multiple runs of pipeline using same reference vs different queries and for each reference CNEE collect orthology status in each run + summary stats, optionally extract FASTA sequence of shared one2one orthologs.

MODULARIZE INTO MODES:
    - `summary`: generate orthology table and summary plot + stats
    - `sequence`: extract FASTA sequence of shared orthologs
    - `cluster`: use tree to calculate phylogenetic signal (D statistic) of CNEE loss

## Input
### `summary` mode
Sample sheet of output directories from multiple runs (one per line)
    - column 1: path to outdir
    - column 2 (optional): path to query FASTA/2bit genome
    - column 3 (optional): query prefix

#### Species prefix labeling
For each run, in order of preference:
1. The manual label from the list file's third column, if given.
2. Otherwise, peek at query_ortho_status.tsv and use the query_id_prefix
   baked into its query_cnee_id values -- unless that's "q" (the pipeline's
   generic default, uninformative across runs) or there's no usable
   (non-"NA") query id at all.
3. Otherwise, the output directory's own folder name.

### `sequence` mode
Sample sheet (as above)
Orthology table (output from `summary` mode)
Reference BED file
Reference genome (FASTA/2bit)
Reference prefix -- also used as the reference sequence's FASTA header

### `cluster` mode
Orthology table (output from `summary` mode)
Tree file -- phylogenetic tree of species in Newick format; tip labels must exactly match
the species labels used as columns in the orthology table
Reference prefix -- used for the output filenames


## Options
### All modes
`--out-dir [./]`: path to desired output directory

### `sequence` mode
`--max-lost-frac [0.3]`: maximum fraction of species CNEE can be lost in
`--max-multi-frac [0.1]`: maximum fraction of species CNEE can be multi-projection in
If EITHER of the flags violated (i.e. lost fraction > max-lost-frac OR multi-proj fraction > max-multi-fraction), omit that CNEE
`--min-length [40]`: minimum CNEE length to output as FASTA (checks length reference sequence)
`--max-seqs [300000]`: maximum number of CNEEs to extract to FASTA files

### `cluster` mode
`--min-lost [3]`: only compute D for a ref CNEE lost in at least this many species (always
skipped if lost in 0 or all species -- D is undefined for an invariant trait)
`--n-perms [10]`: number of permutations/simulations per distinct lost-species-count group,
for each of the two null distributions (random, Brownian motion)
`--seed [0]`: random seed, for reproducible permutations/simulations

## Output
### `summary` mode
TSV of orthology relations, where each row is ref CNEE and each column is a query species. E.g.:

```
ref_ID  sppA    sppB    sppC    ...
ce1 one2one one2one one2one ...
ce2 one2one uncertain_loss  one2many    ...
```

#### Summary figures
Density plot: counts of CNEE loss frequency across samples (i.e. X cnees lost in 1 spp, Y lost in 2 spp, etc.)
Density plot: counts of CNEE multi-projection frequency across samples

### `sequence` mode
Using the orthology table that was created from `summary` mode, for each ref CNEE that is one2one in sufficient samples (determined by `--max-lost-frac` and `--max-multi-frac`), extract the FASTA sequence for each ortholog from every sample (using the query BED file for coordinates) and combine into file `[ref CNEE ID]_orthologs.fasta` 

Organization: as this will create a lot of FASTA files, output them to a subdir `one2one_ortholog_fasta`; create sub-subdirs in that subdir named `batch1, batch2, etc` that holds at max 10000 FASTA files (i.e. if current `batch` dir is full, create the next one and output there) 

Efficiency: reuses each run's own pipeline-built FASTA index
(`<outdir>/work/query_fasta.idx`) instead of re-indexing, so column 2's
genome path must be the same FASTA/2bit file that run was originally given.
Because this needs query_ortho_status.tsv's actual data (not just its
prefix) to locate each one2one query ortholog, a manual label (column 3) no
longer lets a run skip having that file when FASTA output is on -- it's
still required, and a run missing it is skipped with a warning (same as a
missing ref_ortho_status.tsv).

## Cluster module (implemented)
Goal: calculate the Fritz & Purvis (2010) D statistic of phylogenetic signal for a binary
trait (CNEE present / CNEE lost) on a given tree, per ref CNEE.

Uses `toytree` (tree traversal, `tree.pcm.simulate_continuous_bm` for the Brownian null)
and `pandas`/`numpy`. Implementation: `run_phylo_d_analysis` / `_compute_d` /
`_simulate_random_d_values` / `_simulate_brownian_d_values` in `cross_orthology.py`.

```
#FXN: calculate d (sister clade distance)
##node values
traverse internal tree nodes (starting from lowest, i.e. postorder):
    get child nodes with trait values
    calculate node value = (unweighted) mean of values of child nodes
##sister clade distances
for each branch in tree:
    calculate d = abs(difference between connected node values)
sum d over every branch

import and index newick tree once; validate tree tip labels == orthology table's species columns

for each ref CNEE in the orthology table:
    get present (trait value = 1) or absent (value = 0) for each species
    if any species is "NA" (missing data): drop this CNEE
    n_lost = count of absent/lost species
    if n_lost == 0 or n_lost == n_species: skip (invariant trait, D undefined)
    if n_lost < --min-lost: skip

    #Calculate d(obs)
    label tree tip nodes with the observed trait values; calculate d(obs)

    #Calculate d(random) and d(brownian) -- ONCE PER DISTINCT n_lost VALUE, not per CNEE
    # (both null distributions depend only on the COUNT of lost species, not which species
    # are lost, so this collapses what looks like a per-CNEE-expensive step into one
    # precomputation per distinct n_lost group -- the key to tractability at CNEE-table scale)
    if n_lost not already cached:
        #d(random): no phylogenetic structure
        for --n-perms permutations:
            randomly shuffle which n_lost tips are "lost" (count fixed, identity random)
            calculate d; collect into d(random) distribution
        #d(brownian): trait evolves under Brownian motion
        for --n-perms simulations (toytree.pcm.simulate_continuous_bm, tips_only=True):
            threshold the simulated continuous tip values so that the n_lost LOWEST become
            "lost" (0) and the rest "present" (1) -- i.e. quantile-matched to the OBSERVED
            prevalence for this n_lost group, NOT a fixed q=0.5 median. (A fixed median would
            only be valid for traits lost in ~50% of species; real CNEE losses are usually a
            rare-event trait, lost in a small fraction of species, so the null must match that
            prevalence to be meaningful.)
            calculate d; collect into d(brownian) distribution
        cache mean(d(random)), mean(d(brownian)) for this n_lost value

    #CALCULATE D FOR THIS CNEE
    D = [d(obs) - mean(d(brownian))] / [mean(d(random)) - mean(d(brownian))]
```

### `cluster` mode output
TSV `[ref-prefix]_phylo_d.tsv`: columns `ref_cnee_id`, `n_lost`, `D`
Density plot (`figures/[ref-prefix].phylo_d_density.png`) of D values across all computed CNEEs,
with reference lines at D=0 (clumped/BM-like) and D=1 (random, no phylogenetic signal)

