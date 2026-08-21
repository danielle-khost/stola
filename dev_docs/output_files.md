# Pipeline output

Everything below is written under a single `-o/--outdir` (default: current
directory) -- there's no way to redirect an individual file to a different
location. Intermediates (chain index, candidate TSV, FASTA byte-offset
indexes) live under `<outdir>/work/`.

## Main output
- `overlap_stats.tsv` -- overlap stats (raw)
- `overlap_stats.summary.tsv` -- overlap stats (summary)
- `query_cnees.bed` -- BED file of query projections (final)
- FASTA of query projections -- not yet implemented
- `query_ortho_status.tsv` -- TSV query to ref key (all projections)
  - query ID, orthologous ref ID, query orthology category
- `ref_ortho_status.tsv` -- TSV ref CNEEs
  - ref ID, ref orthology category, number projections (final), number pruned projections, pre-pruning ref orthology category
- `pruned_query_projections.tsv` -- per-projection log of synteny-verification prunes

## "Development output"
Output from development/exploratory output, keep all this for now:
`overlap_stats.strict.tsv`, `overlap_stats.strict.summary.tsv`,
`overlap_stats.synteny.tsv`, `overlap_stats.synteny_neighbors.tsv`. 