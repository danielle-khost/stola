# Pipeline output

Everything below is written under a single `-o/--outdir` (default: current
directory) -- there's no way to redirect an individual file to a different
location. Intermediates (chain index, candidate TSV, FASTA/2bit byte-offset
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
- `<cnee_bed_stem>_named.bed` -- only written when the input reference CNEE
  BED had no name column; a BED4 copy naming each row with the cnee_uid the
  pipeline synthesized for it, so the same ids can be recognized across
  separate runs of the same reference BED (e.g. against different query
  species)

## "Development output"
Exploratory/diagnostic output, not considered main pipeline output -- lives
under `<outdir>/work/` alongside the other intermediates rather than in
`--outdir` itself: `work/overlap_stats.synteny.tsv`,
`work/overlap_stats.synteny_neighbors.tsv`.

"Strict" (solid-projections-only) output was removed -- it had fallen out of
use and wasn't adding signal beyond the all-projections output.