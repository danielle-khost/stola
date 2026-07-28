# Objective
After building graph using all projections and clustering to query loci, verify projections in certain categories using synteny criteria. Projections that meet verification criteria are kept in the graph, failed projections are pruned (and query clusters re-defined if necessary). Then call orthology on the pruned ref -> query graph.

Defining "solid" with same logic as rest of script, i.e. `sequence ID > 0.7, gap fraction < 0.1`

## Categories to verify
- Single ref -> single query where projection isn't solid (i.e. putative "one2one(weak)")
- Multi-projection cases (i.e. mutative "one2many"/"many2one"/"complex")

## Verification criteria
For a projection to be retained, it must meet at least 2 of the 4 following criteria:
- solid projection
- synteny matches on left flank
- synteny matches on right flank
- alignment to nearest flanking neighbors is solid

"Synteny matches" means that the candidate/nearest neighbor query CNEE is orthologous to the reference CNEE  