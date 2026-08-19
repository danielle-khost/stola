# Objective
After building graph using all projections and clustering to query loci, verify projections in certain categories using synteny criteria. Projections that meet verification criteria are kept in the graph, failed projections are pruned (and query clusters re-defined if necessary). Then call orthology on the pruned ref -> query graph.

Defining "solid" with same logic as rest of script, i.e. `sequence ID > 0.7, gap fraction < 0.1`

## Categories to verify
- Single ref -> single query where projection isn't solid (i.e. putative "one2one(weak)")
- Multi-projection cases (i.e. mutative "one2many"/"many2one"/"complex")
- Lost

"Synteny matches" means that the candidate/nearest neighbor query CNEE is orthologous to the reference CNEE  


## Verification criteria v2
For reference CNEEs that are weak one2ones or multi-projection, look at each projection that CNEE has and give it two scores (stored as edge attributes):

### Alignment score
```
align_score = sequence identity / max(sequence identity other connected projections)
```

- For one2one(weak), align_score = 1.0 (as only projection)
- For lost, align score = 0.0 (no aligned bases)
- For multi-projection cases, calculate as above
  - Use projections from all connected references 
    - e.g. for a many2one, look at projections from other refs that touch same query locus
    - e.g. for complex, do same as many2one and also projections from other refs that touch other queries

### Synteny score
```
Fraction of synteny matching nearest 5 neighboring CNEEs on each flank (trusted neighbors only)
```
"Trusted" neighbors = refs that are NOT themselves a verification target, i.e. already a solid one2one -- this avoids counting a soon-to-be-pruned neighbor's projection against (or for) the candidate. In practice this is the only pool of refs that can never be pruned, since a solid lone one2one edge is never re-examined.

The mean-flank-sequence-identity metric explored alongside this one turned out to be unhelpful in practice and is not used.

### For "lost" CNEEs
Treated differently from other categories, and **not currently verified** -- the criteria below only apply to one2one(weak)/one2many/many2one/complex/many2many. Lost CNEEs need different logic entirely (from lost CNEE in ref, number of bases out in ref to first solid 1:1, vs orthologous locus in query calc the number of bases out to first solid 1:1) and are left alone for now.

## Projection pruning (implemented)
For every projection in a category needing verification (i.e. NOT lost), retain it only if BOTH:
- `align_score > 0.5`
- `neighbor_fraction > 0.5`, where neighbor_fraction is the synteny score above computed over the nearest 5 trusted neighbors on EACH flank (10 checked total). A projection with no trusted neighbor on either flank fails this criterion outright -- there's no synteny evidence either way, so it can't be confirmed.

Culled projections are logged to a `culled_projections.tsv` (projection name, reference ortholog, initial pre-pruning orthology category, query coordinates, sequence identity, gap fraction) for tracking. A ref left with no surviving projection at all ends up `uncertain_loss` in the final orthology output, same as any other unverifiable ref.

Once verification has run, the one2one(weak) label is dropped from final output: any surviving lone ref/query edge -- whether it was always solid or only survived by passing verification -- is reported as plain `one2one`. The weak/solid distinction only still appears in the *initial* category recorded per culled projection.

Implemented in `src/toga4cnee/orthology.py` (`_verify_and_prune_projections`, `_score_projection`); thresholds and neighbor count are exposed as `--min-align-score`/`--min-neighbor-fraction`/`--verification-num-neighbors`.

