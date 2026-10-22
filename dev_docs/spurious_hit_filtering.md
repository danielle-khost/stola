# Goal
For all multi-projection CNEEs (one2many, many2one, complex), identify which reference projections are spurious (likely repeat-driven) hits, incorporating comparisons to "competing chains" overlapping query loci.

Example: refA is a potential one2many with chains to query_locus1 and query_locus2. At each query locus check if other chains overlap and how their alignment metrics compare to candidate projections.

### Filtering algo

Modified relative scoring algo:
```
for each query locus window:
    get chains with aligned blocks overlapping CNEE
    soft depth filter: if aligned block depth over window > N (20):
         sort chains by chain score
         retain only the top N best-scoring chains
    if query locus only covered by single chain:
        retain if 50% flanking sequence is aligned
    slice chains to get CNEE window and 10kb query flanks
    get "local score" from slices: seq identity over CNEE x aligned base fraction
    calc relative score for each slice: score / (score max scoring slice)
    for each slice:
        if relative score > 0.5, retain projection
```