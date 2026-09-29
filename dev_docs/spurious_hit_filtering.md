# Goal
For all multi-projection CNEEs (one2many, many2one, complex), identify which reference projections are spurious (likely repeat-driven) hits, incorporating comparisons to "competing chains" overlapping query loci.

Example: refA is a potential one2many with chains to query_locus1 and query_locus2. At each query locus check if other chains overlap and how their alignment metrics compare to candidate projections.

# Exploration Workflow
1. Hard filter step: remove any "CNEE-only" projections, i.e. where fraction of chain length comprised of CNEE is > 0.1 (as likely represents repeat driven alignment and doens't have enough info to evaluate alignment)
2. For remaining multi-projections:
   1.  Calculate aligned block coverage over CNEE and flanking regions (10kb either side)
       1.  If flanks of chain are too short, treat "missing" bases as zero coverage (e.g. if chain has 5kb on each flank of 1.0 coverage, total covered fraction would be 0.5)
   2.  Calculate mean sequence identity over CNEE and flanking regions
3.  At each query locus that a multi-proj maps to, identify "competing chains", i.e. any other chains that overlap the same locus
    1.  Calculate the same stats
    2.  DEPTH FILTER: to avoid repeat-driven pileup, discard projections where locus depth >= 20
4.  Normalize the aligned block coverage and mean sequence identity stats for each chain (over CNEE and flanks) at each query locus relative to the highest score at that locus

# Output
- how many projections removed by the hard filter
- raw stats for each multi-proj
- stats normalized by competing chain stats
- distribution plots of normalized scores 

# Implementation Workflow 
OVERALL LOGIC: are there better alignments over a given candidate CNEE locus, as determined by coverage and alignment score
```
For each multi-proj CNEE
    For each query loci
        Get aligned block depth over candidate CNEE locus
            If depth over locus > filter threshold (default 20, tuneable), discard that projection (repeat-driven)
            If depth over locus = 1 (i.e. CNEE candidate chain is the ONLY chain with aligned blocks over candidate CNEE locus):
                If aligned block coverage fraction on flanks (10kb each side) > 0.5 (ANY coverage, does not have to be from this chain), accept projection
                Else low/no coverage on flanks, reject projection
            If 1 < aligned block depth over locus < depth filter threshold (i.e. there are other chains competing at this locus):
                Calculate alignment score over CNEE for candidate chain (`S(c)`)and competitors (`S(i...j)`)
                Get max score over CNEE (`S(max)`)
                Normalize scores by dividing each by `S(max)` to get score relative to best alignment (`R(c), R(i...j)`)
                    If `R(c)` = 1 and `R(i...j)` < 0.7: accept candidate
                    Else if `R(c)` < 0.7: reject candidate
                    Else if `R(c)` = 1 and `R(i...j)` >= 0.7 OR `R(c)` >= 0.7: test flanks
                        Calculate alignment score of flanks for candidate chain and competitors (missing seq counts as gaps)
                        If flank align score of candidate > competitors, accept projection
                        Else discard projection
```