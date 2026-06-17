# GOALS
- Project annotation of conserved nonexonic elements (CNEEs) from reference species with annotation onto query species without annotation using genome alignment
- Inspired by [TOGA2](https://github.com/hillerlab/TOGA2) which uses genome alignment to transfer gene annotations in a similar fashion
  - Use this as the "overall framework" to base this project on
- Want to output query CNEE annotation as well as categorize orthology classes of CNEEs in query relative to reference, similar to TOGA
  - I.e. one2one, one2many, many2one, missing, etc.
- Also want to classify mutations in query CNEEs vs reference

# INPUTS
- Reference and query genomes (FASTA or 2bit format)
- Whole genome alignment of query and reference in chain file format
- Annotation of CNEEs in reference species in BED format

# WORKFLOW LOGIC
- (Based off of TOGA2's workflow)
- *Build index* of chain file for fast lookup
  - TOGA uses BST index
- Get *candidate chains* x CNEE pairs:
  - Intersect CNEE annotation with chain alignment to identify candidate chains, i.e. chains that overlap CNEE annotation
- For each candidate chain pairs, calculate *overlap stats*:
  - position and size of CNEE overlap
  - fraction of CNEE covered by aligned chain
  - fraction of CNEE in gaps in chain
  - flanking region of CNEE coverage
  - chain score
  - Number of other chains overlapping same CNEE (maybe?)

OLD:
~ - *Classify* each candidate chain pairs based on overlap stats
  - This is what TOGA does, but not sure if relevant for CNEEs or not
    - Will need to consider what makes sense here!
  - ML classifier separates "orthologous" from "non-orthologous", then separates non-ortho into "paralogs", "spanning" or "pseudogene" (not relevant to CNEEs)
- For classified candidate chain pairs, perform *local realignment* within chain
  - TOGA aligns exons, intron, flanking regions, which are not relevant to CNEEs
  - For us: local realign CNEES and regions flanking CNEEs to get synteny info
  - Characterize mutation, etc.
- Use realignment to *create CNEE projections* to query species ~
  
- Realignment step not necessary for CNEEs: just use the initial chain alignment info
  - Instead, will need to have some kind of filtering/decision making step
    - Do not implement this yet, for now just report all alignments
- ~Construct graph between reference and query CNEEs, use connections to~ *Infer orthology* relations btwn reference and query CNEEs based on overlap stats with TBD filtering logic
  - one2one, one2many (duplicated in query), many2one (collapsed in query), missing, etc.
  - Also infor mutations (substitutions, indels, etc) in query CNEEs

# OUTPUT
- BED file of CNEEs in query
- TSV file matching query CNEEs to reference CNEEs, with ortho relationships
- TSV of mutations in query CNEEs