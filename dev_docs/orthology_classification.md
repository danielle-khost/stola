# GOAL
For each candidate projection to query, give it a query ID and classify it into different orthology relation categories between each reference and query CNEEs

# Filtering criteria
Overlap fraction, sequence similarity, number of projections for the CNEE
- provides an "identity threshold" 

# Orthology categories:
## Ref CNEE categories
- one2one: CNEE single copy in both query and ref; passes filtering criteria
- one2many: single copy in ref, multicopy in query; query projections pass filtering criteria
- many2one: multiple ref CNEEs align to same query locus; ref alignments pass filtering criteria
- complex: complex cases where multiple ref CNEEs overlap multiple the same query loci
- unaligned: ref CNEE has NO chain overlaps
- probable loss: ref CNEE overlaps a chain but only in chain gap AND chain gap not a contig boundary
- assembly gap loss: ref CNEE overlaps a chain but only in chain gap and gap is a contig boundary
- ALT: missing/loss: if cannot infer contig boundaries, collapse "probable loss" and "assembly gap loss" into this category
- uncertain loss: ref CNEE has only pruned projections

## Query CNEE categories
- one2one: CNEE single copy in both query and ref; passes filtering criteria
- one2many: single copy in ref, multicopy in query; query projections pass filtering criteria
- many2one: multiple ref CNEEs align to same query locus; ref alignments pass filtering criteria
- complex: complex cases where multiple ref CNEEs overlap multiple the same query loci
- pruned: query projection fails verification filter

## Loss sub-categorization (implemented)
Whether a chain gap is a contig boundary is inferred from the query FASTA
alone (no AGP): scan for N-runs, and see whether one (or a small handful of)
length(s) clearly dominate the runs of at least `--min-gap-length` bases --
see `assembly_gaps.py`. Three outcomes:
- **high confidence** (one dominant length): only that length counts as a
  boundary marker.
- **medium confidence** (a small set of 2-3 lengths together dominate):
  any of those lengths counts as a boundary marker -- some assemblies
  legitimately use more than one fixed-length gap type (e.g. contig-gap vs
  scaffold-gap).
- **none** (too spread out, or too few qualifying N-runs, or no query FASTA
  supplied at all): the probable_loss/assembly_gap_loss distinction can't be
  made confidently, so both collapse into `missing/loss` instead of
  guessing. A given ref's own chain gap can also individually be ambiguous
  (no aligned block on one side within its own chain) and falls back to
  `missing/loss` even when the genome's overall convention is confident.

This sub-categorization (plus per-ref pruning counts and pre-pruning
category) is independent of synteny verification and always runs when a
query FASTA is supplied, written to the ref-CNEE TSV described in
output_files.md. It does not change out_orthology's own "lost"/"unaligned"
labels.

# Implementation
- Create network graph where:
  - each node = query or ref CNEE 
- Draw an edge between ref node to query node (i.e. each alignment)
  - Store edge attributes:
    - sequence_identity
    - fraction_gap
  - A "solid" edge passes the current filtering thresholds
- for each query node:
  - get connected ref nodes
  - use info re: connected nodes to assign orthology category and ID for that query node (i.e CNEE):
    - if query node has at least 1 solid edge:
      - if connected to single ref node that has no other connections, then: "one2one"
        - set query ID as "q_[ref CNEE name]"
      - if connected single ref node and:
        - ref node connects multiple query nodes AND those query nodes only connect this reference
        - then: "one2many"
          - set query ID as "q_[ref DNEE name]_[n]"
      - if connects to multiple ref nodes and those ref nodes only connect this query node, then: "many2one"
      - for now, if none of these criteria met, then: "complex"
    - else if query node only has non-solid edge(s), then: "uncertain_loss"
- for each reference node with 0 chain-overlap candidates:
  - set query CNEE ID to "NA"
  - set orthology category as: "missing"
- for each reference node with chain-overlap candidates but no aligned query bases:
  - set query CNEE ID to "NA"
  - set orthology category as: "lost"

# Output
Tab separated list:
- query CNEE w/ ID
- partner ref CNEE(s) as comma separated list
- orthology category

# Orthology flowchart

"Solid" match for ref vs query projection = sequence identity > 0.7 AND gap fraction < 0.1.
If genome FASTAs are not provided, sequence identity is unavailable and the
solid-match decision uses gap fraction only.

```mermaid
flowchart TD

A["For each ref CNEE"] --> B{"Does it overlap a chain?"}
B -- "No" --> C["Missing"]
B -- "Yes" --> D{"Does it have aligned query bases?"}
D -- "No" --> E["Lost"]
D -- "Yes" --> F{"Does it have >=1 solid query projection?"}
F -- "No" --> O["Uncertain Loss"]
F -- "Yes" --> G{"Exactly 1 solid query?"}
G -- "Yes" --> H{"Does the 1 query match other refs?"}
H -- "No" --> I["one2one"]
H -- "Yes" --> J{"Do those other refs match other queries?"}
J -- "No" --> K["many2one"]
J -- "Yes" --> L["complex"]
G -- "No" --> M{"Do all multiple query matches match only this ref?"}
M -- "Yes" --> N["one2many"]
M -- "No" --> L
```

# Orthology flowchart v2
```mermaid
flowchart TD

A["For each ref CNEE"] --> B{"Does it overlap a chain?"}
B -- "No" --> C["Missing"]
B -- "Yes" --> D{"Does it have aligned query bases?"}
D -- "No" --> E["Lost"]
D -- "Yes" --> F{"Does it have exactly 1 query projection?"}
F -- "Yes" --> G{"Does the 1 query match other refs?"}
G -- "No" --> H{"Is it a solid projection?"}
H -- "Yes" --> H1["one2one"]
H -- "No" --> H2["one2one(weak)"]
G -- "Yes" --> I{"Do those other refs match other queries?"}
I -- "No" --> I1["many2one"]
I -- "Yes" --> I2["complex"]
F -- "No" --> J["Do all the multiple queries match only this ref?"]
J -- "No" --> I2
J -- "Yes" --> K["one2many"]

```
