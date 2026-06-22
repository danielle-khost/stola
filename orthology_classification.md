# GOAL
For each candidate projection to query, classify it into different orthology relation categories between each reference and query CNEEs

# Filtering criteria
Overlap fraction, sequence similarity, number of projections for the CNEE
- provides an "identity threshold" 

# Possible categories:
- one2one: CNEE single copy in both query and ref; passes filtering criteria
- one2many: single copy in ref, multicopy in query; query projections pass filtering criteria
- many2one: multiple ref CNEEs align to same query locus; ref alignments pass filtering criteria
- missing: ref CNEE has NO chain overlaps
- "lost/diverged": ref CNEE has ONLY low quality alignment(s) to query; BELOW filtering threshold
- many2many/"complex": complex cases where multiple ref CNEEs overlap multiple the same query loci
  - i.e. a "combo" of many2one and one2many

# Implementation
- Create network graph where:
  - each node = query or ref CNEE 
- Draw an edge between ref node to query node (i.e. each alignment)
  - Weight the edge based on whether pass filtering criteria
    - passes threshold, weight = 1
    - below threshold, weight = 0.5 (placeholder)
- for each query node:
  - get connected ref nodes
    - check if connected ref nodes have edges with other query nodes
  - use info re: connected nodes to assign orthology category for that query node (i.e CNEE):
    - if query node has connections (edges) with weight = 1 (i.e. pass the filtering threshold):
      - if connected to single ref node that has no other connections, then: one2one
      - if connected single ref node and:
        - ref node connects multiple query nodes AND those query nodes only connect this reference
        - then: one2many
      - if connects to multiple ref nodes and those ref nodes only connect this query node, then: many2one