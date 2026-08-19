# Goal
Using chain x CNEE overlap, get the sequence for CNEE from reference and query FASTA files, calculate similarity scores and number of mismatches. Use the BioPython package.

# Pseudocode
library BioPython
from Bio import Align
alignment = Align.parse(*.chain, "chain")
reference = ref.fasta
query = query.fasta

for each candidate chain x CNEE pair
    get CNEE sequence (reference and query)
    report mismatches, gaps etc. in each aligned pair using `count` method (or similar)
    calculate alignment score for aligned pair using `aligner.score` method (or similar)
