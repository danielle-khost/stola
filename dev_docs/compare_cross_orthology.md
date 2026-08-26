# Utility script: compare_cross_orthology.py
Goal: take output across multiple runs of pipeline using same reference vs different queries and for each reference CNEE collect orthology status in each run.

## Input
List file of output directories from multiple runs (one per line)
    - column 1: path to outdir
    - column 2 (optional): path to query FASTA/2bit genome
    - column 3 (optional): query prefix
Reference BED file
Reference genome (FASTA/2bit) -- only needed for FASTA output (see below);
required unless --no-fasta is given
Reference prefix -- also used as the reference sequence's FASTA header

## Column label
For each run, in order of preference:
1. The manual label from the list file's third column, if given.
2. Otherwise, peek at query_ortho_status.tsv and use the query_id_prefix
   baked into its query_cnee_id values -- unless that's "q" (the pipeline's
   generic default, uninformative across runs) or there's no usable
   (non-"NA") query id at all.
3. Otherwise, the output directory's own folder name.

## Intended output
TSV of orthology relations, where each row is ref CNEE and each column is a query species. E.g.:

```
ref_ID  sppA    sppB    sppC    ...
ce1 one2one one2one one2one ...
ce2 one2one uncertain_loss  one2many    ...
```

### Output FASTA
Is true by default and requires 2nd column in input key with path to genome. Disable using `--no-fasta` flag.

Once collected list is created, for each ref CNEE that is one2one in **all** samples, extract the FASTA sequence for each ortholog from every sample (using the query BED file for coordinates) and combine into file `[ref CNEE ID]_orthologs.fasta` 

Organization: as this will create a lot of FASTA files, output them to a subdir `one2one_ortholog_fasta`; create sub-subdirs in that subdir named `batch1, batch2, etc` that holds at max 10000 FASTA files (i.e. if current `batch` dir is full, create the next one and output there) 

Efficiency: reuses each run's own pipeline-built FASTA index
(`<outdir>/work/query_fasta.idx`) instead of re-indexing, so column 2's
genome path must be the same FASTA/2bit file that run was originally given.
Because this needs query_ortho_status.tsv's actual data (not just its
prefix) to locate each one2one query ortholog, a manual label (column 3) no
longer lets a run skip having that file when FASTA output is on -- it's
still required, and a run missing it is skipped with a warning (same as a
missing ref_ortho_status.tsv).

## Pseudocode
```
for each output directory:
    peek at query_ortho_status.tsv to get species prefix (or generate ID if prefix is generic)
    for each ref CNEE in ref BED:
        get ortho status from ref_ortho_status.tsv in outdir; save it

output summary table [ref prefix]_combined_orthology.tsv

```

