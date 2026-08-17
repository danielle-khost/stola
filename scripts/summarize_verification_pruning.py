#!/usr/bin/env python
"""Summarize how synteny verification/pruning re-categorized each CNEE.

For every reference CNEE whose *initial* (pre-verification, pre-pruning)
connected-component shape would have made it one2one(weak)/one2many/
many2one/complex -- i.e. the categories classify_orthology's verification
step actually checks -- looks up what category it ended up as in the
*final*, already-computed overlap_stats.orthology.tsv, and cross-tabulates
initial category x final category as a table.

This does NOT re-run verification; it just re-derives the "would-be" initial
category the same way orthology.py's classify_orthology does internally
(from the raw, unpruned graph), and compares it against whatever orthology.tsv
you already have on disk. Categories are matched to whichever
--min-sequence-identity/--max-fraction-gap were used to generate that file --
using different thresholds here will produce a nonsensical comparison.

Example:
    PYTHONPATH=src python scripts/summarize_verification_pruning.py \\
        --cnee-bed tests/gGal7_cnee_final_data_filtered.fixed.bed \\
        --stats tests/test_runs/colLivi/overlap_stats.tsv \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from toga4cnee.orthology import (
    _build_orthology_graph,
    _cluster_query_projections,
    _edge_is_solid,
    _node_ref_id,
    _read_projection_evidence,
)

INITIAL_CATEGORIES = ["one2one(weak)", "one2many", "many2one", "complex/many2many"]
FINAL_CATEGORIES = [
    "one2one",
    "one2one(weak)",
    "one2many",
    "many2one",
    "complex/many2many",
    "uncertain_loss",
    "lost",
    "missing",
]


def _initial_category_by_ref(
    nx, graph, min_sequence_identity: float, max_fraction_gap: float
) -> dict[str, str]:
    # Mirrors classify_orthology's internal "does this component need
    # verification" shape check, but keeps every ref's would-be category
    # instead of collapsing straight to a keep/prune decision.
    category_by_ref: dict[str, str] = {}
    for component in nx.connected_components(graph):
        ref_nodes = [node for node in component if str(node).startswith("ref:")]
        query_nodes = [node for node in component if str(node).startswith("query:")]
        if not ref_nodes or not query_nodes:
            continue
        if len(ref_nodes) == 1 and len(query_nodes) == 1:
            edge_data = graph.edges[ref_nodes[0], query_nodes[0]]
            if not _edge_is_solid(edge_data, min_sequence_identity, max_fraction_gap):
                category_by_ref[_node_ref_id(ref_nodes[0])] = "one2one(weak)"
            # A solid lone edge is plain one2one -- not a category to verify.
        elif len(ref_nodes) == 1:
            category_by_ref[_node_ref_id(ref_nodes[0])] = "one2many"
        elif len(query_nodes) == 1:
            for node in ref_nodes:
                category_by_ref[_node_ref_id(node)] = "many2one"
        else:
            for node in ref_nodes:
                category_by_ref[_node_ref_id(node)] = "complex/many2many"
    return category_by_ref


def _final_category_by_ref(orthology_path: str | Path) -> dict[str, str]:
    category_by_ref: dict[str, str] = {}
    with Path(orthology_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_cnee_id", "partner_ref_cnees", "orthology_category"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Orthology file is missing columns: {', '.join(missing)}")
        for row in reader:
            category = row["orthology_category"]
            for ref_id in row["partner_ref_cnees"].split(","):
                if ref_id:
                    category_by_ref[ref_id] = category
    return category_by_ref


def _print_table(counts: dict[str, dict[str, int]], out_path: Path | None) -> None:
    present_final_categories = [
        category
        for category in FINAL_CATEGORIES
        if any(counts[initial].get(category, 0) for initial in INITIAL_CATEGORIES)
    ]
    header = ["initial_category", "initial_count"] + present_final_categories
    rows = []
    for initial in INITIAL_CATEGORIES:
        row_counts = counts[initial]
        initial_count = sum(row_counts.values())
        if initial_count == 0:
            continue
        row = [initial, str(initial_count)] + [
            str(row_counts.get(category, 0)) for category in present_final_categories
        ]
        rows.append(row)

    widths = [max(len(header[i]), *(len(row[i]) for row in rows)) for i in range(len(header))]
    def _fmt(cells: list[str]) -> str:
        return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths))

    print(_fmt(header))
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print(_fmt(row))

    if out_path is not None:
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        print(f"\nWrote {out_path}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnee-bed", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path, help="overlap_stats.tsv (all-mode, raw)")
    parser.add_argument(
        "--orthology", required=True, type=Path, help="Already-computed overlap_stats.orthology.tsv"
    )
    parser.add_argument("--min-sequence-identity", type=float, default=0.7)
    parser.add_argument("--max-fraction-gap", type=float, default=0.1)
    parser.add_argument("--out", type=Path, default=None, help="Optional TSV output path")
    args = parser.parse_args(argv)

    try:
        import networkx as nx
    except ImportError as exc:
        raise SystemExit("networkx is required (activate the biopython conda environment).") from exc

    from toga4cnee.bed import parse_bed

    ref_ids = [interval.cnee_uid for interval in parse_bed(args.cnee_bed)]

    print("Reading overlap stats and building the initial graph...", file=sys.stderr)
    evidence = _read_projection_evidence(args.stats)
    initial_clusters = _cluster_query_projections(evidence.projections)
    initial_graph = _build_orthology_graph(nx, ref_ids, initial_clusters)
    initial_category_by_ref = _initial_category_by_ref(
        nx, initial_graph, args.min_sequence_identity, args.max_fraction_gap
    )

    print("Reading final orthology categories...", file=sys.stderr)
    final_category_by_ref = _final_category_by_ref(args.orthology)

    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    unmatched = 0
    for ref_id, initial_category in initial_category_by_ref.items():
        final_category = final_category_by_ref.get(ref_id)
        if final_category is None:
            unmatched += 1
            continue
        counts[initial_category][final_category] += 1

    if unmatched:
        print(
            f"warning: {unmatched} ref CNEE(s) with an initial category had no entry in "
            "--orthology (bed/orthology mismatch?)",
            file=sys.stderr,
        )

    _print_table(counts, args.out)


if __name__ == "__main__":
    main()
