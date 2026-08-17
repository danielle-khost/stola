#!/usr/bin/env python
"""Exploratory 2D histograms of proposed synteny-verification scores.

For each category that would need synteny verification (one2one(weak),
one2many, many2one, complex/many2many, lost), samples up to --sample-size
reference CNEEs, scores every one of their candidate projections, and plots
align_score vs. each of the two synteny-score components as a 2-panel PNG.

This is a standalone exploration tool -- it does NOT read an existing
orthology.tsv. Categories are re-derived directly from the raw overlap stats
(the same "initial, unpruned graph" shape orthology.py's classify_orthology
uses internally), since the point is to evaluate new verification criteria
*before* any pruning decision has been baked in.

Score definitions (see synteny_verification.md on the verification-dev branch):

  align_score = this projection's sequence identity / max sequence identity
                among every projection (including this one) in the same
                connected component (the whole component -- all refs and
                queries reachable from this one, not just the immediate
                edge). Always <= 1.0; the component's best projection lands
                at exactly 1.0. one2one(weak) with a single projection => 1.0.
                lost => 0.0 (no aligned bases).

  synteny_score (neighbors) = fraction of the nearest 5 *solid* reference
                neighbor CNEEs per side whose own query ortholog(s) fall
                anywhere among this candidate's nearest 5 query neighbors
                (looked up against a *solid-only* query-cluster index, to
                avoid nearby non-solid noise shadowing the true ortholog).

  synteny_score (flank align) = mean sequence identity of the reference
                flank out to the farthest solid neighbor found on each side
                (using this candidate's own chain's aligned blocks). For
                "lost" candidates (no aligned bases at all), this is instead
                the mean bp distance to the nearest aligned chain block on
                each side.

Requires matplotlib (`pip install matplotlib` or `pip install .[viz]`) and
networkx, in addition to the toga4cnee package itself.

Example:
    PYTHONPATH=src python scripts/plot_verification_scores.py \\
        --cnee-bed tests/gGal7_cnee_final_data_filtered.fixed.bed \\
        --stats tests/test_runs/colLivi/overlap_stats.tsv \\
        --chain-index tests/test_runs/colLivi/work/chain_index.json \\
        --ref-fasta tests/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna \\
        --query-fasta tests/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna \\
        --out-dir tests/test_runs/colLivi/verification_scores \\
        --sample-size 2000
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from toga4cnee.bed import parse_bed
from toga4cnee.chain import iter_chain_blocks
from toga4cnee.index import load_index
from toga4cnee.orthology import (
    Projection,
    _build_orthology_graph,
    _cluster_query_projections,
    _edge_is_solid,
    _node_ref_id,
    _query_neighbor_index_from_clusters,
    _read_projection_evidence,
    _ref_node,
)
from toga4cnee.sequence_compare import SequenceComparator, aligned_segments_for_interval
from toga4cnee.summary import is_solid_match
from toga4cnee.synteny import (
    BlockIndex,
    QueryNeighborIndex,
    RefNeighborIndex,
    _aligned_bases,
    _flanking_blocks,
    _query_flank_anchors,
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except ImportError as exc:  # pragma: no cover - environment guard, not test logic
    raise SystemExit(
        "matplotlib is required for this script (`pip install matplotlib` or "
        "`pip install .[viz]`)."
    ) from exc

CATEGORIES = ["one2one(weak)", "one2many", "many2one", "complex/many2many", "lost"]
CATEGORY_FILENAMES = {
    "one2one(weak)": "one2one_weak",
    "one2many": "one2many",
    "many2one": "many2one",
    "complex/many2many": "complex",
    "lost": "lost",
}

# Sequential single-hue ramp (blue, light->dark) per the dataviz skill's
# palette -- this is a magnitude (point-density) encoding, never a rainbow.
SEQUENTIAL_BLUE = LinearSegmentedColormap.from_list(
    "sequential_blue",
    ["#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#0d366b"],
)

NUM_NEIGHBORS = 5


@dataclass
class ScoredPoint:
    ref_id: str
    label: str
    align_score: float
    synteny_neighbors: float | None
    synteny_flank: float | None


# --------------------------------------------------------------------------
# Category derivation from the raw, unpruned graph (mirrors orthology.py's
# internal "does this component need verification" logic, but keeps every
# ref's category rather than only a keep/prune decision).
# --------------------------------------------------------------------------


def _category_by_ref_from_components(
    nx, graph, min_sequence_identity: float, max_fraction_gap: float
) -> dict[str, str]:
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


def _component_id_by_ref(nx, graph) -> dict[str, int]:
    component_id_by_ref: dict[str, int] = {}
    for component_id, component in enumerate(nx.connected_components(graph)):
        for node in component:
            if str(node).startswith("ref:"):
                component_id_by_ref[_node_ref_id(node)] = component_id
    return component_id_by_ref


# --------------------------------------------------------------------------
# Solid-only neighbor lookups (ref side skips non-solid CNEEs while walking;
# query side is restricted to a cluster index built from solid projections
# only -- both to avoid the "noisy initial cluster" false-negative bug).
# --------------------------------------------------------------------------


def _ref_has_solid_map(
    projections: list[Projection], min_sequence_identity: float, max_fraction_gap: float
) -> dict[str, bool]:
    has_solid: dict[str, bool] = defaultdict(bool)
    for projection in projections:
        if is_solid_match(
            projection.sequence_identity,
            projection.fraction_gap,
            min_sequence_identity,
            max_fraction_gap,
        ):
            has_solid[projection.ref_id] = True
    return has_solid


def _solid_ref_neighbor_chain(
    ref_index: RefNeighborIndex, ref_has_solid: dict[str, bool], ref_uid: str, side: str, n: int
) -> list:
    # Walk outward on the reference genome, skipping any CNEE that doesn't
    # itself have a solid projection, collecting up to n qualifying BedIntervals.
    chrom = ref_index.chrom_by_uid.get(ref_uid)
    if chrom is None:
        return []
    intervals = ref_index.intervals_by_chrom[chrom]
    position = ref_index.position_by_uid[ref_uid]
    found = []
    if side == "left":
        index = position - 1
        while index >= 0 and len(found) < n:
            candidate = intervals[index]
            if ref_has_solid.get(candidate.cnee_uid):
                found.append(candidate)
            index -= 1
    else:
        index = position + 1
        while index < len(intervals) and len(found) < n:
            candidate = intervals[index]
            if ref_has_solid.get(candidate.cnee_uid):
                found.append(candidate)
            index += 1
    return found


def _query_neighbor_chain(
    query_index: QueryNeighborIndex,
    record_by_node_id: dict[str, object],
    chrom: str,
    anchor: int,
    direction: str,
    n: int,
) -> list[str]:
    found = []
    current_anchor = anchor
    for _ in range(n):
        node_id, _distance = query_index.neighbor(chrom, current_anchor, direction)
        if node_id is None:
            break
        found.append(node_id)
        # Re-anchor on that cluster's own far edge to keep walking outward.
        match = record_by_node_id[node_id]
        current_anchor = match.start if direction == "before" else match.end
    return found


def _anchors_for_own_interval(query_chrom: str, query_start: int, query_end: int, query_strand: str):
    forward = query_strand == "+"
    if forward:
        return {"left": (query_chrom, query_start, "before"), "right": (query_chrom, query_end, "after")}
    return {"left": (query_chrom, query_end, "after"), "right": (query_chrom, query_start, "before")}


def _synteny_neighbors_score(
    ref_uid: str,
    anchors: dict[str, tuple[str, int, str] | None],
    ref_index: RefNeighborIndex,
    ref_has_solid: dict[str, bool],
    query_index: QueryNeighborIndex,
    record_by_node_id: dict[str, object],
    ref_to_solid_query_nodes: dict[str, set[str]],
) -> float | None:
    total_checked = 0
    total_matched = 0
    for side in ("left", "right"):
        anchor = anchors.get(side)
        ref_neighbors = _solid_ref_neighbor_chain(ref_index, ref_has_solid, ref_uid, side, NUM_NEIGHBORS)
        if not ref_neighbors:
            continue
        if anchor is None:
            total_checked += len(ref_neighbors)
            continue
        chrom, point, direction = anchor
        query_neighbor_nodes = set(
            _query_neighbor_chain(query_index, record_by_node_id, chrom, point, direction, NUM_NEIGHBORS)
        )
        for neighbor_interval in ref_neighbors:
            total_checked += 1
            neighbor_query_nodes = ref_to_solid_query_nodes.get(neighbor_interval.cnee_uid, set())
            if neighbor_query_nodes & query_neighbor_nodes:
                total_matched += 1
    if total_checked == 0:
        return None
    return total_matched / total_checked


def _window_identity(
    ref_chrom: str,
    window_start: int,
    window_end: int,
    record,
    blocks: BlockIndex,
    sequence_comparator: SequenceComparator,
) -> float | None:
    if window_end <= window_start:
        return None
    overlapping_blocks = list(blocks.overlapping(window_start, window_end))
    if not overlapping_blocks:
        return None
    metrics = sequence_comparator.compare(
        ref_chrom,
        window_start,
        window_end,
        record.query_chrom,
        record.query_strand,
        aligned_segments_for_interval(overlapping_blocks, window_start, window_end, record.query_strand),
    )
    compared = int(metrics["seq_compared_bases"])
    if compared == 0:
        return None
    return metrics["seq_matches"] / compared


def _synteny_flank_align_score(
    ref_chrom: str,
    ref_start: int,
    ref_end: int,
    left_neighbors: list,
    right_neighbors: list,
    record,
    blocks: BlockIndex,
    sequence_comparator: SequenceComparator,
) -> float | None:
    values = []
    if left_neighbors:
        window_start = left_neighbors[-1].start
        identity = _window_identity(ref_chrom, window_start, ref_start, record, blocks, sequence_comparator)
        if identity is not None:
            values.append(identity)
    if right_neighbors:
        window_end = right_neighbors[-1].end
        identity = _window_identity(ref_chrom, ref_end, window_end, record, blocks, sequence_comparator)
        if identity is not None:
            values.append(identity)
    if not values:
        return None
    return sum(values) / len(values)


def _synteny_flank_distance_score(
    ref_start: int, ref_end: int, blocks: BlockIndex
) -> float | None:
    preceding, following = _flanking_blocks(blocks, ref_start, ref_end)
    values = []
    if preceding is not None:
        values.append(max(0, ref_start - preceding.ref_end))
    if following is not None:
        values.append(max(0, following.ref_start - ref_end))
    if not values:
        return None
    return sum(values) / len(values)


# --------------------------------------------------------------------------
# Main.
# --------------------------------------------------------------------------


def _sample(items: list[str], size: int, rng: random.Random) -> list[str]:
    ordered = sorted(items)
    if len(ordered) <= size:
        return ordered
    return rng.sample(ordered, size)


def _plot_category(category: str, points: list[ScoredPoint], out_dir: Path, dpi: int) -> Path:
    align = [p.align_score for p in points]
    neighbors = [p.synteny_neighbors for p in points if p.synteny_neighbors is not None]
    align_for_neighbors = [p.align_score for p in points if p.synteny_neighbors is not None]
    flank = [p.synteny_flank for p in points if p.synteny_flank is not None]
    align_for_flank = [p.align_score for p in points if p.synteny_flank is not None]

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(f"Verification scores: {category}  (n={len(points)})", fontsize=13, fontweight="bold")

    if align_for_neighbors:
        h = ax_left.hist2d(align_for_neighbors, neighbors, bins=25, cmap=SEQUENTIAL_BLUE, cmin=1)
        fig.colorbar(h[3], ax=ax_left, label="count")
    ax_left.set_xlabel("align_score")
    ax_left.set_ylabel("synteny_score (neighbors)")
    ax_left.set_title(f"n={len(neighbors)}", fontsize=9, color="#555555")
    for spine in ("top", "right"):
        ax_left.spines[spine].set_visible(False)

    flank_label = "synteny_score (flank align, mean bp to nearest aligned block)" if category == "lost" else "synteny_score (flank align, mean sequence identity)"
    if align_for_flank:
        if category == "lost":
            # bp-distance spans several orders of magnitude (a handful of
            # candidates sit in very sparse chains) -- a linear axis squashes
            # the informative bulk of the distribution into a sliver, so bin
            # and display this one on a log scale instead.
            import numpy as np

            floored = [max(1.0, v) for v in flank]
            y_bins = np.geomspace(min(floored), max(floored), 26)
            h = ax_right.hist2d(align_for_flank, floored, bins=[25, y_bins], cmap=SEQUENTIAL_BLUE, cmin=1)
            ax_right.set_yscale("log")
        else:
            h = ax_right.hist2d(align_for_flank, flank, bins=25, cmap=SEQUENTIAL_BLUE, cmin=1)
        fig.colorbar(h[3], ax=ax_right, label="count")
    ax_right.set_xlabel("align_score")
    ax_right.set_ylabel(flank_label)
    ax_right.set_title(f"n={len(flank)}", fontsize=9, color="#555555")
    for spine in ("top", "right"):
        ax_right.spines[spine].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_path = out_dir / f"{CATEGORY_FILENAMES[category]}.png"
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnee-bed", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path, help="overlap_stats.tsv (all-mode, raw)")
    parser.add_argument("--chain-index", required=True, type=Path)
    parser.add_argument("--ref-fasta", required=True, type=Path)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path, default=None, help="FASTA index cache dir [default: alongside --chain-index]")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-sequence-identity", type=float, default=0.7)
    parser.add_argument("--max-fraction-gap", type=float, default=0.1)
    parser.add_argument("--sample-size", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args(argv)

    try:
        import networkx as nx
    except ImportError as exc:
        raise SystemExit(
            "networkx is required (activate the biopython conda environment)."
        ) from exc

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    ref_intervals = parse_bed(args.cnee_bed)
    ref_by_uid = {interval.cnee_uid: interval for interval in ref_intervals}
    ref_ids = [interval.cnee_uid for interval in ref_intervals]
    ref_index = RefNeighborIndex.from_bed(args.cnee_bed)

    print("Reading overlap stats...", file=sys.stderr)
    evidence = _read_projection_evidence(args.stats)
    refs_with_any_projection = {p.ref_id for p in evidence.projections}
    lost_ref_ids = evidence.refs_with_candidates - refs_with_any_projection

    initial_clusters = _cluster_query_projections(evidence.projections)
    initial_graph = _build_orthology_graph(nx, ref_ids, initial_clusters)
    category_by_ref = _category_by_ref_from_components(
        nx, initial_graph, args.min_sequence_identity, args.max_fraction_gap
    )
    component_id_by_ref = _component_id_by_ref(nx, initial_graph)

    projections_by_component: dict[int, list[Projection]] = defaultdict(list)
    for projection in evidence.projections:
        component_id = component_id_by_ref.get(projection.ref_id)
        if component_id is not None:
            projections_by_component[component_id].append(projection)

    ref_has_solid = _ref_has_solid_map(evidence.projections, args.min_sequence_identity, args.max_fraction_gap)
    solid_projections = [
        p
        for p in evidence.projections
        if is_solid_match(p.sequence_identity, p.fraction_gap, args.min_sequence_identity, args.max_fraction_gap)
    ]
    solid_clusters = _cluster_query_projections(solid_projections)
    solid_query_neighbor_index = _query_neighbor_index_from_clusters(solid_clusters)
    solid_cluster_by_node_id = {cluster.node_id: cluster for cluster in solid_clusters}
    ref_to_solid_query_nodes: dict[str, set[str]] = defaultdict(set)
    for cluster in solid_clusters:
        for ref_id in cluster.refs:
            ref_to_solid_query_nodes[ref_id].add(cluster.node_id)

    refs_by_category: dict[str, list[str]] = defaultdict(list)
    for ref_id, category in category_by_ref.items():
        refs_by_category[category].append(ref_id)
    refs_by_category["lost"] = sorted(lost_ref_ids)

    sampled_by_category = {
        category: _sample(refs_by_category.get(category, []), args.sample_size, rng)
        for category in CATEGORIES
    }
    for category in CATEGORIES:
        print(
            f"{category}: {len(refs_by_category.get(category, []))} available, "
            f"sampled {len(sampled_by_category[category])}",
            file=sys.stderr,
        )

    work_dir = args.work_dir if args.work_dir is not None else args.chain_index.parent
    sequence_comparator = SequenceComparator(args.ref_fasta, args.query_fasta, work_dir)
    payload, records = load_index(args.chain_index)
    chain_records_by_id = {record.chain_id: record for record in records}
    chain_path = payload["indexed_chain_path"]
    block_cache: dict[str, BlockIndex] = {}

    points_by_category: dict[str, list[ScoredPoint]] = defaultdict(list)

    try:
        with Path(chain_path).open("rb") as chain_handle:

            def _blocks_for(chain_id: str) -> BlockIndex:
                if chain_id not in block_cache:
                    record = chain_records_by_id[chain_id]
                    chain_handle.seek(record.start_byte)
                    chain_text = chain_handle.read(record.byte_length).decode("utf-8")
                    block_cache[chain_id] = BlockIndex.from_blocks(iter_chain_blocks(chain_text))
                return block_cache[chain_id]

            # --- non-lost categories: score every projection of every sampled ref ---
            sampled_non_lost = {
                ref_id
                for category in ("one2one(weak)", "one2many", "many2one", "complex/many2many")
                for ref_id in sampled_by_category[category]
            }
            projections_by_ref: dict[str, list[Projection]] = defaultdict(list)
            for projection in evidence.projections:
                if projection.ref_id in sampled_non_lost:
                    projections_by_ref[projection.ref_id].append(projection)

            done = 0
            for ref_id in sorted(sampled_non_lost):
                category = category_by_ref[ref_id]
                ref_interval = ref_by_uid[ref_id]
                left_neighbors = _solid_ref_neighbor_chain(ref_index, ref_has_solid, ref_id, "left", NUM_NEIGHBORS)
                right_neighbors = _solid_ref_neighbor_chain(ref_index, ref_has_solid, ref_id, "right", NUM_NEIGHBORS)
                component_id = component_id_by_ref[ref_id]
                component_projections = projections_by_component[component_id]

                for projection in projections_by_ref[ref_id]:
                    # max includes this projection itself, so the best-scoring
                    # projection in the component always lands at exactly 1.0.
                    max_identity_in_component = max(
                        (p.sequence_identity or 0.0) for p in component_projections
                    )
                    if max_identity_in_component <= 0:
                        align_score = 1.0
                    else:
                        align_score = (projection.sequence_identity or 0.0) / max_identity_in_component

                    anchors = _anchors_for_own_interval(
                        projection.query_chrom, projection.query_start, projection.query_end, projection.query_strand
                    )
                    neighbors_score = _synteny_neighbors_score(
                        ref_id,
                        anchors,
                        ref_index,
                        ref_has_solid,
                        solid_query_neighbor_index,
                        solid_cluster_by_node_id,
                        ref_to_solid_query_nodes,
                    )
                    record = chain_records_by_id[projection.chain_id]
                    blocks = _blocks_for(projection.chain_id)
                    flank_score = _synteny_flank_align_score(
                        ref_interval.chrom,
                        projection.ref_start,
                        projection.ref_end,
                        left_neighbors,
                        right_neighbors,
                        record,
                        blocks,
                        sequence_comparator,
                    )
                    points_by_category[category].append(
                        ScoredPoint(
                            ref_id=ref_id,
                            label=projection.projection_name,
                            align_score=align_score,
                            synteny_neighbors=neighbors_score,
                            synteny_flank=flank_score,
                        )
                    )
                done += 1
                if done % 500 == 0:
                    print(f"  scored {done}/{len(sampled_non_lost)} non-lost refs...", file=sys.stderr)

            # --- lost: score every zero-aligned candidate row of every sampled ref ---
            lost_sample = set(sampled_by_category["lost"])
            if lost_sample:
                print("Scanning stats for lost candidates...", file=sys.stderr)
                with Path(args.stats).open("r", encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    done = 0
                    for row in reader:
                        if row["cnee_uid"] not in lost_sample or int(row["aligned_bases"]) > 0:
                            continue
                        ref_id = row["cnee_uid"]
                        ref_interval = ref_by_uid[ref_id]
                        ref_start, ref_end = int(row["ref_start"]), int(row["ref_end"])
                        left_neighbors = _solid_ref_neighbor_chain(ref_index, ref_has_solid, ref_id, "left", NUM_NEIGHBORS)
                        right_neighbors = _solid_ref_neighbor_chain(ref_index, ref_has_solid, ref_id, "right", NUM_NEIGHBORS)
                        record = chain_records_by_id[row["chain_id"]]
                        blocks = _blocks_for(row["chain_id"])
                        anchors_raw = _query_flank_anchors(row, record, blocks)
                        anchors = {
                            side: ((record.query_chrom, point, direction) if point is not None else None)
                            for side, (point, direction) in anchors_raw.items()
                        }
                        neighbors_score = _synteny_neighbors_score(
                            ref_id,
                            anchors,
                            ref_index,
                            ref_has_solid,
                            solid_query_neighbor_index,
                            solid_cluster_by_node_id,
                            ref_to_solid_query_nodes,
                        )
                        flank_score = _synteny_flank_distance_score(ref_start, ref_end, blocks)
                        points_by_category["lost"].append(
                            ScoredPoint(
                                ref_id=ref_id,
                                label=row.get("query_projection_name", ref_id),
                                align_score=0.0,
                                synteny_neighbors=neighbors_score,
                                synteny_flank=flank_score,
                            )
                        )
                        done += 1
                        if done % 500 == 0:
                            print(f"  scored {done} lost candidates...", file=sys.stderr)
    finally:
        sequence_comparator.close()

    for category in CATEGORIES:
        points = points_by_category.get(category, [])
        if not points:
            print(f"warning: no points to plot for {category}", file=sys.stderr)
            continue
        out_path = _plot_category(category, points, args.out_dir, args.dpi)
        print(f"Wrote {len(points)} scored projections -> {out_path}")


if __name__ == "__main__":
    main()
