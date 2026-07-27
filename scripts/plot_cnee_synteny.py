#!/usr/bin/env python
"""Render "synteny cartoon" PNGs for spot-checking reference CNEE orthology calls.

For each requested reference CNEE, draws a Ref: track (the CNEE and its
immediate left/right ref-genome neighbors) and a Query: track (the same for
its query-side projection(s)), colored by orthology category and shaped by
projection strength (rectangle = strong/solid, ellipse = weak). One2many and
many2one CNEEs get one row per query/ref partner respectively; uncertain_loss
and lost CNEEs fan out over every candidate alignment (or, for lost, every
candidate chain's interpolated query anchor).

Requires matplotlib (`pip install matplotlib`, or `pip install .[viz]` from
the repo root) in addition to the toga4cnee package itself.

Example:
    PYTHONPATH=src python scripts/plot_cnee_synteny.py \\
        --cnee-bed tests/gGal7_cnee_final_data_filtered.fixed.bed \\
        --query-bed tests/test_runs/colLivi/query_cnees.bed \\
        --stats tests/test_runs/colLivi/overlap_stats.tsv \\
        --orthology tests/test_runs/colLivi/overlap_stats.orthology.tsv \\
        --chain-index tests/test_runs/colLivi/work/chain_index.json \\
        --out-dir tests/test_runs/colLivi/synteny_cartoons \\
        ce3 ce44 ce6
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from toga4cnee.bed import parse_bed
from toga4cnee.chain import iter_chain_blocks
from toga4cnee.index import load_index
from toga4cnee.summary import stats_row_is_solid
from toga4cnee.synteny import (
    BlockIndex,
    OrthologyIndex,
    QueryBedRecord,
    QueryNeighborIndex,
    RefNeighborIndex,
    _flanking_blocks,
    _query_flank_anchors,
    _read_query_bed,
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse, Rectangle
except ImportError as exc:  # pragma: no cover - environment guard, not test logic
    raise SystemExit(
        "matplotlib is required for this script (`pip install matplotlib` or "
        "`pip install .[viz]`)."
    ) from exc

ALIGNED_CATEGORIES = {"one2one", "one2many", "many2one", "complex/many2many"}

CATEGORY_COLORS = {
    "one2one": "#4CAF50",
    "one2many": "#9C27B0",
    "many2one": "#00BCD4",
    "complex/many2many": "#9E9E9E",
    "uncertain_loss": "#FFC107",
    "lost": "#F44336",
    "missing": "#F5F5F5",
}
CATEGORY_LABELS = {
    "one2one": "one2one",
    "one2many": "one2many",
    "many2one": "many2one",
    "complex/many2many": "complex",
    "uncertain_loss": "uncertain loss",
    "lost": "lost",
    "missing": "missing",
}
CATEGORY_ORDER = [
    "one2one",
    "one2many",
    "many2one",
    "complex/many2many",
    "uncertain_loss",
    "lost",
    "missing",
]


# --------------------------------------------------------------------------
# Data model for one drawable box.
# --------------------------------------------------------------------------


@dataclass
class NodeBox:
    label: str
    category: str
    strong: bool | None  # None only for "missing" / "lost", where it's moot
    sublabel: str = ""

    @property
    def shape(self) -> str:
        if self.category == "missing":
            return "missing"
        if self.category == "lost":
            return "rect"
        return "rect" if self.strong else "ellipse"


@dataclass
class Row:
    center: NodeBox
    # Nearest-to-farthest chains of (box, distance-to-previous-box-in-chain).
    left_chain: list[tuple[NodeBox, int | None]] = field(default_factory=list)
    right_chain: list[tuple[NodeBox, int | None]] = field(default_factory=list)


@dataclass
class CartoonJob:
    title: str
    filename: str
    category: str
    ref_rows: list[Row] = field(default_factory=list)
    query_rows: list[Row] = field(default_factory=list)
    query_note: str = ""
    ref_note: str = ""


# --------------------------------------------------------------------------
# Topology: reference/query neighbor lookups + orthology-graph reconstruction.
# --------------------------------------------------------------------------


def _resolve_ref_ids(
    requested: list[str],
    ref_by_uid: dict[str, object],
    ids_by_display_name: dict[str, list[str]],
) -> list[str]:
    resolved: list[str] = []
    for name in requested:
        if name in ref_by_uid:
            resolved.append(name)
            continue
        candidates = ids_by_display_name.get(name, [])
        if len(candidates) == 1:
            resolved.append(candidates[0])
        elif len(candidates) > 1:
            print(
                f"warning: '{name}' matches {len(candidates)} distinct CNEEs "
                f"({', '.join(candidates)}); specify the cnee_uid directly. Skipping.",
                file=sys.stderr,
            )
        else:
            print(f"warning: '{name}' not found in --cnee-bed. Skipping.", file=sys.stderr)
    return resolved


def _build_ref_to_queries(
    ref_partners_by_query_id: dict[str, frozenset[str]],
) -> dict[str, set[str]]:
    ref_to_queries: dict[str, set[str]] = defaultdict(set)
    for query_id, ref_ids in ref_partners_by_query_id.items():
        for ref_id in ref_ids:
            ref_to_queries[ref_id].add(query_id)
    return ref_to_queries


def _discover_component(
    ref_uid: str,
    ref_to_queries: dict[str, set[str]],
    query_to_refs: dict[str, frozenset[str]],
) -> tuple[set[str], set[str]]:
    # BFS the ref<->query bipartite graph implied by the orthology TSV so a
    # many2one/complex partner (e.g. r6 alongside a requested r2) is pulled in
    # automatically, matching how orthology.py grouped them into one component.
    seen_refs = {ref_uid}
    seen_queries: set[str] = set()
    queue: deque[tuple[str, str]] = deque([("ref", ref_uid)])
    while queue:
        kind, node_id = queue.popleft()
        if kind == "ref":
            for query_id in ref_to_queries.get(node_id, ()):
                if query_id not in seen_queries:
                    seen_queries.add(query_id)
                    queue.append(("query", query_id))
        else:
            for other_ref in query_to_refs.get(node_id, ()):
                if other_ref not in seen_refs:
                    seen_refs.add(other_ref)
                    queue.append(("ref", other_ref))
    return seen_refs, seen_queries


def _ref_neighbor_chain(
    ref_index: RefNeighborIndex, cnee_uid: str, side: str, n: int
) -> list[tuple[str, int | None]]:
    # Walk up to n neighbors out from cnee_uid on one side, each paired with
    # the (local) gap to the previous box in the chain -- not the cumulative
    # distance back to cnee_uid -- so every connecting line in the cartoon
    # reads as "gap between these two boxes."
    chrom = ref_index.chrom_by_uid.get(cnee_uid)
    if chrom is None:
        return []
    intervals = ref_index.intervals_by_chrom[chrom]
    position = ref_index.position_by_uid[cnee_uid]
    current = intervals[position]
    chain: list[tuple[str, int | None]] = []
    if side == "left":
        edge = current.start
        index = position - 1
        while index >= 0 and len(chain) < n:
            neighbor = intervals[index]
            chain.append((neighbor.cnee_uid, max(0, edge - neighbor.end)))
            edge = neighbor.start
            index -= 1
    else:
        edge = current.end
        index = position + 1
        while index < len(intervals) and len(chain) < n:
            neighbor = intervals[index]
            chain.append((neighbor.cnee_uid, max(0, neighbor.start - edge)))
            edge = neighbor.end
            index += 1
    return chain


def _query_neighbor_chain(
    query_index: QueryNeighborIndex,
    query_by_id: dict[str, QueryBedRecord],
    chrom: str,
    anchor: int,
    direction: str,
    n: int,
) -> list[tuple[str, int | None]]:
    # Same idea as _ref_neighbor_chain but over the query BED: each step's
    # search anchor is re-based on the previously found neighbor's far edge.
    chain: list[tuple[str, int | None]] = []
    current_anchor = anchor
    for _ in range(n):
        neighbor_id, distance = query_index.neighbor(chrom, current_anchor, direction)
        if neighbor_id is None:
            break
        chain.append((neighbor_id, distance))
        record = query_by_id[neighbor_id]
        current_anchor = record.start if direction == "before" else record.end
    return chain


# --------------------------------------------------------------------------
# Stats loading (filtered to only the CNEEs a job actually needs).
# --------------------------------------------------------------------------


def _read_needed_stats(
    stats_path: str | Path, needed_ref_uids: set[str]
) -> dict[str, list[dict[str, str]]]:
    stats_by_ref: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not needed_ref_uids:
        return stats_by_ref
    with Path(stats_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cnee_uid = row["cnee_uid"]
            if cnee_uid in needed_ref_uids:
                stats_by_ref[cnee_uid].append(row)
    return stats_by_ref


def _ref_is_strong(
    ref_uid: str,
    query_chrom: str,
    query_start: int,
    query_end: int,
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    # True iff any of this ref's own candidate rows lands (overlaps) the given
    # query interval with a solid (strong) match -- an edge property, not a
    # node property, since one2many/many2one rows can differ row to row.
    for row in stats_by_ref.get(ref_uid, []):
        if row["query_chrom"] != query_chrom:
            continue
        start_text, end_text = row["query_aligned_start"], row["query_aligned_end"]
        if start_text == "." or end_text == ".":
            continue
        start, end = int(start_text), int(end_text)
        if start < query_end and end > query_start:
            if stats_row_is_solid(
                row,
                min_sequence_identity=min_sequence_identity,
                max_fraction_gap=max_fraction_gap,
            ):
                return True
    return False


def _query_node_is_strong(
    query_id: str,
    query_by_id: dict[str, QueryBedRecord],
    query_to_refs: dict[str, frozenset[str]],
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> bool:
    record = query_by_id[query_id]
    return any(
        _ref_is_strong(
            ref_id,
            record.chrom,
            record.start,
            record.end,
            stats_by_ref,
            min_sequence_identity,
            max_fraction_gap,
        )
        for ref_id in query_to_refs.get(query_id, ())
    )


# --------------------------------------------------------------------------
# Node builders.
# --------------------------------------------------------------------------


def _ref_node_box(
    ref_uid: str,
    ref_by_uid: dict[str, object],
    category_by_ref_uid: dict[str, str],
    ref_to_queries: dict[str, set[str]],
    query_by_id: dict[str, QueryBedRecord],
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> NodeBox:
    interval = ref_by_uid[ref_uid]
    category = category_by_ref_uid.get(ref_uid, "missing")
    strong: bool | None = None
    if category in ALIGNED_CATEGORIES:
        strong = any(
            _ref_is_strong(
                ref_uid,
                query_by_id[query_id].chrom,
                query_by_id[query_id].start,
                query_by_id[query_id].end,
                stats_by_ref,
                min_sequence_identity,
                max_fraction_gap,
            )
            for query_id in ref_to_queries.get(ref_uid, ())
        )
    elif category == "uncertain_loss":
        strong = False  # never solid by construction (see orthology.py)
    return NodeBox(
        label=interval.cnee_id,
        category=category,
        strong=strong,
        sublabel=f"{interval.chrom}:{interval.start:,}-{interval.end:,}",
    )


def _query_node_box(
    query_id: str,
    query_by_id: dict[str, QueryBedRecord],
    category_by_query_id: dict[str, str],
    query_to_refs: dict[str, frozenset[str]],
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> NodeBox:
    record = query_by_id[query_id]
    category = category_by_query_id.get(query_id, "missing")
    strong = _query_node_is_strong(
        query_id, query_by_id, query_to_refs, stats_by_ref, min_sequence_identity, max_fraction_gap
    )
    return NodeBox(
        label=query_id,
        category=category,
        strong=strong,
        sublabel=f"{record.chrom}:{record.start:,}-{record.end:,}",
    )


def _ref_neighbor_box(
    ref_uid: str | None,
    ref_by_uid: dict[str, object],
    category_by_ref_uid: dict[str, str],
    ref_to_queries: dict[str, set[str]],
    query_by_id: dict[str, QueryBedRecord],
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> NodeBox | None:
    if ref_uid is None:
        return None
    return _ref_node_box(
        ref_uid,
        ref_by_uid,
        category_by_ref_uid,
        ref_to_queries,
        query_by_id,
        stats_by_ref,
        min_sequence_identity,
        max_fraction_gap,
    )


def _query_neighbor_box(
    query_id: str | None,
    query_by_id: dict[str, QueryBedRecord],
    category_by_query_id: dict[str, str],
    query_to_refs: dict[str, frozenset[str]],
    stats_by_ref: dict[str, list[dict[str, str]]],
    min_sequence_identity: float,
    max_fraction_gap: float,
) -> NodeBox | None:
    if query_id is None:
        return None
    return _query_node_box(
        query_id,
        query_by_id,
        category_by_query_id,
        query_to_refs,
        stats_by_ref,
        min_sequence_identity,
        max_fraction_gap,
    )


# --------------------------------------------------------------------------
# Needed-CNEE bookkeeping so the (potentially huge) stats file is scanned for
# only the CNEEs a requested job's diagram will actually show.
# --------------------------------------------------------------------------


@dataclass
class JobTopology:
    category: str
    component_refs: list[str]
    component_queries: list[str]
    ref_left_chain: dict[str, list[tuple[str, int | None]]] = field(default_factory=dict)
    ref_right_chain: dict[str, list[tuple[str, int | None]]] = field(default_factory=dict)
    query_left_chain: dict[str, list[tuple[str, int | None]]] = field(default_factory=dict)
    query_right_chain: dict[str, list[tuple[str, int | None]]] = field(default_factory=dict)


def _plan_aligned_job(
    category: str,
    component_refs: set[str],
    component_queries: set[str],
    ref_by_uid: dict[str, object],
    ref_index: RefNeighborIndex,
    query_by_id: dict[str, QueryBedRecord],
    query_index: QueryNeighborIndex,
    query_to_refs: dict[str, frozenset[str]],
    needed_ref_uids: set[str],
    num_neighbors: int,
) -> JobTopology:
    topo = JobTopology(
        category=category,
        component_refs=sorted(
            component_refs, key=lambda r: (ref_by_uid[r].chrom, ref_by_uid[r].start)
        ),
        component_queries=sorted(
            component_queries, key=lambda q: (query_by_id[q].chrom, query_by_id[q].start)
        ),
    )
    for ref_uid in topo.component_refs:
        needed_ref_uids.add(ref_uid)
        left_chain = _ref_neighbor_chain(ref_index, ref_uid, "left", num_neighbors)
        right_chain = _ref_neighbor_chain(ref_index, ref_uid, "right", num_neighbors)
        topo.ref_left_chain[ref_uid] = left_chain
        topo.ref_right_chain[ref_uid] = right_chain
        needed_ref_uids.update(neighbor_id for neighbor_id, _ in left_chain)
        needed_ref_uids.update(neighbor_id for neighbor_id, _ in right_chain)
    for query_id in topo.component_queries:
        needed_ref_uids.update(query_to_refs.get(query_id, ()))
        record = query_by_id[query_id]
        left_chain = _query_neighbor_chain(
            query_index, query_by_id, record.chrom, record.start, "before", num_neighbors
        )
        right_chain = _query_neighbor_chain(
            query_index, query_by_id, record.chrom, record.end, "after", num_neighbors
        )
        topo.query_left_chain[query_id] = left_chain
        topo.query_right_chain[query_id] = right_chain
        for neighbor_id, _ in left_chain:
            needed_ref_uids.update(query_to_refs.get(neighbor_id, ()))
        for neighbor_id, _ in right_chain:
            needed_ref_uids.update(query_to_refs.get(neighbor_id, ()))
    return topo


def _plan_orphan_job(
    category: str,
    ref_uid: str,
    ref_index: RefNeighborIndex,
    needed_ref_uids: set[str],
    num_neighbors: int,
) -> JobTopology:
    topo = JobTopology(category=category, component_refs=[ref_uid], component_queries=[])
    needed_ref_uids.add(ref_uid)
    left_chain = _ref_neighbor_chain(ref_index, ref_uid, "left", num_neighbors)
    right_chain = _ref_neighbor_chain(ref_index, ref_uid, "right", num_neighbors)
    topo.ref_left_chain[ref_uid] = left_chain
    topo.ref_right_chain[ref_uid] = right_chain
    needed_ref_uids.update(neighbor_id for neighbor_id, _ in left_chain)
    needed_ref_uids.update(neighbor_id for neighbor_id, _ in right_chain)
    return topo


# --------------------------------------------------------------------------
# uncertain_loss / lost query-side fan-out (needs raw stats, and for "lost"
# the chain blocks flanking the gap -- the one place this script reads the
# chain file, exactly as the user specified).
# --------------------------------------------------------------------------


@dataclass
class FanRow:
    label: str
    chrom: str
    anchor_start: int | None
    anchor_end: int | None
    left_chain: list[tuple[str, int | None]]
    right_chain: list[tuple[str, int | None]]
    inferred: bool


def _uncertain_loss_fan(
    ref_uid: str,
    stats_by_ref: dict[str, list[dict[str, str]]],
    query_index: QueryNeighborIndex,
    query_by_id: dict[str, QueryBedRecord],
    num_neighbors: int,
) -> list[FanRow]:
    rows = [
        row
        for row in stats_by_ref.get(ref_uid, [])
        if int(row["aligned_bases"]) > 0
        and row["query_aligned_start"] != "."
        and row["query_aligned_end"] != "."
    ]
    rows.sort(key=lambda row: (row["query_chrom"], int(row["query_aligned_start"])))
    fans: list[FanRow] = []
    for row in rows:
        chrom = row["query_chrom"]
        start, end = int(row["query_aligned_start"]), int(row["query_aligned_end"])
        left_chain = _query_neighbor_chain(query_index, query_by_id, chrom, start, "before", num_neighbors)
        right_chain = _query_neighbor_chain(query_index, query_by_id, chrom, end, "after", num_neighbors)
        fans.append(
            FanRow(
                label=row.get("query_projection_name", row["cnee_id"]),
                chrom=chrom,
                anchor_start=start,
                anchor_end=end,
                left_chain=left_chain,
                right_chain=right_chain,
                inferred=False,
            )
        )
    return fans


def _lost_fan(
    ref_uid: str,
    stats_by_ref: dict[str, list[dict[str, str]]],
    query_index: QueryNeighborIndex,
    query_by_id: dict[str, QueryBedRecord],
    chain_records_by_id: dict[str, object],
    chain_handle,
    block_cache: dict[str, BlockIndex],
    num_neighbors: int,
) -> list[FanRow]:
    rows = sorted(stats_by_ref.get(ref_uid, []), key=lambda row: row["chain_id"])
    fans: list[FanRow] = []
    for row in rows:
        chain_id = row["chain_id"]
        record = chain_records_by_id[chain_id]
        if chain_id not in block_cache:
            chain_handle.seek(record.start_byte)
            chain_text = chain_handle.read(record.byte_length).decode("utf-8")
            block_cache[chain_id] = BlockIndex.from_blocks(iter_chain_blocks(chain_text))
        blocks = block_cache[chain_id]
        anchors = _query_flank_anchors(row, record, blocks)
        left_anchor, left_dir = anchors["left"]
        right_anchor, right_dir = anchors["right"]
        left_chain: list[tuple[str, int | None]] = []
        right_chain: list[tuple[str, int | None]] = []
        if left_anchor is not None:
            left_chain = _query_neighbor_chain(
                query_index, query_by_id, record.query_chrom, left_anchor, left_dir, num_neighbors
            )
        if right_anchor is not None:
            right_chain = _query_neighbor_chain(
                query_index, query_by_id, record.query_chrom, right_anchor, right_dir, num_neighbors
            )
        points = [p for p in (left_anchor, right_anchor) if p is not None]
        anchor_start = min(points) if points else None
        anchor_end = max(points) if points else None
        fans.append(
            FanRow(
                label=f"{row['cnee_id']} (inferred)",
                chrom=record.query_chrom,
                anchor_start=anchor_start,
                anchor_end=anchor_end,
                left_chain=left_chain,
                right_chain=right_chain,
                inferred=True,
            )
        )
    return fans


# --------------------------------------------------------------------------
# Drawing.
# --------------------------------------------------------------------------

BOX_W, BOX_H = 1.5, 0.85
SLOT_DX = 2.3
ROW_DY = 1.7
SECTION_GAP = 1.3
LEGEND_ROW_DY = 0.5
LEGEND_TOP_GAP = 0.9
TITLE_HEIGHT = 0.55
SUBTITLE_HEIGHT = 0.45
TOP_MARGIN = 0.3
BOTTOM_MARGIN = 0.3
TITLE_WRAP_CHARS = 78
NOTE_HEIGHT = 0.5
CHARS_PER_UNIT_LABEL = 15 / 2.3
CHARS_PER_UNIT_SUBLABEL = 24 / 2.3
DEFAULT_NUM_NEIGHBORS = 5


def _format_bp(distance: int | None) -> str:
    if distance is None:
        return "(edge of contig)"
    if distance >= 1_000_000:
        return f"{distance / 1_000_000:.2f} Mb"
    if distance >= 1_000:
        return f"{distance / 1_000:.1f} kb"
    return f"{distance:,} bp"


def _fit_fontsize(text: str, base_fontsize: float, base_chars: int, min_fontsize: float) -> float:
    # Ellipses/rectangles have a fixed data-unit width, so long labels (a
    # 9-ref many2one component, a raw chain projection name, ...) need a
    # smaller font rather than overflowing past the box edge.
    if len(text) <= base_chars:
        return base_fontsize
    return max(min_fontsize, base_fontsize * base_chars / len(text))


def _draw_box(ax, x: float, y: float, box: NodeBox | None) -> None:
    if box is None:
        return
    color = CATEGORY_COLORS.get(box.category, "#FFFFFF")
    shape = box.shape
    # Ellipses taper, so they have less usable horizontal text room than a
    # rectangle of the same nominal width.
    width_factor = 0.72 if shape == "ellipse" else 1.0
    if shape == "ellipse":
        patch = Ellipse(
            (x, y), BOX_W, BOX_H, facecolor=color, edgecolor="black", linewidth=1.2, zorder=2
        )
    elif shape == "missing":
        patch = Rectangle(
            (x - BOX_W / 2, y - BOX_H / 2),
            BOX_W,
            BOX_H,
            facecolor=color,
            edgecolor="#666666",
            linewidth=1.2,
            linestyle="--",
            hatch="//",
            zorder=2,
        )
    else:
        patch = Rectangle(
            (x - BOX_W / 2, y - BOX_H / 2),
            BOX_W,
            BOX_H,
            facecolor=color,
            edgecolor="black",
            linewidth=1.2,
            zorder=2,
        )
    ax.add_patch(patch)
    text_color = "black" if box.category in {"uncertain_loss", "missing"} else "white"
    label_fontsize = _fit_fontsize(
        box.label, 9, base_chars=int(CHARS_PER_UNIT_LABEL * BOX_W * width_factor), min_fontsize=5.0
    )
    ax.text(
        x,
        y + 0.12,
        box.label,
        ha="center",
        va="center",
        fontsize=label_fontsize,
        fontweight="bold",
        color=text_color,
        zorder=3,
    )
    if box.sublabel:
        sublabel_fontsize = _fit_fontsize(
            box.sublabel,
            6.3,
            base_chars=int(CHARS_PER_UNIT_SUBLABEL * BOX_W * width_factor),
            min_fontsize=4.2,
        )
        ax.text(
            x,
            y - 0.22,
            box.sublabel,
            ha="center",
            va="center",
            fontsize=sublabel_fontsize,
            color=text_color,
            zorder=3,
        )


def _draw_connector(ax, x_from: float, x_to: float, y: float, distance: int | None) -> None:
    ax.plot(
        [x_from + BOX_W / 2, x_to - BOX_W / 2],
        [y, y],
        color="#1f2937",
        linewidth=1.2,
        zorder=1,
    )
    ax.text(
        (x_from + x_to) / 2,
        y + 0.16,
        _format_bp(distance),
        ha="center",
        va="bottom",
        fontsize=7.0,
        color="#1f2937",
        zorder=3,
    )


def _draw_row(ax, y: float, row: Row) -> None:
    _draw_box(ax, 0.0, y, row.center)
    prev_x = 0.0
    for box, distance in row.left_chain:
        next_x = prev_x - SLOT_DX
        _draw_box(ax, next_x, y, box)
        _draw_connector(ax, next_x, prev_x, y, distance)
        prev_x = next_x
    prev_x = 0.0
    for box, distance in row.right_chain:
        next_x = prev_x + SLOT_DX
        _draw_box(ax, next_x, y, box)
        _draw_connector(ax, prev_x, next_x, y, distance)
        prev_x = next_x


def _draw_legend(ax, y_top: float, x_left: float, x_right: float) -> float:
    swatch_w, swatch_h = 0.5, 0.35
    x0 = x_left
    y = y_top
    for category in CATEGORY_ORDER:
        color = CATEGORY_COLORS[category]
        style = {"facecolor": color, "edgecolor": "black", "linewidth": 1.0}
        if category == "missing":
            style.update(edgecolor="#666666", linestyle="--", hatch="//")
        ax.add_patch(Rectangle((x0, y - swatch_h / 2), swatch_w, swatch_h, **style))
        ax.text(x0 + swatch_w + 0.15, y, CATEGORY_LABELS[category], va="center", fontsize=8)
        y -= LEGEND_ROW_DY

    shape_x0 = x_right - 2.2
    y = y_top
    ax.add_patch(
        Rectangle(
            (shape_x0, y - swatch_h / 2),
            swatch_w,
            swatch_h,
            facecolor="white",
            edgecolor="black",
            linewidth=1.0,
        )
    )
    ax.text(shape_x0 + swatch_w + 0.15, y, "strong", va="center", fontsize=8)
    y -= LEGEND_ROW_DY
    ax.add_patch(
        Ellipse(
            (shape_x0 + swatch_w / 2, y), swatch_w, swatch_h, facecolor="white", edgecolor="black", linewidth=1.0
        )
    )
    ax.text(shape_x0 + swatch_w + 0.15, y, "weak", va="center", fontsize=8)
    return y - LEGEND_ROW_DY


def render_job(
    job: CartoonJob, out_dir: Path, dpi: int, thresholds_note: str, num_neighbors: int
) -> Path:
    import textwrap

    n_ref_rows = max(len(job.ref_rows), 1)
    n_query_rows = max(len(job.query_rows), 1)
    title_lines = textwrap.wrap(job.title, TITLE_WRAP_CHARS) or [job.title]
    ref_note_h = NOTE_HEIGHT if job.ref_note else 0.0
    # A truncation note on the query side only needs its own slot when rows
    # are also being drawn -- the "no rows at all" case already reserves a
    # slot for job.query_note via the n_query_rows=1 floor above.
    query_note_h = NOTE_HEIGHT if (job.query_note and job.query_rows) else 0.0

    # Compute every vertical extent up front so the axes' ylim always matches
    # what's actually drawn -- a fixed guess here caused the legend (and
    # title/subtitle) to silently clip off in earlier iterations.
    title_block_h = TITLE_HEIGHT * len(title_lines)
    ref_block_h = n_ref_rows * ROW_DY + ref_note_h
    query_block_h = n_query_rows * ROW_DY + query_note_h
    legend_h = len(CATEGORY_ORDER) * LEGEND_ROW_DY

    total_h = (
        TOP_MARGIN
        + title_block_h
        + SUBTITLE_HEIGHT
        + ref_block_h
        + SECTION_GAP
        + query_block_h
        + LEGEND_TOP_GAP
        + legend_h
        + BOTTOM_MARGIN
    )

    # Horizontal extent scales with how many neighbors are shown per side --
    # center slot (index 0) plus num_neighbors slots on each side.
    x_left = -(num_neighbors * SLOT_DX) - BOX_W / 2 - 1.0
    x_right = (num_neighbors * SLOT_DX) + BOX_W / 2 + 1.0
    section_label_x = x_left - 1.0
    fig_w = max(10.0, (2 * num_neighbors + 1) * 1.55)

    fig, ax = plt.subplots(figsize=(fig_w, max(total_h * 0.62, 4.0)))
    ax.set_xlim(section_label_x - 1.4, x_right + 0.6)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    x_mid = 0.0
    y = total_h - TOP_MARGIN
    for i, line in enumerate(title_lines):
        ax.text(
            x_mid,
            y - TITLE_HEIGHT * (i + 0.5),
            line,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
    y -= title_block_h
    ax.text(
        x_mid,
        y - SUBTITLE_HEIGHT / 2,
        thresholds_note,
        ha="center",
        va="center",
        fontsize=8,
        color="#555555",
    )
    y -= SUBTITLE_HEIGHT

    # Ref section: each row is drawn centered in its own ROW_DY-tall slot,
    # walking top-down from y.
    ax.text(section_label_x, y - ref_block_h / 2, "Ref:", fontsize=11, fontweight="bold", va="center")
    row_y = y - ROW_DY / 2
    for row in job.ref_rows:
        _draw_row(ax, row_y, row)
        row_y -= ROW_DY
    if job.ref_note:
        ax.text(
            x_mid,
            row_y + ROW_DY / 2 - NOTE_HEIGHT / 2,
            job.ref_note,
            ha="center",
            va="center",
            fontsize=8.5,
            style="italic",
            color="#666666",
        )
    y -= ref_block_h + SECTION_GAP

    # Query section.
    ax.text(section_label_x, y - query_block_h / 2, "Query:", fontsize=11, fontweight="bold", va="center")
    if job.query_rows:
        row_y = y - ROW_DY / 2
        for row in job.query_rows:
            _draw_row(ax, row_y, row)
            row_y -= ROW_DY
        if job.query_note:
            ax.text(
                x_mid,
                row_y + ROW_DY / 2 - NOTE_HEIGHT / 2,
                job.query_note,
                ha="center",
                va="center",
                fontsize=8.5,
                style="italic",
                color="#666666",
            )
    else:
        ax.text(
            x_mid,
            y - query_block_h / 2,
            job.query_note or "no query data",
            ha="center",
            va="center",
            fontsize=9,
            style="italic",
            color="#666666",
        )
    y -= query_block_h + LEGEND_TOP_GAP

    _draw_legend(ax, y, x_left, x_right)

    out_path = out_dir / job.filename
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------
# Main.
# --------------------------------------------------------------------------


MAX_ROWS_PER_SECTION = 25
MAX_TITLE_LABELS = 8
MAX_FILENAME_LABELS = 4


def _clean_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def _title_label_list(labels: list[str]) -> str:
    if len(labels) <= MAX_TITLE_LABELS:
        return ", ".join(labels)
    shown = ", ".join(labels[:MAX_TITLE_LABELS])
    return f"{shown}, +{len(labels) - MAX_TITLE_LABELS} more"


def _safe_filename(labels: list[str]) -> str:
    # A many2one/complex component can chain into hundreds of refs in
    # repetitive regions; joining every id verbatim can exceed the OS's
    # filename length limit, so fall back to a truncated name + content hash.
    joined = _clean_filename("_".join(labels))
    if len(joined) <= 140:
        return joined + ".png"
    digest = hashlib.md5(",".join(labels).encode("utf-8")).hexdigest()[:8]
    head = _clean_filename("_".join(labels[:MAX_FILENAME_LABELS]))
    return f"{head}_plus{len(labels) - MAX_FILENAME_LABELS}more_{digest}.png"


def build_jobs(
    requested_ref_uids: list[str],
    category_by_ref_uid: dict[str, str],
    ref_to_queries: dict[str, set[str]],
    query_to_refs: dict[str, frozenset[str]],
    skip_categories: set[str] = frozenset(),
) -> "OrderedDict[frozenset, tuple[str, set[str], set[str]]]":
    from collections import OrderedDict

    jobs: "OrderedDict[frozenset, tuple[str, set[str], set[str]]]" = OrderedDict()
    for ref_uid in requested_ref_uids:
        category = category_by_ref_uid.get(ref_uid, "missing")
        if category in skip_categories:
            print(
                f"skipping {ref_uid}: category '{category}' is excluded "
                "(pass --include-complex to include it)",
                file=sys.stderr,
            )
            continue
        if category in ALIGNED_CATEGORIES:
            component_refs, component_queries = _discover_component(
                ref_uid, ref_to_queries, query_to_refs
            )
        else:
            component_refs, component_queries = {ref_uid}, set()
        key = frozenset(component_refs)
        jobs.setdefault(key, (category, component_refs, component_queries))
    return jobs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cnee_ids", nargs="*", help="Reference CNEE cnee_id/cnee_uid values to plot")
    parser.add_argument("--cnee-id-file", type=Path, default=None, help="File with one CNEE id per line")
    parser.add_argument("--cnee-bed", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path, help="All-mode graph-assigned query_cnees.bed")
    parser.add_argument("--stats", required=True, type=Path, help="overlap_stats.tsv (all-mode, raw)")
    parser.add_argument("--orthology", required=True, type=Path, help="overlap_stats.orthology.tsv (all-mode)")
    parser.add_argument(
        "--chain-index",
        type=Path,
        default=None,
        help="chain_index.json; required only if a requested/neighboring CNEE is 'lost'",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-sequence-identity", type=float, default=0.7)
    parser.add_argument("--max-fraction-gap", type=float, default=0.1)
    parser.add_argument(
        "--num-neighbors",
        type=int,
        default=DEFAULT_NUM_NEIGHBORS,
        help=f"Ref/query neighbors to show on each side of the candidate [default: {DEFAULT_NUM_NEIGHBORS}]",
    )
    parser.add_argument(
        "--include-complex",
        action="store_true",
        help="Also generate cartoons for requested CNEEs whose category is complex/many2many "
        "(skipped by default -- these components tend to be too cluttered to be useful)",
    )
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args(argv)

    requested_names = list(args.cnee_ids)
    if args.cnee_id_file is not None:
        requested_names.extend(
            line.strip()
            for line in args.cnee_id_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    if not requested_names:
        parser.error("no CNEE ids given (positional args or --cnee-id-file)")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    ref_intervals = parse_bed(args.cnee_bed)
    ref_by_uid = {iv.cnee_uid: iv for iv in ref_intervals}
    ids_by_display_name: dict[str, list[str]] = defaultdict(list)
    for iv in ref_intervals:
        ids_by_display_name[iv.cnee_id].append(iv.cnee_uid)
    ref_index = RefNeighborIndex.from_bed(args.cnee_bed)

    orth = OrthologyIndex.from_path(args.orthology)
    ref_to_queries = _build_ref_to_queries(orth.ref_partners_by_query_id)

    query_by_chrom = _read_query_bed(args.query_bed)
    query_by_id: dict[str, QueryBedRecord] = {
        rec.query_id: rec for records in query_by_chrom.values() for rec in records
    }
    query_index = QueryNeighborIndex.from_records(query_by_chrom)

    resolved_ref_uids = _resolve_ref_ids(requested_names, ref_by_uid, ids_by_display_name)
    if not resolved_ref_uids:
        print("error: none of the requested CNEE ids could be resolved", file=sys.stderr)
        raise SystemExit(1)

    skip_categories = set() if args.include_complex else {"complex/many2many"}
    jobs = build_jobs(
        resolved_ref_uids,
        orth.category_by_ref_uid,
        ref_to_queries,
        orth.ref_partners_by_query_id,
        skip_categories=skip_categories,
    )
    if not jobs:
        print("error: nothing left to plot after category filtering", file=sys.stderr)
        raise SystemExit(1)

    needed_ref_uids: set[str] = set()
    topologies: dict[frozenset, JobTopology] = {}
    for key, (category, component_refs, component_queries) in jobs.items():
        if category in ALIGNED_CATEGORIES:
            topologies[key] = _plan_aligned_job(
                category,
                component_refs,
                component_queries,
                ref_by_uid,
                ref_index,
                query_by_id,
                query_index,
                orth.ref_partners_by_query_id,
                needed_ref_uids,
                args.num_neighbors,
            )
        else:
            (ref_uid,) = component_refs
            topologies[key] = _plan_orphan_job(
                category, ref_uid, ref_index, needed_ref_uids, args.num_neighbors
            )

    stats_by_ref = _read_needed_stats(args.stats, needed_ref_uids)

    # Second pass: uncertain_loss/lost fan-out rows need their own stats (just
    # scanned above) plus, for "lost" only, the chain blocks flanking the gap.
    needs_chain = any(topo.category == "lost" for topo in topologies.values())
    chain_records_by_id: dict[str, object] = {}
    chain_handle = None
    block_cache: dict[str, BlockIndex] = {}
    if needs_chain:
        if args.chain_index is None:
            print(
                "error: at least one requested/neighboring CNEE is 'lost', which needs "
                "--chain-index to interpolate its query position",
                file=sys.stderr,
            )
            raise SystemExit(1)
        payload, records = load_index(args.chain_index)
        chain_records_by_id = {record.chain_id: record for record in records}
        chain_handle = Path(payload["indexed_chain_path"]).open("rb")

    fan_by_ref: dict[str, list[FanRow]] = {}
    additional_needed: set[str] = set()
    for topo in topologies.values():
        if topo.category not in {"uncertain_loss", "lost"}:
            continue
        (ref_uid,) = topo.component_refs
        if topo.category == "uncertain_loss":
            fans = _uncertain_loss_fan(
                ref_uid, stats_by_ref, query_index, query_by_id, args.num_neighbors
            )
        else:
            fans = _lost_fan(
                ref_uid,
                stats_by_ref,
                query_index,
                query_by_id,
                chain_records_by_id,
                chain_handle,
                block_cache,
                args.num_neighbors,
            )
        fan_by_ref[ref_uid] = fans
        for fan in fans:
            for neighbor_id, _ in fan.left_chain:
                additional_needed.update(orth.ref_partners_by_query_id.get(neighbor_id, ()))
            for neighbor_id, _ in fan.right_chain:
                additional_needed.update(orth.ref_partners_by_query_id.get(neighbor_id, ()))

    additional_needed -= needed_ref_uids
    if additional_needed:
        stats_by_ref.update(_read_needed_stats(args.stats, additional_needed))

    if chain_handle is not None:
        chain_handle.close()

    thresholds_note = (
        f"solid-match thresholds: sequence identity > {args.min_sequence_identity}, "
        f"gap fraction < {args.max_fraction_gap}"
    )

    written: list[Path] = []
    for key, (category, component_refs, component_queries) in jobs.items():
        topo = topologies[key]

        def ref_box(ref_uid: str | None) -> NodeBox | None:
            return _ref_neighbor_box(
                ref_uid,
                ref_by_uid,
                orth.category_by_ref_uid,
                ref_to_queries,
                query_by_id,
                stats_by_ref,
                args.min_sequence_identity,
                args.max_fraction_gap,
            )

        def query_box(query_id: str | None) -> NodeBox | None:
            return _query_neighbor_box(
                query_id,
                query_by_id,
                orth.category_by_query_id,
                orth.ref_partners_by_query_id,
                stats_by_ref,
                args.min_sequence_identity,
                args.max_fraction_gap,
            )

        ref_note = ""
        ref_uids_to_draw = topo.component_refs
        if len(ref_uids_to_draw) > MAX_ROWS_PER_SECTION:
            dropped = len(ref_uids_to_draw) - MAX_ROWS_PER_SECTION
            print(
                f"warning: component for {topo.component_refs[0]} has "
                f"{len(ref_uids_to_draw)} ref CNEEs; showing only the first "
                f"{MAX_ROWS_PER_SECTION} (dropping {dropped})",
                file=sys.stderr,
            )
            ref_uids_to_draw = ref_uids_to_draw[:MAX_ROWS_PER_SECTION]
            ref_note = f"+{dropped} more ref CNEE(s) in this component not shown"

        ref_rows = [
            Row(
                center=ref_box(ref_uid),
                left_chain=[
                    (ref_box(neighbor_id), dist)
                    for neighbor_id, dist in topo.ref_left_chain.get(ref_uid, [])
                ],
                right_chain=[
                    (ref_box(neighbor_id), dist)
                    for neighbor_id, dist in topo.ref_right_chain.get(ref_uid, [])
                ],
            )
            for ref_uid in ref_uids_to_draw
        ]

        query_rows: list[Row] = []
        query_note = ""
        if category in ALIGNED_CATEGORIES:
            query_ids_to_draw = topo.component_queries
            if len(query_ids_to_draw) > MAX_ROWS_PER_SECTION:
                dropped = len(query_ids_to_draw) - MAX_ROWS_PER_SECTION
                print(
                    f"warning: component for {topo.component_refs[0]} has "
                    f"{len(query_ids_to_draw)} query CNEEs; showing only the first "
                    f"{MAX_ROWS_PER_SECTION} (dropping {dropped})",
                    file=sys.stderr,
                )
                query_ids_to_draw = query_ids_to_draw[:MAX_ROWS_PER_SECTION]
                query_note = f"+{dropped} more query CNEE(s) in this component not shown"
            for query_id in query_ids_to_draw:
                query_rows.append(
                    Row(
                        center=query_box(query_id),
                        left_chain=[
                            (query_box(neighbor_id), dist)
                            for neighbor_id, dist in topo.query_left_chain.get(query_id, [])
                        ],
                        right_chain=[
                            (query_box(neighbor_id), dist)
                            for neighbor_id, dist in topo.query_right_chain.get(query_id, [])
                        ],
                    )
                )
        elif category == "missing":
            query_note = "missing: no candidate chain overlaps this CNEE at all"
        else:
            (ref_uid,) = component_refs
            fans = fan_by_ref.get(ref_uid, [])
            if not fans:
                query_note = f"{category}: no candidate alignment found"
            for fan in fans:
                center = NodeBox(
                    label=fan.label,
                    category=category,
                    strong=False if category == "uncertain_loss" else None,
                    sublabel=(
                        f"{fan.chrom}:{fan.anchor_start:,}-{fan.anchor_end:,}"
                        if fan.anchor_start is not None and fan.anchor_end is not None
                        else fan.chrom
                    ),
                )
                query_rows.append(
                    Row(
                        center=center,
                        left_chain=[
                            (query_box(neighbor_id), dist) for neighbor_id, dist in fan.left_chain
                        ],
                        right_chain=[
                            (query_box(neighbor_id), dist) for neighbor_id, dist in fan.right_chain
                        ],
                    )
                )

        ref_labels = [ref_by_uid[r].cnee_id for r in topo.component_refs]
        title = f"CNEE {_title_label_list(ref_labels)}: {CATEGORY_LABELS[category]}"
        filename = _safe_filename(ref_labels)
        job = CartoonJob(
            title=title,
            filename=filename,
            category=category,
            ref_rows=ref_rows,
            query_rows=query_rows,
            query_note=query_note,
            ref_note=ref_note,
        )
        out_path = render_job(job, args.out_dir, args.dpi, thresholds_note, args.num_neighbors)
        written.append(out_path)
        print(f"Wrote {out_path}")

    print(f"Wrote {len(written)} PNG(s) for {len(resolved_ref_uids)} requested CNEE id(s) -> {args.out_dir}")


if __name__ == "__main__":
    main()
