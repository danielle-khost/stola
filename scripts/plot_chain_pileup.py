#!/usr/bin/env python
"""IGV-style chain pileups over multi-projection query loci, before (top
panel) and after (bottom panel) multi-projection filtering, from a completed
`stola` run.

    python scripts/plot_chain_pileup.py RUN_DIR --ref-genome ref.fa --query-genome query.2bit \\
        [--random 12 | --locus chrQ:1000-2000 ...]

A locus is the merged query span of the multi-projection (one2many, many2one,
complex) candidate projections landing there -- the same loci the pipeline
filters. Each figure is centered on it and extends --flank bases each side.
Every chain with aligned blocks in that view is one row: aligned blocks are
boxes, chain gaps dotted lines (running to the panel edge where the chain
continues beyond the view), and each candidate projection's own span is
outlined in black. Rows are colored by what filtering did:
  - kept: carries a candidate projection that survived
  - depth / no_flank_support / low_score: carries a candidate projection
    pruned for that reason (see pruned_query_projections.tsv)
  - competitor: scored by the filter at this locus, but not a candidate
  - other: aligned in the view but not scored (only in the flanks, or
    outside the soft depth filter's top chains)
The note right of each row gives the chain's local score components at the
locus: identity over the locus, aligned fraction over the locus +/- flanks,
and relative score R.

Scores and verdicts are recomputed with the pipeline's own score_locus /
decide_locus, then checked against the run's pruned_query_projections.tsv;
pass the same filter settings the run used (defaults match the pipeline's).
Rows are capped at --max-rows per panel.

--random N picks N loci at random (seeded), cycling across the initial
categories present so one2many, many2one and complex loci all appear.
Figures go to --out-dir (default RUN_DIR/chain_pileups/) with an index.tsv.

Requires matplotlib (`pip install matplotlib`, or `pip install .[viz]` from
the repo root) in addition to the stola package itself.
"""

from __future__ import annotations

import argparse
import csv
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from stola.index import load_index  # noqa: E402
from stola.orthology import (  # noqa: E402
    Projection,
    _iter_overlapping_projection_groups,
    _read_projection_evidence,
)
from stola.sequence_compare import SequenceComparator  # noqa: E402
from stola.spurious_filter import (  # noqa: E402
    DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD,
    DEFAULT_FLANK_WINDOW,
    DEFAULT_RELATIVE_SCORE_THRESHOLD,
    DEFAULT_SOFT_DEPTH_LIMIT,
    ChainEnvelopeCache,
    build_query_envelopes,
    decide_locus,
    score_locus,
)

MULTI_CATEGORIES = ("one2many", "many2one", "complex/many2many")
# Envelope-overlapping chains parsed per view when looking for chains with
# aligned blocks in it, best chain score first.
MAX_CHAINS_PARSED = 500

STATUS_COLORS = {
    "kept": "#2a78d6",
    "low_score": "#eb6834",
    "no_flank_support": "#1baf7a",
    "depth": "#4a3aa7",
    "competitor": "#6b7280",
    "other": "#c4c9d0",
}
STATUS_LABELS = {
    "kept": "kept projection",
    "low_score": "pruned: low_score",
    "no_flank_support": "pruned: no_flank_support",
    "depth": "pruned: depth",
    "competitor": "scored competitor (not a candidate)",
    "other": "other aligned chain",
}
STATUS_ORDER = list(STATUS_COLORS)


def parse_region(text: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"(.+):([\d,]+)-([\d,]+)", text)
    if not match:
        raise argparse.ArgumentTypeError(f"expected chrom:start-end, got {text!r}")
    chrom, start, end = match.groups()
    return chrom, int(start.replace(",", "")), int(end.replace(",", ""))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chains_aligned_in_view(query_chrom, view_start, view_end, envelopes, cache):
    """(chain_ids with aligned bases in the view, best chain score first;
    number of chains whose envelope overlaps the view)."""
    envelope = envelopes[query_chrom]
    mask = (envelope["starts"] < view_end) & (envelope["ends"] > view_start)
    ids = envelope["chain_ids"][mask]
    scores = envelope["scores"][mask]
    order = sorted(range(len(ids)), key=lambda i: -scores[i])[:MAX_CHAINS_PARSED]
    aligned = []
    for i in order:
        chain_id = str(ids[i])
        _, blocks = cache.get(chain_id)
        if blocks.aligned_bases(view_start, view_end) > 0:
            aligned.append(chain_id)
    return aligned, int(ids.size)


def draw_panel(ax, rows, view_start, view_end, window, cache, title, max_rows):
    """rows: list of (chain_id, status, label, note, projection spans)."""
    import matplotlib.patches as mpatches

    shown = rows[:max_rows]
    ax.axvspan(window[0], window[1], color="#f6d860", alpha=0.5, zorder=0)
    ax.annotate("", xy=((window[0] + window[1]) / 2, 1.0), xycoords=("data", "axes fraction"),
                xytext=(0, 9), textcoords="offset points",
                arrowprops={"arrowstyle": "-|>", "color": "#0b0b0b", "lw": 1.2})
    for y, (chain_id, status, _label, note, spans) in enumerate(shown):
        _, blocks = cache.get(chain_id)
        color = STATUS_COLORS[status]
        segments = [
            (max(b.query_start, view_start), min(b.query_end, view_end))
            for b in blocks.overlapping(view_start, view_end)
        ]
        # Dotted links between consecutive blocks, plus to the panel edge
        # where the chain continues beyond the view.
        intervals = list(segments)
        if blocks.query_starts[0] < view_start:
            intervals.insert(0, (view_start, view_start))
        if blocks.query_ends[-1] > view_end:
            intervals.append((view_end, view_end))
        for (_, gap_start), (gap_end, _) in zip(intervals, intervals[1:]):
            if gap_end > gap_start:
                ax.plot([gap_start, gap_end], [y, y], linestyle=(0, (1, 2)), color=color, linewidth=1, zorder=1)
        for seg_start, seg_end in segments:
            ax.add_patch(mpatches.Rectangle(
                (seg_start, y - 0.32), seg_end - seg_start, 0.64, facecolor=color, edgecolor="none", zorder=2,
            ))
        for span_start, span_end in spans:
            ax.add_patch(mpatches.Rectangle(
                (span_start, y - 0.42), span_end - span_start, 0.84, facecolor="none",
                edgecolor="#0b0b0b", linewidth=1.2, zorder=3,
            ))
        ax.text(view_end, y, f"  {note}", va="center", ha="left", fontsize=7, color="#52514e", clip_on=False)

    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([label for _, _, label, _, _ in shown], fontsize=7)
    ax.set_ylim(len(shown) - 0.5 if shown else 0.5, -0.6)
    ax.set_xlim(view_start, view_end)
    hidden = len(rows) - len(shown)
    ax.set_title(title + (f"  (+{hidden} more not shown)" if hidden > 0 else ""), fontsize=10, loc="left")
    if not shown:
        ax.text(0.5, 0.5, "no projections retained at this locus", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="#52514e")
    ax.tick_params(axis="x", labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x", color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)


def plot_locus(locus, locus_index, initial_by_ref, final_by_ref, pruned_names, records_by_id, envelopes, cache,
               comparator, args, out_dir):
    import matplotlib.pyplot as plt

    scores = score_locus(
        locus, records_by_id, envelopes, cache, comparator.ref, comparator.query,
        args.flank, args.soft_depth_limit,
    )
    verdicts = decide_locus(locus, scores, args.relative_score_threshold, args.flank_any_coverage_threshold)
    mismatches = [
        p.projection_name for p, v in zip(locus, verdicts) if v.accepted == (p.projection_name in pruned_names)
    ]

    chrom = locus[0].query_chrom
    window = (scores.window_start, scores.window_end)
    view_start, view_end = scores.slice_start, scores.slice_end
    local = {s.chain_id: scores.local_score(s) for s in scores.slices}
    s_max = max(local.values(), default=0.0)
    slice_len = scores.slice_end - scores.slice_start
    stats_by_chain = {
        s.chain_id: (
            f"id {s.window_matches / max(1, s.window_aligned_bases):.2f} · af {s.slice_aligned_bases / slice_len:.2f}"
            + (f" · R {local[s.chain_id] / s_max:.2f}" if s_max > 0 else "")
        )
        for s in scores.slices
    }

    projections_by_chain: dict[str, list[tuple[Projection, object]]] = defaultdict(list)
    for p, v in zip(locus, verdicts):
        projections_by_chain[p.chain_id].append((p, v))

    def chain_status(chain_id: str) -> str:
        entries = projections_by_chain.get(chain_id)
        if entries:
            if any(v.accepted for _, v in entries):
                return "kept"
            return entries[0][1].outcome
        return "competitor" if chain_id in local else "other"

    def row(chain_id: str):
        record, blocks = cache.get(chain_id)
        in_view = blocks.overlapping(view_start, view_end)
        ref_pos = f"{record.ref_chrom}:{in_view[0].ref_start:,}" if in_view else record.ref_chrom
        label = f"{chain_id} · {ref_pos} ({record.query_strand})"
        entries = projections_by_chain.get(chain_id, [])
        refs = sorted({p.ref_id for p, _ in entries})
        note = stats_by_chain.get(chain_id, "")
        if refs:
            # Names are cut at "__" (e.g. "__nearest_<gene>" suffixes); index.tsv has them in full.
            short = [r.split("__")[0] for r in refs]
            shown_refs = ", ".join(short[:2]) + (f" +{len(short) - 2}" if len(short) > 2 else "")
            note = f"{note}  [{shown_refs}]".strip()
        spans = [(p.query_start, p.query_end) for p, _ in entries]
        return chain_id, chain_status(chain_id), label, note, spans

    aligned, envelope_count = chains_aligned_in_view(chrom, view_start, view_end, envelopes, cache)
    for chain_id in projections_by_chain:  # candidate chains always appear
        if chain_id not in aligned:
            aligned.append(chain_id)
    score_rank = {cid: rank for rank, cid in enumerate(aligned)}
    before = sorted(
        (row(cid) for cid in aligned),
        key=lambda r: (STATUS_ORDER.index(r[1]), score_rank[r[0]]),
    )
    after = [r for r in before if r[1] == "kept"]

    refs = sorted({p.ref_id for p in locus})
    initial = sorted({initial_by_ref.get(r, "?") for r in refs})
    finals = sorted({final_by_ref.get(r, "?") for r in refs})
    n_rows = min(len(before), args.max_rows) + max(1, min(len(after), args.max_rows))
    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, sharex=True, figsize=(16, 2.4 + 0.24 * n_rows),
        gridspec_kw={"height_ratios": [max(1, min(len(before), args.max_rows)), max(1, min(len(after), args.max_rows))]},
    )
    draw_panel(
        ax_top, before, view_start, view_end, window, cache,
        f"Before filtering: {len(before)} chains aligned in view"
        + (f" ({envelope_count:,} chain spans overlap it)" if envelope_count > len(before) else ""),
        args.max_rows,
    )
    draw_panel(ax_bottom, after, view_start, view_end, window, cache,
               f"After filtering: {len(after)} chains with retained projections", args.max_rows)
    # Ticks every 2.5 kb, labelled as kb offsets from the locus center.
    center = (window[0] + window[1]) / 2
    step = 2500
    ticks = [center + k * step for k in range(-int((center - view_start) // step), int((view_end - center) // step) + 1)]
    ax_bottom.set_xticks(ticks)
    ax_bottom.set_xticklabels([f"{(x - center) / 1000:+g} kb" if x != center else "0" for x in ticks])
    ax_bottom.set_xlabel(
        f"{chrom} position relative to the CNEE locus (arrow, shaded); view = locus ±{args.flank:,} bp",
        fontsize=9,
    )
    fig.suptitle(
        f"{chrom}:{window[0]:,}-{window[1]:,} ({locus[0].query_strand})   "
        f"{len(refs)} ref CNEE(s) projected here; their categories over all loci: "
        f"initial {', '.join(initial)} -> final {', '.join(finals)}",
        x=0.01, ha="left", fontsize=11,
    )
    import matplotlib.patches as mpatches

    present = [s for s in STATUS_ORDER if any(r[1] == s for r in before)]
    fig.legend(
        handles=[mpatches.Patch(color=STATUS_COLORS[s], label=STATUS_LABELS[s]) for s in present],
        loc="upper right", ncol=min(3, len(present)), fontsize=8, frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 0.99, 0.94))
    tag = initial[0].replace("complex/many2many", "complex") if len(initial) == 1 else "mixed"
    out_path = out_dir / f"{locus_index:03d}_{tag}_{chrom}_{window[0]}-{window[1]}.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return {
        "figure": out_path.name, "query_locus": f"{chrom}:{window[0]}-{window[1]}", "strand": locus[0].query_strand,
        "ref_cnees": ",".join(refs), "initial_categories": ",".join(initial), "final_categories": ",".join(finals),
        "chains_in_view": len(before), "envelope_depth": scores.envelope_depth,
        "projections": len(locus), "kept": sum(v.accepted for v in verdicts),
        "verdict_mismatches_vs_run": ",".join(mismatches) or ".",
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("run_dir", type=Path, help="Output directory of a completed stola run")
    parser.add_argument("--ref-genome", required=True, type=Path, help="Reference genome the run used")
    parser.add_argument("--query-genome", required=True, type=Path, help="Query genome the run used")
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--random", type=int, metavar="N", help="Plot N random multi-projection loci")
    target.add_argument("--locus", type=parse_region, nargs="+", metavar="CHROM:START-END",
                        help="Plot every multi-projection locus overlapping these query regions")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", type=Path, help="[default: RUN_DIR/chain_pileups]")
    parser.add_argument("--max-rows", type=int, default=30, help="Chains drawn per panel [default: 30]")
    parser.add_argument("--flank", type=int, default=DEFAULT_FLANK_WINDOW,
                        help=f"The run's --flank-window [default: {DEFAULT_FLANK_WINDOW}]")
    parser.add_argument("--soft-depth-limit", type=int, default=DEFAULT_SOFT_DEPTH_LIMIT)
    parser.add_argument("--relative-score-threshold", type=float, default=DEFAULT_RELATIVE_SCORE_THRESHOLD)
    parser.add_argument("--flank-any-coverage-threshold", type=float, default=DEFAULT_FLANK_ANY_COVERAGE_THRESHOLD)
    args = parser.parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")

    run_dir = args.run_dir
    out_dir = args.out_dir or run_dir / "chain_pileups"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading run outputs ...", flush=True)
    ref_rows = read_tsv(run_dir / "ref_ortho_status.tsv")
    initial_by_ref = {r["ref_cnee_id"]: r["pre_pruning_category"] for r in ref_rows}
    final_by_ref = {r["ref_cnee_id"]: r["orthology_category"] for r in ref_rows}
    pruned_names = {r["projection_name"] for r in read_tsv(run_dir / "pruned_query_projections.tsv")}
    projections = _read_projection_evidence(run_dir / "overlap_stats.tsv").projections
    eligible = [p for p in projections if initial_by_ref.get(p.ref_id) in MULTI_CATEGORIES]
    loci = list(_iter_overlapping_projection_groups(eligible))
    print(f"  {len(loci):,} multi-projection loci", flush=True)

    if args.locus:
        selected = [
            locus for locus in loci
            if any(locus[0].query_chrom == chrom and min(p.query_start for p in locus) < end
                   and max(p.query_end for p in locus) > start for chrom, start, end in args.locus)
        ]
        if not selected:
            sys.exit("No multi-projection locus overlaps the given region(s)")
    else:
        rng = random.Random(args.seed)
        by_kind: dict[str, list] = defaultdict(list)
        for locus in loci:
            kinds = {initial_by_ref[p.ref_id] for p in locus}
            by_kind[kinds.pop() if len(kinds) == 1 else "mixed"].append(locus)
        pools = {kind: rng.sample(pool, len(pool)) for kind, pool in sorted(by_kind.items())}
        selected = []
        while len(selected) < args.random and any(pools.values()):
            for pool in pools.values():
                if pool and len(selected) < args.random:
                    selected.append(pool.pop())

    print("Loading chain index ...", flush=True)
    payload, records = load_index(run_dir / "work" / "chain_index.json")
    records_by_id = {r.chain_id: r for r in records}
    envelopes = build_query_envelopes(records)
    cache = ChainEnvelopeCache(payload["indexed_chain_path"], records_by_id)
    comparator = SequenceComparator(args.ref_genome, args.query_genome, run_dir / "work")
    rows = []
    try:
        for index, locus in enumerate(selected, start=1):
            result = plot_locus(
                locus, index, initial_by_ref, final_by_ref, pruned_names, records_by_id, envelopes, cache,
                comparator, args, out_dir,
            )
            rows.append(result)
            flag = "" if result["verdict_mismatches_vs_run"] == "." else "  ** verdicts differ from run **"
            print(f"  [{index}/{len(selected)}] {result['figure']}{flag}", flush=True)
    finally:
        cache.close()
        comparator.close()

    with (out_dir / "index.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} figures + index.tsv -> {out_dir}")


if __name__ == "__main__":
    main()
