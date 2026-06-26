"""Write lifted query-side CNEE projections as BED12."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def _bed_score(chain_score: object) -> int:
    # Clamp a chain score into the BED score range.
    try:
        score = int(chain_score)
    except (TypeError, ValueError):
        return 0
    return min(1000, max(0, score))


def row_to_bed12(row: dict[str, object]) -> str | None:
    """Convert an overlap-stats row into a query-side BED12 projection."""
    # Convert one aligned overlap-stat row into a legacy per-projection BED12 row.
    if int(row["aligned_bases"]) <= 0:
        return None
    query_start = row["query_aligned_start"]
    query_end = row["query_aligned_end"]
    if query_start == "." or query_end == ".":
        return None
    start = int(query_start)
    end = int(query_end)
    name = str(row["query_projection_name"])
    return "\t".join(
        [
            str(row["query_chrom"]),
            str(start),
            str(end),
            name,
            str(_bed_score(row["chain_score"])),
            str(row["query_strand"]),
            str(start),
            str(end),
            "0",
            str(row["query_block_count"]),
            str(row["query_block_sizes"]),
            str(row["query_block_starts"]),
        ]
    )


def write_query_bed(rows: Iterable[dict[str, object]], out_path: str | Path) -> int:
    # Write legacy BED12 rows for aligned overlap-stat projections.
    count = 0
    with Path(out_path).open("w", encoding="utf-8") as handle:
        for row in rows:
            line = row_to_bed12(row)
            if line is None:
                continue
            handle.write(line + "\n")
            count += 1
    return count
