"""Half-open interval helpers."""

from __future__ import annotations


def overlap_length(start1: int, end1: int, start2: int, end2: int) -> int:
    """Return the positive overlap length between two half-open intervals."""
    # Clamp the half-open intersection to zero for non-overlapping intervals.
    return max(0, min(end1, end2) - max(start1, start2))


def fraction(numerator: int, denominator: int) -> float:
    """Return a guarded floating-point fraction."""
    # Avoid division errors for zero-length intervals or empty spans.
    return numerator / denominator if denominator else 0.0
