"""Small shared data models for TOGA4CNEE."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class ChainRecord:
    chain_id: str
    score: int
    ref_chrom: str
    ref_size: int
    ref_strand: str
    ref_start: int
    ref_end: int
    query_chrom: str
    query_size: int
    query_strand: str
    query_start: int
    query_end: int
    start_byte: int
    byte_length: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ChainRecord":
        return cls(
            chain_id=str(data["chain_id"]),
            score=int(data["score"]),
            ref_chrom=str(data["ref_chrom"]),
            ref_size=int(data["ref_size"]),
            ref_strand=str(data["ref_strand"]),
            ref_start=int(data["ref_start"]),
            ref_end=int(data["ref_end"]),
            query_chrom=str(data["query_chrom"]),
            query_size=int(data["query_size"]),
            query_strand=str(data["query_strand"]),
            query_start=int(data["query_start"]),
            query_end=int(data["query_end"]),
            start_byte=int(data["start_byte"]),
            byte_length=int(data["byte_length"]),
        )


@dataclass(frozen=True)
class BedInterval:
    cnee_uid: str
    cnee_id: str
    chrom: str
    start: int
    end: int
    line_number: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class ChainBlock:
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int

    @property
    def ref_length(self) -> int:
        return self.ref_end - self.ref_start
