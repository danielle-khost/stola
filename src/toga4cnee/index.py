"""Build and load transparent chain indexes."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path
from typing import Any

from .chain import parse_chain_header
from .models import ChainRecord

FORMAT_VERSION = 1


def _plain_chain_path_for_gzip(chain_path: Path, index_path: Path) -> Path:
    name = chain_path.name[:-3] if chain_path.name.endswith(".gz") else chain_path.name
    return index_path.parent / name


def _ensure_plain_chain(chain_file: str | Path, index_path: str | Path) -> Path:
    chain_path = Path(chain_file)
    if chain_path.suffix != ".gz":
        return chain_path.resolve()
    plain_path = _plain_chain_path_for_gzip(chain_path, Path(index_path))
    with gzip.open(chain_path, "rb") as src, plain_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return plain_path.resolve()


def scan_chain_file(chain_path: str | Path) -> list[ChainRecord]:
    records: list[ChainRecord] = []
    current_header: str | None = None
    current_start = 0
    byte_position = 0
    with Path(chain_path).open("rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b"chain "):
                if current_header is not None:
                    records.append(
                        parse_chain_header(
                            current_header,
                            start_byte=current_start,
                            byte_length=byte_position - current_start,
                        )
                    )
                current_header = raw_line.decode("utf-8").rstrip("\n")
                current_start = byte_position
            byte_position += len(raw_line)

    if current_header is not None:
        records.append(
            parse_chain_header(
                current_header,
                start_byte=current_start,
                byte_length=byte_position - current_start,
            )
        )
    if not records:
        raise ValueError(f"No chains found in {chain_path}")
    return records


def build_index(chain_file: str | Path, out_path: str | Path) -> dict[str, Any]:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    indexed_chain_path = _ensure_plain_chain(chain_file, out)
    records = scan_chain_file(indexed_chain_path)
    payload: dict[str, Any] = {
        "format_version": FORMAT_VERSION,
        "source_chain_path": str(Path(chain_file).resolve()),
        "indexed_chain_path": str(indexed_chain_path),
        "chain_count": len(records),
        "records": [record.to_dict() for record in records],
    }
    with out.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, separators=(",", ":"), sort_keys=True)
        handle.write("\n")
    return payload


def load_index(index_path: str | Path) -> tuple[dict[str, Any], list[ChainRecord]]:
    with Path(index_path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if int(payload.get("format_version", -1)) != FORMAT_VERSION:
        raise ValueError(
            f"Unsupported index format version: {payload.get('format_version')}"
        )
    records = [ChainRecord.from_dict(item) for item in payload["records"]]
    return payload, records
