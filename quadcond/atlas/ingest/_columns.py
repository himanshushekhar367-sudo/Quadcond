"""Tolerant column matching for third-party supplementary tables.

Supplementary files are inconsistent: ``Tm``, ``T_m``, ``Tm (°C)``, ``melting
temperature`` all mean the same thing.  Rather than hard-coding one spelling
and failing on a file we have never seen, adapters declare *aliases* and this
resolver picks the best match, then reports exactly what it matched so the
ingest log is auditable.
"""
from __future__ import annotations

import re
from typing import Mapping, Sequence


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def resolve(
    columns: Sequence[str],
    spec: Mapping[str, Sequence[str]],
    *,
    required: Sequence[str] = (),
) -> dict[str, str | None]:
    """Map logical field -> actual column name.

    ``spec`` maps a logical field to candidate aliases, tried in order:
    exact-normalised match first, then substring match.
    """
    norm = {_norm(c): c for c in columns}
    out: dict[str, str | None] = {}
    for field, aliases in spec.items():
        hit = None
        for a in aliases:
            na = _norm(a)
            if na in norm:
                hit = norm[na]
                break
        if hit is None:
            for a in aliases:
                na = _norm(a)
                for nc, c in norm.items():
                    if na and na in nc:
                        hit = c
                        break
                if hit:
                    break
        out[field] = hit
    missing = [f for f in required if out.get(f) is None]
    if missing:
        raise ValueError(
            f"could not find columns for {missing}; available columns were {list(columns)}"
        )
    return out


def read_table(path):
    """Read csv/tsv/xlsx/xls transparently."""
    import pandas as pd
    from pathlib import Path

    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls", ".xlsm"}:
        return pd.read_excel(p)
    if p.suffix.lower() in {".tsv", ".tab"}:
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p)
