"""Adapter: G4ShapePredictor topology dataset.

Liew, D., Lim, Z. W. & Yong, E. H. (2024) *Machine learning-based prediction of
DNA G-quadruplex folding topology with G4ShapePredictor*. Scientific Reports
14:24238.  doi:10.1038/s41598-024-74826-2
File: ``g4sp supplementary/G4 Dataset.xlsx`` in github.com/donn-liew/G4ShapePredictor

1,005 G4-forming sequences with experimentally assigned folding topology
(0 = parallel, 1 = antiparallel, 2 = hybrid).

Condition caveat, recorded on every row: the dataset was assembled from K+
buffer studies but does **not** report per-sequence salt, pH or temperature.
We therefore write the K+-buffer reference condition and flag every condition
field as imputed, plus a ``condition_not_reported`` QC flag.  This is exactly
the limitation QuadCond exists to expose -- a topology model trained here
cannot legitimately be queried at 10 mM K+ or pH 5.5 without an explicit
out-of-domain warning.
"""
from __future__ import annotations

import re
from pathlib import Path

from ...conditions import Condition
from ..db import Record

SOURCE = "g4sp_topology"
DOI = "10.1038/s41598-024-74826-2"
URL = "https://github.com/donn-liew/G4ShapePredictor"
TOPOLOGY_MAP = {0: "parallel", 1: "antiparallel", 2: "hybrid"}

_PDB = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _demojibake(text: str) -> str:
    """Repair UTF-8 text that was decoded once as cp1252/latin-1.

    The citation column of the source spreadsheet arrives with en dashes and
    accented author names mangled by exactly one bad decode, so author
    provenance reads as noise. Reversing that decode byte by byte -- cp1252 for
    the 0x80-0x9F specials, latin-1 for the rest -- and re-decoding as UTF-8
    recovers the original. Anything that does not round-trip cleanly is
    returned untouched, so ASCII citations and genuine non-Latin text are safe.
    """
    if not text or all(ord(c) < 0x80 for c in text):
        return text
    raw = bytearray()
    for ch in text:
        o = ord(ch)
        if o < 0x100:
            raw.append(o)
        else:
            try:
                raw += ch.encode("cp1252")
            except UnicodeEncodeError:
                return text
    try:
        fixed = bytes(raw).decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        return text
    return fixed.replace("\u00a0", " ")


def _method(comment: str) -> str:
    c = (comment or "").strip()
    if _PDB.match(c):
        return "PDB structure (NMR/X-ray)"
    if "in-house" in c.lower() and "cd" in c.lower():
        return "CD (in-house, G4SP authors)"
    if c:
        return "literature (CD/NMR)"
    return "unspecified"


def _source_id(comment: str) -> str:
    c = (comment or "").strip()
    if _PDB.match(c):
        return c.upper()
    m = re.match(r"^(PQS\d+)", c)
    if m:
        return m.group(1)
    return c[:120]


def load(path: str | Path) -> list[Record]:
    import pandas as pd

    df = pd.read_excel(path, index_col=0)
    required = {"Sequence", "Topology"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    # K+ reference buffer, every field flagged as imputed.
    cond = Condition.from_mapping({}, track_imputed=True)

    out: list[Record] = []
    for _, r in df.iterrows():
        topo = TOPOLOGY_MAP.get(int(r["Topology"]))
        if topo is None:
            continue
        comment = _demojibake(str(r.get("Comments", "") or ""))
        out.append(
            Record(
                sequence=str(r["Sequence"]),
                kind="G4",
                source=SOURCE,
                evidence_tier="experimental",
                label_class="biophysical",
                condition=cond,
                folded=1,  # every entry is an experimentally observed G4
                topology=topo,
                method=_method(comment),
                source_doi=DOI,
                source_id=_source_id(comment),
                organism="mixed (predominantly human)",
                qc_flags=["condition_not_reported", "k_buffer_assumed"],
            )
        )
    return out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="G4ShapePredictor curated G4 topology dataset (1,005 sequences)",
        doi=DOI,
        url=URL,
        evidence_tier="experimental",
        licence="see source repository",
        notes="Topology from CD/NMR/X-ray. Per-sequence buffer conditions NOT reported; "
              "K+ reference condition imputed on every row.",
    )
