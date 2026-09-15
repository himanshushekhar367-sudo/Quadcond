"""Adapter: i-motif transitional pH (pH_T) from the iM-Seeker supplementary data.

Yang, B. et al. (2024) *Prediction of DNA i-motifs via machine learning.*
Nucleic Acids Research 52(5):2188-2202. doi:10.1093/nar/gkae092
File: ``Supplementary Data Set.xlsx`` inside ``gkae092_supplemental_files.zip``

Two sheets:

``Supplementary Data Set 1``  171 constructs with a measured transitional pH.
``Supplementary Data Set 2``  the 120 of those whose putative i-motif the
                              authors' searcher could locate, given as the
                              detected sub-sequence.

**Only Sheet 1 is ingested as measurements.** Sheet 2 is the same 120
measurements re-expressed as extracted motifs; ingesting both would enter each
pH_T value twice under two different sequence spellings, inflating n and putting
near-identical rows on both sides of a cross-validation split. The detected
motif is instead attached to its parent record as a QC flag, so nothing is lost.

Buffer, from the paper's methods: 10 mM sodium cacodylate with 100 mM KCl,
titrated across pH 4-8.

**The record's pH is deliberately left unset, not set to pH_T.** An earlier
version stored ``condition.ph = ph_t`` on the reasoning that pH_T is the pH at
which the structure is half folded. That is true and it is also target leakage:
the pH_T head then reads its own answer out of a condition feature, and reports
R^2 = 0.99 for doing so. A pH titration does not happen *at* one pH -- it sweeps
pH 4 to 8 and reports the midpoint. So the condition here is the salt, the label
is the midpoint, and the titration range is recorded as a QC flag.

Sequences written with the authors' deletion notation (a Δ character marking a
deleted base) are rejected rather than stripped: removing the Δ would silently
invent a sequence that was never measured.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from ...conditions import Condition
from ..db import Record
from ._columns import read_table, resolve

SOURCE = "im_pht_biophysical"
DOI = "10.1093/nar/gkae092"
URL = "https://github.com/YANGB1/iM-Seeker"

# Methods buffer for the whole panel.
BUFFER = {"k": 100.0, "na": 10.0, "li_nh4": 0.0, "mg": 0.0, "temperature": 25.0}

SHEET_MEASUREMENTS = "Supplementary Data Set 1"
SHEET_MOTIFS = "Supplementary Data Set 2"

SPEC = {
    "sequence": ["sequences", "sequence", "seq", "dna sequence", "i-motif sequence"],
    "ph_t": ["pht", "ph_t", "ph t", "transitional ph", "transition ph"],
    "name": ["new_name", "original name", "name", "id", "label"],
}
MOTIF_SPEC = {
    "source_item": ["source data items", "source", "parent"],
    "motif": ["putative im detected by putative-im-searcher", "putative im", "detected"],
    "ph_t": ["pht", "ph_t"],
}

_DELETION = re.compile(r"[Δ∆]")


def _clean_sequence(raw) -> tuple[str | None, list[str]]:
    if raw is None:
        return None, ["sequence_missing"]
    s = str(raw).strip().upper()
    if not s or s == "NAN":
        return None, ["sequence_missing"]
    if _DELETION.search(s):
        return None, ["sequence_uses_deletion_notation"]
    s = re.sub(r"\s+", "", s).replace("U", "T")
    if not re.fullmatch(r"[ACGT]+", s):
        return None, ["sequence_has_non_nucleotide_characters"]
    return s, []


def _motif_index(path: str | Path) -> dict[str, str]:
    """Map construct name -> the detected putative i-motif, from Sheet 2."""
    try:
        import pandas as pd

        df = pd.read_excel(path, sheet_name=SHEET_MOTIFS, header=1)
    except Exception:
        return {}
    cols = resolve(df.columns, MOTIF_SPEC)
    if not cols.get("source_item") or not cols.get("motif"):
        return {}
    out: dict[str, str] = {}
    for _, r in df.iterrows():
        item = str(r[cols["source_item"]])
        motif = str(r[cols["motif"]]).strip().upper()
        # "ACA(C1/9=T)_2|TGTTCCC..." -> name before the pipe
        name = item.split("|")[0].strip()
        if name and motif and motif != "NAN":
            out[name] = motif
    return out


def load(path: str | Path, *, dry_run: bool = False, sheet: str | None = None):
    """Parse Sheet 1 into pH_T records, annotated with Sheet 2's detected motifs."""
    import pandas as pd

    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls", ".xlsm"}:
        try:
            df = pd.read_excel(path, sheet_name=sheet or SHEET_MEASUREMENTS, header=1)
        except ValueError:  # single-sheet file, or renamed sheet
            df = pd.read_excel(path, header=1)
    else:
        df = read_table(path)

    cols = resolve(df.columns, SPEC, required=["sequence", "ph_t"])
    if dry_run:
        return cols

    motifs = _motif_index(path) if path.suffix.lower().startswith(".xls") else {}

    out: list[Record] = []
    stats = {"rows": len(df), "rejected_sequence": 0, "rejected_pht": 0,
             "kept": 0, "with_detected_motif": 0}

    for i, r in df.iterrows():
        seq, seq_flags = _clean_sequence(r[cols["sequence"]])
        if seq is None:
            stats["rejected_sequence"] += 1
            continue
        try:
            pht = float(r[cols["ph_t"]])
        except (TypeError, ValueError):
            stats["rejected_pht"] += 1
            continue
        if math.isnan(pht):
            stats["rejected_pht"] += 1
            continue

        flags = list(seq_flags)
        if not (3.0 <= pht <= 9.0):
            flags.append("pht_outside_measured_range")

        name = str(r[cols["name"]]).strip() if cols.get("name") else f"row{i}"
        motif = motifs.get(name)
        if motif:
            flags.append(f"detected_putative_im={motif}")
            stats["with_detected_motif"] += 1

        # pH is NOT set from pH_T -- see the module docstring. The measurement is
        # a sweep, and storing its midpoint as the condition would let the pH_T
        # head predict its own target.
        flags.append("ph_titration_range=4.0-8.0")
        cond = Condition.from_mapping(BUFFER, track_imputed=True)
        out.append(
            Record(
                sequence=seq,
                kind="iM",
                source=SOURCE,
                evidence_tier="experimental",
                condition=cond,
                label_class="biophysical",
                folded=1,          # pH_T is measured on a sequence that does fold
                ph_t=pht,
                method="UV/CD pH titration; 10 mM Na cacodylate + 100 mM KCl, pH 4-8",
                source_doi=DOI,
                source_id=name,
                organism="in vitro (synthetic oligonucleotide)",
                qc_flags=flags,
            )
        )
        stats["kept"] += 1

    load.last_stats = stats
    return out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="i-motif transitional pH (pH_T): 171 constructs, iM-Seeker Supplementary "
              "Data Set 1",
        doi=DOI,
        url=URL,
        evidence_tier="experimental",
        notes="Buffer: 10 mM sodium cacodylate + 100 mM KCl, titrated across pH 4-8. "
              "The record's condition pH is deliberately LEFT UNSET and flagged imputed, "
              "with the titration range recorded as a QC flag. An earlier version of this "
              "note claimed the pH was set to the measured pH_T; the adapter never does "
              "that, and must not -- pH_T is the label, so writing it into a condition "
              "feature lets the pH_T head read its own answer and report R^2 = 0.99 for "
              "doing so. Sheet 2's extracted motifs are attached as flags, not ingested "
              "separately.",
    )
