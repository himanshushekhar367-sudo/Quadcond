"""Adapter: your own lab measurements.

This is the differentiating input.  Public tables are collapsed onto a handful
of buffers; a lab that measures the *same* oligo across a cation/pH grid
supplies exactly the axis every published predictor is missing.

Expected table (csv/tsv/xlsx), one row per measurement per condition::

    sequence,name,kind,method,k,na,li_nh4,mg,ph,temperature,strand_conc,
    crowder_pct,tm,ph_t,dg,folded,topology,replicate,notes

Only ``sequence`` and at least one label column are required.  ``kind`` is
inferred from base composition when absent.  Everything else is optional and
missing conditions are flagged as imputed rather than silently defaulted.

    quadcond atlas ingest lab my_cd_melts.xlsx --source deep_lab_cd_2026
"""
from __future__ import annotations

from pathlib import Path

from ...conditions import Condition
from ...motifs import clean
from ..db import Record
from ._columns import read_table, resolve

SPEC = {
    "sequence": ["sequence", "seq", "oligo", "oligonucleotide"],
    "name": ["name", "id", "label", "construct", "gene"],
    "kind": ["kind", "structure", "class", "type"],
    "method": ["method", "technique", "assay", "instrument"],
    "k": ["k", "k+", "kcl", "potassium"],
    "na": ["na", "na+", "nacl", "sodium"],
    "li_nh4": ["li", "linh4", "li/nh4", "nh4", "lithium", "ammonium"],
    "mg": ["mg", "mg2+", "mgcl2", "magnesium"],
    "ph": ["ph"],
    "temperature": ["temperature", "temp"],
    "strand_conc": ["strandconc", "strand concentration", "oligo conc", "concentration", "um"],
    "crowder_pct": ["crowder", "peg", "crowding"],
    "tm": ["tm", "melting temperature", "t_m"],
    "ph_t": ["pht", "ph_t", "transitional ph"],
    "dg": ["dg", "deltag", "free energy", "g37"],
    "folded": ["folded", "forms", "positive", "status"],
    "topology": ["topology", "conformation", "fold"],
    "replicate": ["replicate", "rep", "run"],
    "notes": ["notes", "comment", "comments"],
}

_TRUE = {"1", "true", "yes", "y", "folded", "positive"}
_FALSE = {"0", "false", "no", "n", "unfolded", "negative"}


def _infer_kind(seq: str) -> str:
    s = clean(seq)
    return "iM" if s.count("C") > s.count("G") else "G4"


def load(path: str | Path, *, source: str, doi: str = "", dry_run: bool = False,
         label_class: str = "biophysical"):
    df = read_table(path)
    cols = resolve(df.columns, SPEC, required=["sequence"])
    if dry_run:
        return cols
    labels = [c for c in ("tm", "ph_t", "dg", "folded", "topology") if cols.get(c)]
    if not labels:
        raise ValueError(
            "no label column found; need at least one of tm / ph_t / dg / folded / topology"
        )
    out: list[Record] = []
    for i, r in df.iterrows():
        seq = str(r[cols["sequence"]]).strip()
        if not seq or seq.lower() == "nan":
            continue
        kind = str(r[cols["kind"]]).strip() if cols.get("kind") else ""
        if kind.lower() in {"g4", "gq", "quadruplex", "g-quadruplex"}:
            kind = "G4"
        elif kind.lower() in {"im", "i-motif", "imotif"}:
            kind = "iM"
        else:
            kind = _infer_kind(seq)

        cond = Condition.from_mapping(
            {
                key: (r[cols[key]] if cols.get(key) else None)
                for key in ("k", "na", "li_nh4", "mg", "ph", "temperature",
                            "crowder_pct", "strand_conc")
            },
            track_imputed=True,
        )

        def num(field):
            if not cols.get(field):
                return None
            try:
                v = float(r[cols[field]])
            except (TypeError, ValueError):
                return None
            return None if v != v else v

        folded = None
        if cols.get("folded"):
            raw = str(r[cols["folded"]]).strip().lower()
            if raw in _TRUE:
                folded = 1
            elif raw in _FALSE:
                folded = 0

        topo = str(r[cols["topology"]]).strip().lower() if cols.get("topology") else None
        if topo in {"", "nan", "none"}:
            topo = None

        flags = ["lab_measurement"]
        if cond.imputed:
            flags.append("condition_partially_imputed:" + ",".join(cond.imputed))
        if cols.get("replicate"):
            flags.append(f"replicate={r[cols['replicate']]}")

        out.append(
            Record(
                sequence=seq,
                kind=kind,
                source=source,
                evidence_tier="experimental",
                label_class=label_class,
                condition=cond,
                folded=folded,
                topology=topo,
                tm=num("tm"),
                dg=num("dg"),
                ph_t=num("ph_t"),
                method=str(r[cols["method"]]) if cols.get("method") else "lab (unspecified)",
                source_doi=doi or None,
                source_id=str(r[cols["name"]]) if cols.get("name") else f"row{i}",
                qc_flags=flags,
            )
        )
    return out


def register(atlas, source: str, *, title: str = "", doi: str = "", notes: str = "") -> None:
    atlas.register_source(
        source,
        title=title or f"Lab dataset: {source}",
        doi=doi,
        evidence_tier="experimental",
        notes=notes or "In-house measurements ingested via the lab adapter.",
    )
