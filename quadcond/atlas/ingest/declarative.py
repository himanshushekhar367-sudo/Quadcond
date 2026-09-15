"""Declarative JSON adapters for arbitrary supplementary and lab tables.

The Python adapters in this package cover the sources we know about. This one
covers the ones we do not: you describe a spreadsheet in a small JSON file --
which column means what, what units it is in, how its categorical values map --
and QuadCond ingests it without anyone writing code.

The unit handling is the point. Half the frustration of harmonising G4/i-motif
literature is that one paper reports 0.1 M KCl, the next 100 mM, and a third
gives ionic strength only. A declarative adapter records the conversion
explicitly, so the atlas stores millimolar and the file keeps saying whatever it
said.

Example (``my_source.json``)::

    {
      "id_prefix": "SMITH_2025",
      "kind": "G4",
      "columns": {
        "sequence": "oligo_sequence",
        "tm": "melting_temperature",
        "k": "potassium_M",
        "ph": "buffer_pH",
        "folded": "fold_call"
      },
      "units": {"k": "M", "tm": "C"},
      "constants": {"na": 0, "mg": 0, "temperature": 25,
                    "method": "UV melting", "nucleic_acid": "DNA"},
      "value_maps": {"folded": {"yes": 1, "no": 0}},
      "source": "smith_2025_uv_melts",
      "source_doi": "10.1000/example",
      "evidence_tier": "experimental",
      "label_class": "biophysical"
    }

    quadcond atlas ingest-json my_source.json melts.xlsx --dry-run

Concentration units convert to mM; temperature to degrees Celsius. An unknown
unit is an error, never a silent pass-through -- a mis-scaled salt column is the
kind of mistake that quietly poisons a condition-aware model.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ...conditions import CONDITION_FIELDS, Condition
from ..db import Record
from ._columns import read_table

# to millimolar
CONCENTRATION_TO_MM = {
    "M": 1000.0, "mol/L": 1000.0, "molar": 1000.0,
    "mM": 1.0, "millimolar": 1.0,
    "uM": 1e-3, "µM": 1e-3, "micromolar": 1e-3,
    "nM": 1e-6, "nanomolar": 1e-6,
}
CONCENTRATION_FIELDS = {"k", "na", "li_nh4", "mg", "strand_conc"}
TEMPERATURE_UNITS = {"C", "degC", "celsius", "K", "kelvin", "F", "fahrenheit"}
LABEL_FIELDS = {"folded", "topology", "tm", "dg", "ph_t"}


class AdapterError(ValueError):
    """Raised when an adapter spec is malformed or a unit is unrecognised."""


def _convert(field: str, value: float, unit: str | None) -> float:
    if unit is None:
        return value
    if field in CONCENTRATION_FIELDS:
        try:
            return value * CONCENTRATION_TO_MM[unit]
        except KeyError:
            raise AdapterError(
                f"unknown concentration unit {unit!r} for field {field!r}; "
                f"known units: {sorted(CONCENTRATION_TO_MM)}"
            ) from None
    if field in {"temperature", "tm"}:
        if unit in {"C", "degC", "celsius"}:
            return value
        if unit in {"K", "kelvin"}:
            return value - 273.15
        if unit in {"F", "fahrenheit"}:
            return (value - 32.0) * 5.0 / 9.0
        raise AdapterError(
            f"unknown temperature unit {unit!r}; known units: {sorted(TEMPERATURE_UNITS)}"
        )
    if unit in {"percent", "%", "fraction"}:
        return value * 100.0 if unit == "fraction" else value
    raise AdapterError(f"no unit conversion defined for field {field!r} (unit {unit!r})")


def load_spec(path: str | Path) -> dict[str, Any]:
    spec = json.loads(Path(path).read_text())
    if "columns" not in spec or "sequence" not in spec.get("columns", {}):
        raise AdapterError("adapter spec needs a columns.sequence mapping")
    if not spec.get("source"):
        raise AdapterError("adapter spec needs a 'source' key (the dataset id)")
    known = set(CONDITION_FIELDS) | LABEL_FIELDS | {"sequence", "name", "kind",
                                                    "method", "nucleic_acid", "notes"}
    for field in list(spec["columns"]) + list(spec.get("constants", {})):
        if field not in known:
            raise AdapterError(
                f"adapter maps unknown field {field!r}; known fields: {sorted(known)}"
            )
    for field, unit in spec.get("units", {}).items():
        if field not in spec["columns"]:
            raise AdapterError(f"units declared for {field!r} but that column is not mapped")
        _convert(field, 1.0, unit)  # fail fast on an unknown unit
    return spec


def preview(spec: dict[str, Any], path: str | Path) -> dict[str, Any]:
    """What the adapter will do, without touching the database."""
    df = read_table(path)
    cols = spec["columns"]
    missing = [c for c in cols.values() if c not in df.columns]
    return {
        "table": str(path),
        "rows": int(len(df)),
        "table_columns": list(df.columns),
        "mapping": cols,
        "missing_columns": missing,
        "units": spec.get("units", {}),
        "constants": spec.get("constants", {}),
        "value_maps": spec.get("value_maps", {}),
        "source": spec["source"],
        "evidence_tier": spec.get("evidence_tier", "experimental"),
        "label_class": spec.get("label_class", "biophysical"),
        "status": "ready" if not missing else "columns missing from the table",
    }


def load(path: str | Path, *, spec: str | Path | dict, dry_run: bool = False):
    spec = spec if isinstance(spec, dict) else load_spec(spec)
    if dry_run:
        return preview(spec, path)

    df = read_table(path)
    cols = spec["columns"]
    units = spec.get("units", {})
    consts = spec.get("constants", {})
    vmaps = spec.get("value_maps", {})
    missing = [c for c in cols.values() if c not in df.columns]
    if missing:
        raise AdapterError(f"table {path} is missing mapped columns {missing}")

    kind_const = consts.get("kind", spec.get("kind", "G4"))
    out: list[Record] = []
    for i, row in df.iterrows():
        def value(field):
            if field in cols:
                raw = row[cols[field]]
                if field in vmaps:
                    raw = vmaps[field].get(str(raw).strip(), vmaps[field].get(raw, raw))
                return raw
            return consts.get(field)

        seq = value("sequence")
        if seq is None or str(seq).strip().lower() in {"", "nan"}:
            continue

        def number(field):
            v = value(field)
            if v is None:
                return None
            try:
                v = float(v)
            except (TypeError, ValueError):
                return None
            if v != v:
                return None
            return _convert(field, v, units.get(field))

        cond = Condition.from_mapping({f: number(f) for f in CONDITION_FIELDS},
                                      track_imputed=True)
        folded = value("folded")
        try:
            folded = int(folded) if folded is not None and str(folded) != "nan" else None
        except (TypeError, ValueError):
            folded = None
        topo = value("topology")
        topo = str(topo).strip().lower() if topo not in (None, "") else None
        if topo in {"nan", "none", ""}:
            topo = None

        flags = [f"declarative_adapter:{spec['source']}"]
        if spec.get("label_class"):
            flags.append(f"label_class={spec['label_class']}")
        if cond.imputed:
            flags.append("condition_partially_imputed:" + ",".join(cond.imputed))

        out.append(
            Record(
                sequence=str(seq).strip(),
                kind=str(value("kind") or kind_const),
                source=spec["source"],
                evidence_tier=spec.get("evidence_tier", "experimental"),
                condition=cond,
                nucleic_acid=str(value("nucleic_acid") or "DNA"),
                folded=folded,
                topology=topo,
                tm=number("tm"),
                dg=number("dg"),
                ph_t=number("ph_t"),
                method=str(value("method") or "unspecified"),
                source_doi=spec.get("source_doi"),
                source_id=str(value("name") or f"{spec.get('id_prefix', 'ROW')}_{i}"),
                qc_flags=flags,
            )
        )
    return out


def register(atlas, spec: dict[str, Any]) -> None:
    atlas.register_source(
        spec["source"],
        title=spec.get("title", f"Declarative import: {spec['source']}"),
        doi=spec.get("source_doi", ""),
        url=spec.get("source_url", ""),
        evidence_tier=spec.get("evidence_tier", "experimental"),
        notes=spec.get("notes", "Ingested through a declarative JSON adapter. "
                                f"label_class={spec.get('label_class', 'biophysical')}."),
    )
