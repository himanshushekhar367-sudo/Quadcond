"""Adapter: G4STAB Supplementary Table 1 -- experimental G4 melting temperatures.

Liew, D., Dharmatilleke, A. D., See, E. & Yong, E. H. (2025) *G4STAB: a
multi-input deep learning model to predict G-quadruplex thermodynamic stability
based on sequence and salt concentration.* Bioinformatics 41(10):btaf545.
doi:10.1093/bioinformatics/btaf545

2,382 melting measurements curated from 197 primary papers. This is the dataset
the condition-aware Tm head is built on, and the reason QuadCond exists: 446 of
its 1,068 distinct sequences were measured under two or more different buffers,
one of them under 67. That within-sequence condition contrast is what makes a
salt response learnable rather than confounded with sequence identity.

The table is not tidy, and three of its columns need real handling:

**Buffers are prose, not columns.** ``Comments`` holds 419 distinct strings like
``20 mM Kpi + 80 mM KCl`` or ``100 mM K + unknown Li cacodylate (pH 7.2)``.
:mod:`quadcond.buffers` parses these into explicit cation concentrations and
reports what it could not resolve; rows whose buffer contains an unrecognised
quantified term are flagged and excluded from the usable set.

**Not every Tm is a number.** 93 rows carry a censored bound (``>90``, ``<30``,
``low``) or two transition temperatures (``51/71``) for a multiphasic melt.
Coercing those to floats would inject fabricated values into the regression
target, so both classes are recorded with their raw text and a QC flag, and
their ``tm`` is left NULL. They are still ingested -- they are real evidence
that the sequence folds, just not a usable single midpoint -- so they train the
folding head while sitting out the Tm head.

**pH 0 means unknown.** The table encodes missing pH as 0.0. That is imputed
from the buffer text where the text states one, and flagged where it does not.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from ...buffers import coverage_report, parse_buffer
from ...conditions import Condition
from ..db import Record
from ._columns import read_table, resolve

SOURCE = "g4stab_experimental_tm"
DOI = "10.1093/bioinformatics/btaf545"
URL = "https://academic.oup.com/bioinformatics/article/41/10/btaf545/8266696"

SPEC = {
    "sequence": ["sequence", "seq", "dna sequence", "g4 sequence", "oligonucleotide"],
    "tm": ["temperature (deg c)", "tm", "t_m", "melting temperature", "temperature"],
    "error": ["error", "sd", "std", "uncertainty"],
    "ph": ["ph"],
    "buffer": ["comments", "buffer", "conditions", "notes"],
    "doi": ["doi", "reference", "citation", "source"],
    "strand_conc": ["dna conc (um)", "dna conc", "strand concentration", "concentration"],
    "wavelength": ["wavelength (nm)", "wavelength"],
}

_CENSORED = re.compile(r"^\s*(?P<op>[<>~]|greater than|less than)\s*(?P<val>\d+(\.\d+)?)\s*$", re.I)
_MULTIPHASIC = re.compile(r"^\s*(?P<a>\d+(\.\d+)?)\s*/\s*(?P<b>\d+(\.\d+)?)\s*$")
_QUALITATIVE = {"low", "high", "n/a", "na", "unknown", "nd", "-", ""}


def _clean_sequence(raw) -> tuple[str | None, list[str]]:
    """Normalise a sequence cell, or reject it.

    Several entries are line-wrapped in the source and arrive with a space in
    the middle of the sequence; those are joined. Anything still containing a
    non-nucleotide character is rejected rather than silently masked to N,
    because an N in the middle of a G-tract changes what the row means.
    """
    flags: list[str] = []
    if raw is None:
        return None, ["sequence_missing"]
    s = str(raw).strip().upper()
    if not s or s == "NAN":
        return None, ["sequence_missing"]
    if re.search(r"\s", s):
        s = re.sub(r"\s+", "", s)
        flags.append("sequence_whitespace_removed")
    s = s.replace("U", "T")
    if not re.fullmatch(r"[ACGT]+", s):
        return None, flags + ["sequence_has_non_nucleotide_characters"]
    if len(s) < 12:
        flags.append("sequence_shorter_than_12nt")
    return s, flags


def _parse_tm(raw) -> tuple[float | None, list[str]]:
    """Return (tm, flags). A censored bound or a biphasic pair yields no number."""
    if raw is None:
        return None, ["tm_missing"]
    if isinstance(raw, (int, float)) and not (isinstance(raw, float) and math.isnan(raw)):
        v = float(raw)
        return (v, []) if 0.0 <= v <= 120.0 else (None, ["tm_outside_physical_range"])
    text = str(raw).strip()
    if text.lower() in _QUALITATIVE:
        return None, [f"tm_qualitative:{text}"] if text else ["tm_missing"]
    m = _CENSORED.match(text)
    if m:
        return None, [f"tm_censored:{text}"]
    m = _MULTIPHASIC.match(text)
    if m:
        return None, [f"tm_multiphasic:{text}"]
    try:
        v = float(text)
    except ValueError:
        return None, [f"tm_unparsed:{text[:32]}"]
    return (v, []) if 0.0 <= v <= 120.0 else (None, ["tm_outside_physical_range"])


def _number(raw):
    try:
        v = float(raw)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def load(path: str | Path, *, dry_run: bool = False, require_usable_buffer: bool = True):
    """Parse the table into atlas records.

    ``require_usable_buffer`` drops rows whose buffer could not be resolved to a
    cation composition. Keep it True for training: a Tm with an unknown cation
    is not a condition-resolved measurement, and silently defaulting it to
    100 mM K+ would teach the model that the buffer does not matter.
    """
    df = read_table(path)
    cols = resolve(df.columns, SPEC, required=["sequence", "tm"])
    if dry_run:
        return cols

    buffer_col = cols.get("buffer")
    out: list[Record] = []
    stats = {"rows": len(df), "rejected_sequence": 0, "rejected_buffer": 0,
             "tm_censored": 0, "tm_multiphasic": 0, "tm_usable": 0, "kept": 0}
    parsed_buffers = []

    for i, r in df.iterrows():
        seq, seq_flags = _clean_sequence(r[cols["sequence"]])
        if seq is None:
            stats["rejected_sequence"] += 1
            continue

        reported_ph = _number(r[cols["ph"]]) if cols.get("ph") else None
        buf_text = r[buffer_col] if buffer_col else None
        pb = parse_buffer(buf_text, reported_ph=reported_ph)
        parsed_buffers.append(pb)
        if require_usable_buffer and not pb.usable:
            stats["rejected_buffer"] += 1
            continue

        tm, tm_flags = _parse_tm(r[cols["tm"]])
        if any(f.startswith("tm_censored") for f in tm_flags):
            stats["tm_censored"] += 1
        if any(f.startswith("tm_multiphasic") for f in tm_flags):
            stats["tm_multiphasic"] += 1
        if tm is not None:
            stats["tm_usable"] += 1

        cond = Condition.from_mapping(
            {
                "k": pb.k, "na": pb.na, "li_nh4": pb.li_nh4, "mg": pb.mg,
                "ph": pb.ph,
                "temperature": 25.0,   # the measurement is a melt; 25 C is the reference
                "strand_conc": _number(r[cols["strand_conc"]]) if cols.get("strand_conc") else None,
            },
            track_imputed=True,
        )

        flags = seq_flags + tm_flags + list(pb.flags)
        if pb.unparsed:
            flags.append("buffer_unparsed_terms:" + "|".join(pb.unparsed)[:120])
        if pb.assumptions:
            flags.append("buffer_assumptions:" + "; ".join(pb.assumptions)[:200])
        err = _number(r[cols["error"]]) if cols.get("error") else None
        if err is not None:
            flags.append(f"reported_error={err}")
        if cols.get("wavelength"):
            flags.append(f"wavelength={r[cols['wavelength']]}")

        out.append(
            Record(
                sequence=seq,
                kind="G4",
                source=SOURCE,
                evidence_tier="experimental",
                condition=cond,
                label_class="biophysical",
                folded=1,      # every row is a measured melting transition
                tm=tm,
                method="thermal melting (UV/CD), literature curation",
                source_doi=str(r[cols["doi"]]).strip() if cols.get("doi") else DOI,
                source_id=f"g4stab_row{i}",
                organism="in vitro (synthetic oligonucleotide)",
                qc_flags=flags,
            )
        )
        stats["kept"] += 1

    load.last_stats = stats
    load.last_buffer_coverage = coverage_report(parsed_buffers)
    return out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="G4STAB Supplementary Table 1: 2,382 experimental G4 melting temperatures "
              "with buffer composition parsed from the reported conditions",
        doi=DOI,
        url=URL,
        evidence_tier="experimental",
        notes="Cation concentrations parsed from the free-text buffer field by "
              "quadcond.buffers; censored (>90, <30) and multiphasic (51/71) melting "
              "values are ingested as folding evidence with tm=NULL rather than coerced "
              "to numbers.",
    )
