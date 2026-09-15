"""Adapter: 5DUVMA -- the condition-resolved i-motif stability landscape.

*High-throughput measurement and prediction of the i-motif DNA stability
landscape.* Nucleic Acids Research 2026, 54(4):gkag110.
doi:10.1093/nar/gkag110
File: ``gkag110_supplemental_file.pdf`` (the supporting information)

This is the dataset the i-motif side of QuadCond was missing, and it is the
mirror image of the iM-Seeker panel. iM-Seeker gives sequence breadth at a
single buffer: 160 different sequences, one condition. 5DUVMA gives condition
*depth* for two sequences: Tel21C and C9 measured across pH 4.35-6.60 and ten
ionic strengths from 0.01 M to 1 M, by high-throughput UV melting.

Neither alone can support a condition-aware i-motif model. With only iM-Seeker,
``im_pht`` learns sequence and has no basis to move with salt -- the ablation
measured that at +0.001 R^2. With only 5DUVMA it would learn the salt response
of two constructs and generalise to nothing. Together the two panels span both
axes, which is the whole point.

Tables ingested, parsed straight from the PDF text so the extraction is
reproducible rather than hand-copied:

======  ==========================================  ==============  ==================
table   content                                     label           varying condition
======  ==========================================  ==============  ==================
S1      Tel21C T_1/2 avg                            ``tm``          pH x ionic strength
S4      Tel21C pH_1/2 avg                           ``ph_t``        temperature x ionic strength
S9      C9 melting T_1/2 (main iM transition)       ``tm``          pH x ionic strength
S10     C9 annealing T_1/2                          ``tm``          pH x ionic strength
S12     C9 pH_1/2, melting                          ``ph_t``        temperature x ionic strength
S13     C9 pH_1/2, annealing                        ``ph_t``        temperature x ionic strength
======  ==========================================  ==============  ==================

Two conventions, both recorded on every row:

**Ionic strength is entered as K+.** The buffer is Britton-Robinson with KCl
added to set the stated ionic strength, and KCl is the dominant 1:1 contributor,
so ``k`` is set to the tabulated ionic strength in mM. The Britton-Robinson
component itself carries some Na+ that this ignores; the assumption is flagged
rather than buried.

**pH is not set on pH_1/2 rows.** Same rule as the iM-Seeker adapter: pH_1/2 is
the label, and copying it into the condition would let the head read its own
answer. Rows measuring pH_1/2 carry temperature and salt as conditions and leave
pH unset.

Cells marked ``a`` in the source ("the quality of the data is not sufficient to
allow an accurate determination") are skipped, not zero-filled.
"""
from __future__ import annotations

import re
from pathlib import Path

from ...conditions import Condition
from ..db import Record

SOURCE = "duvma_im_stability_landscape"
DOI = "10.1093/nar/gkag110"
URL = "https://academic.oup.com/nar/article/54/4/gkag110/8474393"

# Sequences as given in the paper: Tel21C = [C3TAA]3C3, C9 = [C9T3]3C9.
TEL21C = "CCC" + "TAACCC" * 3          # 21 nt
C9 = ("C" * 9 + "TTT") * 3 + "C" * 9   # 45 nt

# Each table: (page index, sequence, label column, row variable, process note)
TABLES = [
    {"table": "S1", "page": 4, "seq": TEL21C, "name": "Tel21C",
     "label": "tm", "row_var": "ph", "process": "average of melting and annealing"},
    {"table": "S4", "page": 6, "seq": TEL21C, "name": "Tel21C",
     "label": "ph_t", "row_var": "temperature", "process": "average of melting and annealing"},
    {"table": "S9", "page": 11, "seq": C9, "name": "C9",
     "label": "tm", "row_var": "ph", "process": "melting", "multirow": True},
    {"table": "S10", "page": 12, "seq": C9, "name": "C9",
     "label": "tm", "row_var": "ph", "process": "annealing"},
    {"table": "S12", "page": 14, "seq": C9, "name": "C9",
     "label": "ph_t", "row_var": "temperature", "process": "melting"},
    {"table": "S13", "page": 15, "seq": C9, "name": "C9",
     "label": "ph_t", "row_var": "temperature", "process": "annealing"},
]

_HEADER_IS = re.compile(r"(\d+(?:\.\d+)?)\s*M\b")
# A value is a number or the literal 'a'. The source contains at least one
# typo ("64..09"), so a doubled decimal point is tolerated and repaired.
_TOKEN = re.compile(r"(?<![\w.])(\d+\.{1,2}\d+|\d+|a)(?![\w])")


def _clean_number(tok: str) -> float | None:
    if tok == "a":
        return None
    tok = tok.replace("..", ".")
    try:
        return float(tok)
    except ValueError:
        return None


def _ionic_strengths(text: str) -> list[float]:
    """Column header: the ten ionic strengths, in molar."""
    head = text.split("\n")
    for i, line in enumerate(head):
        found = _HEADER_IS.findall(line)
        if len(found) >= 5:
            # the header can wrap over two lines
            nxt = _HEADER_IS.findall(head[i + 1]) if i + 1 < len(head) else []
            vals = [float(v) for v in found + nxt]
            return vals
    # fall back: scan the whole page for "<n> M" tokens before the first data row
    return [float(v) for v in _HEADER_IS.findall(text)][:10]


def _parse_grid(text: str, spec: dict) -> list[tuple[float, float, float]]:
    """Return (row_value, ionic_strength_M, measurement) triples."""
    ionic = _ionic_strengths(text)
    if not ionic:
        return []
    n_cols = len(ionic)
    out: list[tuple[float, float, float]] = []

    for raw_line in text.split("\n"):
        line = raw_line.strip()
        if not line or line.lower().startswith(("table", "all ", "a indicates")):
            continue
        if spec.get("multirow"):
            # C9 melting: "<pH> iM v1 ... v10", then "iM1 ...", "iM0 ...".
            # Only the main iM transition is ingested.
            m = re.match(r"^(\d+\.\d+)\s+iM\s+(.*)$", line)
            if not m:
                continue
            row_val, rest = float(m.group(1)), m.group(2)
        else:
            m = re.match(r"^(\d+(?:\.\d+)?)\s+(.*)$", line)
            if not m:
                continue
            row_val, rest = float(m.group(1)), m.group(2)

        toks = _TOKEN.findall(rest)
        if len(toks) < n_cols:
            continue
        toks = toks[:n_cols]
        for is_m, tok in zip(ionic, toks):
            val = _clean_number(tok)
            if val is not None:
                out.append((row_val, is_m, val))
    return out


def load(path: str | Path, *, dry_run: bool = False):
    """Parse every ingestable grid out of the supporting-information PDF."""
    import pypdf

    reader = pypdf.PdfReader(str(path))
    pages = [(p.extract_text() or "") for p in reader.pages]

    stats: dict[str, int] = {}
    out: list[Record] = []

    for spec in TABLES:
        idx = spec["page"]
        if idx >= len(pages):
            continue
        text = pages[idx]
        if f"Table S{spec['table'][1:]}" not in text.replace(" ", "").replace("TableS", "Table S"):
            # tolerate the source's inconsistent "Table S1 2." spacing
            if spec["table"].lower().replace("s", "s ") not in text.lower():
                pass  # fall through and try to parse anyway
        triples = _parse_grid(text, spec)
        stats[spec["table"]] = len(triples)
        if dry_run:
            continue

        for row_val, ionic_m, value in triples:
            k_mm = ionic_m * 1000.0
            flags = [
                f"source_table={spec['table']}",
                f"construct={spec['name']}",
                f"process={spec['process']}",
                f"ionic_strength_M={ionic_m:g}",
                "ionic_strength_entered_as_K+ (Britton-Robinson buffer with KCl; the "
                "buffer's own Na+ contribution is not separated out)",
            ]
            cond_map = {"k": k_mm, "na": 0.0, "li_nh4": 0.0, "mg": 0.0}
            tm = ph_t = None
            if spec["label"] == "tm":
                cond_map["ph"] = row_val
                cond_map["temperature"] = 25.0
                tm = value
                if not (0.0 <= value <= 120.0):
                    flags.append("tm_outside_physical_range")
            else:
                # pH_1/2 row: temperature is the condition, pH stays unset.
                cond_map["temperature"] = row_val
                ph_t = value
                flags.append("ph_titration; pH not set as a condition (it is the label)")
                if not (3.0 <= value <= 9.0):
                    flags.append("pht_outside_measured_range")

            out.append(
                Record(
                    sequence=spec["seq"],
                    kind="iM",
                    source=SOURCE,
                    evidence_tier="experimental",
                    condition=Condition.from_mapping(cond_map, track_imputed=True),
                    label_class="biophysical",
                    folded=1,
                    tm=tm,
                    ph_t=ph_t,
                    method=f"high-throughput UV melting (5DUVMA), {spec['process']}",
                    source_doi=DOI,
                    source_id=f"{spec['name']}_{spec['table']}_"
                              f"{spec['row_var']}{row_val:g}_IS{ionic_m:g}M",
                    organism="in vitro (synthetic oligonucleotide)",
                    qc_flags=flags,
                )
            )

    load.last_stats = stats
    return stats if dry_run else out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="5DUVMA: i-motif thermal and pH stability of Tel21C and C9 across pH and "
              "ten ionic strengths",
        doi=DOI,
        url=URL,
        evidence_tier="experimental",
        notes="Condition depth for two sequences, complementary to the iM-Seeker panel's "
              "sequence breadth at one buffer. Ionic strength entered as K+ (Britton-"
              "Robinson buffer with KCl). pH_1/2 rows deliberately leave the condition "
              "pH unset, because pH_1/2 is the label.",
    )
