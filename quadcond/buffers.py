"""Parsing buffer descriptions written in prose into cation concentrations.

The G4 thermodynamic literature does not report ionic conditions in columns. It
reports them in a sentence: ``10 mM lithium cacodylate + 100 mM KCl``,
``20 mM Kpi + 80 mM KCl``, ``PBS (pH 7.4)``, ``100 mM K + unknown Li cacodylate
(pH 7.2)``. G4STAB's supplementary table carries 419 distinct such strings
across its 2,382 measurements, and the whole point of a condition-aware model is
that those strings are the input variables.

This module turns them into millimolar concentrations of K⁺, Na⁺, Li⁺/NH₄⁺ and
Mg²⁺, and — just as importantly — says when it could not.

Three rules keep it honest:

1. **Every component is attributed to a named species.** A term the species
   table does not recognise is reported in ``unparsed``, never silently dropped
   and never quietly folded into a cation. A row whose buffer contains an
   unparsed term with a concentration is flagged, so it can be excluded from
   training with one filter.
2. **Stoichiometry is explicit and documented.** Na₂HPO₄ contributes two Na⁺ per
   formula unit, KH₂PO₄ one K⁺. Where a real buffer is a mixture whose exact
   speciation depends on pH — "20 mM potassium phosphate" is some blend of
   KH₂PO₄ and K₂HPO₄ — the convention used is stated in ``PHOSPHATE_K_PER_MM``
   and recorded as an assumption on the parsed result, rather than being hidden
   in a magic number.
3. **"Unknown" means unknown.** ``100 mM K + unknown Li cacodylate`` yields
   K⁺ = 100 with a flag saying the lithium term was unquantified. It does not
   yield Li⁺ = 0, because those are different claims.

The parser is deliberately conservative: it would rather flag a row than guess
at it. Coverage on the G4STAB table is reported by
``scripts/06_ingest_g4stab.py``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

# A "20 mM potassium phosphate buffer, pH 7" is a K2HPO4/KH2PO4 mixture. At
# pH 7.0 the two are present in roughly a 2:3 ratio, giving about 1.4 K+ per
# phosphate; the usual convention in the G4 literature is to quote the total
# phosphate concentration. We use 1.5, state it, and attach the assumption to
# every row that depends on it so the choice is visible and can be changed.
PHOSPHATE_K_PER_MM = 1.5
PHOSPHATE_NA_PER_MM = 1.5
PHOSPHATE_LI_PER_MM = 1.5

# Phosphate-buffered saline, standard 1x formulation.
PBS = {"na": 157.0, "k": 4.45, "li_nh4": 0.0, "mg": 0.0}

# Species table. Order matters: the first pattern that matches a term wins, so
# more specific patterns (disodium phosphate) precede more general ones (sodium).
# Each entry maps a regex to the millimolar cation contribution per millimolar
# of the named species.
SPECIES: list[tuple[str, dict[str, float]]] = [
    # --- explicitly stoichiometric salts ------------------------------------
    (r"na2hpo4|disodium\s*(hydrogen\s*)?phosphate", {"na": 2.0}),
    (r"k2hpo4|dipotassium\s*(hydrogen\s*)?phosphate", {"k": 2.0}),
    (r"nah2po4|monosodium\s*phosphate", {"na": 1.0}),
    (r"kh2po4|monopotassium\s*phosphate", {"k": 1.0}),
    (r"li3po4|lithium\s*orthophosphate", {"li_nh4": 3.0}),
    (r"kno3|potassium\s*nitrate", {"k": 1.0}),
    (r"nano3|sodium\s*nitrate", {"na": 1.0}),
    (r"mgcl2|magnesium\s*chloride|mg\s*\(?2\+?\)?|mg2\+|\bmg\b", {"mg": 1.0}),
    # --- simple chlorides / hydroxides ---------------------------------------
    (r"\bkcl\b|potassium\s*chloride", {"k": 1.0}),
    (r"\bnacl\b|sodium\s*chloride", {"na": 1.0}),
    (r"\blicl\b|lithium\s*chloride", {"li_nh4": 1.0}),
    (r"nh4cl|ammonium\s*chloride|\bnh4\+?\b|ammonium\s*acetate", {"li_nh4": 1.0}),
    (r"\bkoh\b|potassium\s*hydroxide", {"k": 1.0}),
    (r"\bnaoh\b|sodium\s*hydroxide", {"na": 1.0}),
    (r"\blioh\b|lithium\s*hydroxide", {"li_nh4": 1.0}),
    # --- phosphate buffers quoted as a total ----------------------------------
    (r"(k\s*pi|kpi|kp\b|k\s*phosphate|potassium\s*phosphate|phosphate\s*.{0,12}potassium)",
     {"k": PHOSPHATE_K_PER_MM}),
    (r"(na\s*pi|napi|nap\b|na\s*phosphate|sodium\s*phosphate)", {"na": PHOSPHATE_NA_PER_MM}),
    (r"(li\s*pi|lipi|lipo4|li\s*phosphate|lithium\s*phosphate)", {"li_nh4": PHOSPHATE_LI_PER_MM}),
    # --- cacodylates and the arsenate analogue --------------------------------
    (r"(li|lithium)\s*(cacodylate|caco\b|aso4|arsenate)|licaco|lias|cacodyl.{0,20}lioh",
     {"li_nh4": 1.0}),
    (r"(na|sodium)\s*(cacodylate|caco\b)|nacaco|cacodyl.{0,20}naoh", {"na": 1.0}),
    (r"(k|potassium)\s*(cacodylate|caco\b)|kcaco|cacodyl.{0,20}koh", {"k": 1.0}),
    # --- buffers contributing no alkali cation --------------------------------
    (r"tris\s*[-\s]?hcl|\btris\b|hepes|tmaa|edta|boric|acetic|h3po4|dmso|spermidine|"
     r"citrate|mes\b|mops|pipes|bis-?tris|glycine|britton|pb\(no3\)2|nitrate\s*of\s*lead",
     {}),
    # --- bare cation names, last so a salt name is preferred -------------------
    (r"\bk\+?\b|potassium", {"k": 1.0}),
    (r"\bna\+?\b|sodium", {"na": 1.0}),
    (r"\bli\+?\b|lithium", {"li_nh4": 1.0}),
]
_COMPILED = [(re.compile(p, re.I), c) for p, c in SPECIES]

# Species we recognise but deliberately refuse to quantify, each with the reason.
# These are the honest middle ground between "parsed" and "no idea": we know what
# the term is, and we know we cannot turn it into a number we would stand behind.
AMBIGUOUS: list[tuple[str, str]] = [
    (r"phosphate\s*(buffer)?$|^\s*\d[\d.]*\s*m?m?\s*phosphate",
     "phosphate buffer with no named counter-ion (usually K+ in this literature, "
     "but not stated)"),
    (r"cacodylate|caco\b",
     "cacodylate buffer with no named counter-ion"),
    (r"\brb\b|rbcl|rubidium",
     "Rb+ supports G4 formation but is not one of the four cations this model "
     "represents"),
    (r"\bcs\b|cscl|caesium|cesium",
     "Cs+ is not one of the four cations this model represents"),
    (r"\bcpk\b",
     "unidentified buffer abbreviation 'CPK'"),
]
_AMBIGUOUS = [(re.compile(p, re.I), why) for p, why in AMBIGUOUS]

# "10 mM", "0.2 mM", "100mM", "1 M", "150 uM", or a bare number ("100 Na").
_QTY = re.compile(
    r"(?P<val>\d+(?:\.\d+)?)\s*(?P<unit>mM|millimolar|M\b|molar|uM|µM|micromolar|nM)?",
    re.I,
)
_UNIT_TO_MM = {None: 1.0, "": 1.0, "mm": 1.0, "millimolar": 1.0,
               "m": 1000.0, "molar": 1000.0,
               "um": 1e-3, "µm": 1e-3, "micromolar": 1e-3, "nm": 1e-6}

_PH = re.compile(r"p\s*h\s*[:=]?\s*(?P<ph>\d+(?:\.\d+)?)", re.I)
# Note the required slash in "n/a": an earlier version wrote it as `\bn/?a\b`,
# which cheerfully matched the element symbol "Na" and silently discarded every
# bare sodium term in the table. Regressions of that shape are why parse_buffer
# reports its components instead of just returning four numbers.
_UNKNOWN = re.compile(
    r"unknown|unsure|not\s*(reported|stated|given)|(?<![a-z])n/a(?![a-z])|not\s*applicable",
    re.I,
)
# A term that only states a pH ("at pH 7.0") carries no species.
_PH_ONLY = re.compile(r"^\W*(at\s*)?p\s*h\s*[:=]?\s*\d+(\.\d+)?\W*$", re.I)


@dataclass
class ParsedBuffer:
    """Cation concentrations in mM, plus an audit trail."""

    k: float = 0.0
    na: float = 0.0
    li_nh4: float = 0.0
    mg: float = 0.0
    ph: float | None = None
    components: list[tuple[str, float, str]] = field(default_factory=list)
    unparsed: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    source_text: str = ""

    @property
    def monovalent(self) -> float:
        return self.k + self.na + self.li_nh4

    @property
    def usable(self) -> bool:
        """Did we recover a cation composition we would train on?

        Requires at least one identified cation-bearing component and no
        quantified term we failed to attribute. A buffer that is entirely
        non-ionic (pure Tris) is usable only if something else supplied ions.
        """
        return self.monovalent + self.mg > 0 and not self.unparsed

    def to_dict(self) -> dict:
        return {
            "k": round(self.k, 4), "na": round(self.na, 4),
            "li_nh4": round(self.li_nh4, 4), "mg": round(self.mg, 4),
            "ph": self.ph, "usable": self.usable,
            "components": self.components, "unparsed": self.unparsed,
            "flags": self.flags, "assumptions": self.assumptions,
            "source_text": self.source_text,
        }


def _split_terms(text: str) -> list[str]:
    """Split a buffer description into additive terms.

    Splits on '+' and on ' and ', but not on commas inside parentheses, and
    keeps a trailing parenthetical (which usually carries the pH) attached to
    the term it followed.
    """
    text = text.replace("＋", "+")
    parts = re.split(r"\s*\+\s*|\s+and\s+", text)
    return [p.strip() for p in parts if p.strip()]


def parse_buffer(text: str, *, reported_ph: float | None = None) -> ParsedBuffer:
    """Parse one buffer description into cation concentrations (mM)."""
    out = ParsedBuffer(source_text=str(text))
    if text is None or str(text).strip().lower() in {"", "nan", "none"}:
        out.flags.append("no_buffer_description")
        return out
    raw = str(text).strip()

    # pH: prefer the value reported in the table's own column; fall back to a
    # "(pH 7.2)" written into the description. A reported pH of 0 is the
    # table's way of saying unknown, not an acidic buffer.
    if reported_ph is not None and reported_ph > 0:
        out.ph = float(reported_ph)
    else:
        m = _PH.search(raw)
        if m:
            out.ph = float(m.group("ph"))
            out.flags.append("ph_recovered_from_buffer_text")
        else:
            out.flags.append("ph_unknown")

    if re.search(r"\bpbs\b|phosphate[- ]buffered saline", raw, re.I):
        out.k, out.na = PBS["k"], PBS["na"]
        out.components.append(("PBS 1x (standard formulation)", 1.0, "k+na"))
        out.assumptions.append(
            "PBS expanded to the standard 1x formulation: 137 mM NaCl, 2.7 mM KCl, "
            "10 mM Na2HPO4, 1.8 mM KH2PO4"
        )
        return out

    for term in _split_terms(raw):
        # Strip a trailing parenthetical that only carries pH, so it does not
        # look like an unparsed component.
        core = re.sub(r"\((?:[^()]*p\s*h[^()]*)\)", " ", term, flags=re.I).strip()
        if not core or _PH_ONLY.match(core):
            continue

        species_hit = None
        for pat, contrib in _COMPILED:
            if pat.search(core):
                species_hit = (pat.pattern, contrib)
                break

        qty = _QTY.search(core)
        has_number = qty is not None and qty.group("val") is not None

        if species_hit is None:
            # Recognised but not quantifiable is a different answer from
            # unrecognised, and only the latter should disqualify a row.
            for pat, why in _AMBIGUOUS:
                if pat.search(core):
                    out.flags.append(f"ambiguous_component:{core}")
                    out.assumptions.append(f"{core}: {why}; contributes 0 to the "
                                           f"cation vector")
                    break
            else:
                if has_number and not _UNKNOWN.search(core):
                    out.unparsed.append(core)
            continue

        pattern, contrib = species_hit
        if not contrib:  # recognised, contributes no alkali/alkaline cation
            out.components.append((core, 0.0, "non-ionic buffer"))
            continue

        if _UNKNOWN.search(core) or not has_number:
            out.flags.append(f"unquantified_component:{core}")
            continue

        val = float(qty.group("val"))
        unit = (qty.group("unit") or "").lower().replace("μ", "u")
        mm = val * _UNIT_TO_MM.get(unit, 1.0)
        if qty.group("unit") is None:
            out.flags.append(f"unit_assumed_mM:{core}")

        for ion, per in contrib.items():
            setattr(out, ion, getattr(out, ion) + mm * per)
        out.components.append((core, mm, ",".join(f"{i}x{p:g}" for i, p in contrib.items())))
        if any(k in pattern for k in ("pi", "phosphate")) and "kh2po4" not in pattern:
            note = (f"phosphate buffers quoted as a total are expanded at "
                    f"{PHOSPHATE_K_PER_MM} cations per mM")
            if note not in out.assumptions:
                out.assumptions.append(note)

    if out.monovalent + out.mg == 0:
        out.flags.append("no_cation_recovered")
    return out


def parse_series(texts, reported_ph=None) -> list[ParsedBuffer]:
    reported_ph = list(reported_ph) if reported_ph is not None else [None] * len(list(texts))
    return [parse_buffer(t, reported_ph=p) for t, p in zip(texts, reported_ph)]


def coverage_report(parsed: list[ParsedBuffer]) -> dict:
    """How much of a table the parser actually recovered."""
    n = len(parsed)
    usable = sum(1 for p in parsed if p.usable)
    flagged = sum(1 for p in parsed if p.unparsed)
    no_cat = sum(1 for p in parsed if "no_cation_recovered" in p.flags)
    unq = sum(1 for p in parsed if any(f.startswith("unquantified") for f in p.flags))
    no_ph = sum(1 for p in parsed if p.ph is None)
    terms: dict[str, int] = {}
    for p in parsed:
        for u in p.unparsed:
            terms[u] = terms.get(u, 0) + 1
    return {
        "rows": n,
        "usable": usable,
        "usable_fraction": round(usable / n, 4) if n else 0.0,
        "with_unparsed_terms": flagged,
        "no_cation_recovered": no_cat,
        "with_unquantified_component": unq,
        "ph_unknown": no_ph,
        "top_unparsed_terms": sorted(terms.items(), key=lambda kv: -kv[1])[:25],
    }
