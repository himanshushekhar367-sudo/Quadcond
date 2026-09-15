"""Explicit experimental-condition handling.

The central premise of QuadCond is that a non-canonical structure prediction is
meaningless without the buffer it refers to.  Every record in the atlas and
every prediction request therefore carries a :class:`Condition`.

Units are fixed and enforced here so that heterogeneous literature sources can
be harmonised on ingestion:

============ ============================================
field        unit
============ ============================================
k            mM  potassium
na           mM  sodium
li_nh4       mM  lithium / ammonium (the classic "does not
                 support G4" control cations)
mg           mM  magnesium (divalent)
ph           dimensionless
temperature  degrees Celsius
crowder_pct  % w/v of PEG200/300 or equivalent crowder
strand_conc  uM  oligonucleotide strand concentration
============ ============================================
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict, field
from typing import Any, Mapping, Sequence

# Reference conditions used when a source does not report a variable.
# These are *defaults with provenance*, not silent zeros: every ingestion
# adapter records which fields were imputed in `condition_imputed`.
DEFAULTS: dict[str, float] = {
    "k": 100.0,
    "na": 0.0,
    "li_nh4": 0.0,
    "mg": 0.0,
    "ph": 7.0,
    "temperature": 25.0,
    "crowder_pct": 0.0,
    "strand_conc": 5.0,
}

CONDITION_FIELDS: tuple[str, ...] = tuple(DEFAULTS)

# Physiological / commonly-cited reference points, used by the CLI and the app.
PRESETS: dict[str, dict[str, float]] = {
    "physiological": {"k": 140.0, "na": 12.0, "mg": 1.0, "ph": 7.4, "temperature": 37.0},
    "nuclear": {"k": 140.0, "na": 10.0, "mg": 0.5, "ph": 7.2, "temperature": 37.0},
    # Acidic micro-environments: endosome/lysosome-like and the mildly acidic
    # extracellular pH reported for solid tumours.  Relevant to i-motifs.
    "acidic_tumour": {"k": 140.0, "na": 12.0, "mg": 1.0, "ph": 6.5, "temperature": 37.0},
    "endolysosomal": {"k": 50.0, "na": 20.0, "mg": 0.5, "ph": 5.5, "temperature": 37.0},
    # Buffers dominating the biophysical literature.
    "k100": {"k": 100.0, "ph": 7.0, "temperature": 25.0},
    "k50": {"k": 50.0, "ph": 7.0, "temperature": 25.0},
    "na100": {"k": 0.0, "na": 100.0, "ph": 7.0, "temperature": 25.0},
    "li100": {"k": 0.0, "li_nh4": 100.0, "ph": 7.0, "temperature": 25.0},
    # The G4STAB web-database "cancer-like" ionic condition.
    "cancer_like": {"k": 70.0, "na": 30.0, "ph": 7.0, "temperature": 25.0},
    # i-motif reference buffer used by iM-Seeker's biophysical set.
    "im_reference": {"k": 100.0, "na": 10.0, "ph": 5.8, "temperature": 25.0},
}


@dataclass(frozen=True)
class Condition:
    """A fully specified experimental condition."""

    k: float = DEFAULTS["k"]
    na: float = DEFAULTS["na"]
    li_nh4: float = DEFAULTS["li_nh4"]
    mg: float = DEFAULTS["mg"]
    ph: float = DEFAULTS["ph"]
    temperature: float = DEFAULTS["temperature"]
    crowder_pct: float = DEFAULTS["crowder_pct"]
    strand_conc: float = DEFAULTS["strand_conc"]
    imputed: tuple[str, ...] = field(default=(), compare=False)

    # ------------------------------------------------------------------ ctors
    @classmethod
    def from_mapping(cls, m: Mapping[str, Any] | None, *, track_imputed: bool = True) -> "Condition":
        m = dict(m or {})
        vals: dict[str, float] = {}
        imputed: list[str] = []
        for f in CONDITION_FIELDS:
            v = m.get(f)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                vals[f] = DEFAULTS[f]
                imputed.append(f)
            else:
                vals[f] = float(v)
        return cls(**vals, imputed=tuple(imputed) if track_imputed else ())

    @classmethod
    def preset(cls, name: str, **overrides: float) -> "Condition":
        if name not in PRESETS:
            raise KeyError(f"unknown preset {name!r}; known: {sorted(PRESETS)}")
        base = dict(PRESETS[name])
        base.update(overrides)
        return cls.from_mapping(base, track_imputed=False)

    # ------------------------------------------------------------- properties
    @property
    def monovalent(self) -> float:
        """Total monovalent cation concentration (mM)."""
        return self.k + self.na + self.li_nh4

    @property
    def ionic_strength(self) -> float:
        """Approximate ionic strength (mM), 1:1 salts plus a 2:1 Mg salt.

        I = 1/2 sum(c_i z_i^2).  For MgCl2 this gives 3*[Mg].
        """
        return 0.5 * (self.k + self.na + self.li_nh4) * 2 + 3.0 * self.mg

    @property
    def k_fraction(self) -> float:
        """Fraction of monovalent cation that is K+ (the G4-competent cation)."""
        tot = self.monovalent
        return self.k / tot if tot > 0 else 0.0

    @property
    def g4_competent_cation(self) -> float:
        """K+ plus a down-weighted Na+ contribution.

        K+ sits in the centre of the G-tetrad stack with near-ideal
        coordination geometry; Na+ stabilises far less well and Li+/NH4+ are
        conventionally treated as non-supporting controls.  The 0.35 weight is
        a *prior*, not a fitted parameter -- models see the raw ion
        concentrations as well and can learn their own weighting.
        """
        return self.k + 0.35 * self.na

    # ------------------------------------------------------------------ utils
    def to_dict(self) -> dict[str, float]:
        d = asdict(self)
        d.pop("imputed", None)
        return d

    def replace(self, **kw: float) -> "Condition":
        d = self.to_dict()
        d.update(kw)
        return Condition.from_mapping(d, track_imputed=False)

    def label(self) -> str:
        parts = [f"K{self.k:g}", f"Na{self.na:g}"]
        if self.li_nh4:
            parts.append(f"Li/NH4 {self.li_nh4:g}")
        if self.mg:
            parts.append(f"Mg{self.mg:g}")
        parts += [f"pH{self.ph:g}", f"{self.temperature:g}C"]
        if self.crowder_pct:
            parts.append(f"crowd{self.crowder_pct:g}%")
        return " ".join(parts)


def condition_distance(a: Condition, b: Condition) -> float:
    """Scaled distance between two conditions, for atlas nearest-neighbour search.

    Each term is normalised by a scale that reflects how much of that variable
    it takes to matter experimentally: ~50 mM for monovalent cations, ~0.5 pH
    units, ~10 C.  A distance of 1.0 is therefore "about one meaningful step".
    """
    d = 0.0
    d += ((a.k - b.k) / 50.0) ** 2
    d += ((a.na - b.na) / 50.0) ** 2
    d += ((a.li_nh4 - b.li_nh4) / 50.0) ** 2
    d += ((a.mg - b.mg) / 2.0) ** 2
    d += ((a.ph - b.ph) / 0.5) ** 2
    d += ((a.temperature - b.temperature) / 10.0) ** 2
    d += ((a.crowder_pct - b.crowder_pct) / 10.0) ** 2
    return math.sqrt(d)


def condition_matrix(conditions: Sequence[Condition]) -> "list[list[float]]":
    """Feature matrix for a batch of conditions (see features.condition_features)."""
    from .features import condition_features

    return [condition_features(c) for c in conditions]
