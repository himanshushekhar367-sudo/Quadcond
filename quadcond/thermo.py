"""Two-state transition shapes.

These are not fitted models -- they are the standard two-state descriptions of
a melting or a pH titration, used to turn a *point* prediction (Tm, pH_T) into
the quantity a biologist actually wants: the fraction folded at their
temperature and their pH.  Keeping them separate from the ML heads makes the
provenance obvious: the head supplies the midpoint, thermodynamics supplies the
shape, and the shape parameter is a stated assumption you can override.
"""
from __future__ import annotations

import math

R = 1.987204e-3  # kcal / (mol K)

# Typical van 't Hoff enthalpy for an intramolecular G4 unfolding transition.
# Reported values span roughly 30-60 kcal/mol depending on sequence and method;
# 45 is a mid-range default and is exposed as a parameter, not hidden.
DEFAULT_DH_G4 = 45.0

# Hill coefficient for i-motif pH titrations. The C.C+ transition is
# cooperative; n is commonly fitted between 2 and 4.
DEFAULT_HILL_IM = 2.5


def folded_fraction_thermal(temperature_c: float, tm_c: float,
                            dh_kcal: float = DEFAULT_DH_G4) -> float:
    """Two-state van 't Hoff folded fraction at a temperature, given Tm.

    theta = 1 / (1 + exp[-dH/R (1/T - 1/Tm)])
    """
    T = temperature_c + 273.15
    Tm = tm_c + 273.15
    if T <= 0 or Tm <= 0:
        return float("nan")
    x = (dh_kcal / R) * (1.0 / T - 1.0 / Tm)
    x = max(min(x, 60.0), -60.0)
    return 1.0 / (1.0 + math.exp(-x))


def folded_fraction_ph(ph: float, ph_t: float, hill: float = DEFAULT_HILL_IM) -> float:
    """Hill-type folded fraction for an i-motif at a given pH, given pH_T.

    theta = 1 / (1 + 10^[n (pH - pH_T)])
    """
    x = hill * (ph - ph_t)
    x = max(min(x, 30.0), -30.0)
    return 1.0 / (1.0 + 10.0 ** x)


def tm_shift_estimate(tm_ref_c: float, cation_ref_mM: float, cation_mM: float,
                      slope_per_decade: float = 12.0) -> float:
    """Rough Tm at a different monovalent cation concentration.

    G4 Tm rises approximately linearly with log10[K+] over the mM-to-100 mM
    range; the slope is sequence dependent and roughly 10-15 C per decade for
    intramolecular G4s.  This is a stated heuristic used ONLY when no
    condition-aware model is available, and every output that uses it is
    labelled ``heuristic``.
    """
    if cation_ref_mM <= 0 or cation_mM <= 0:
        return float("nan")
    return tm_ref_c + slope_per_decade * math.log10(cation_mM / cation_ref_mM)
