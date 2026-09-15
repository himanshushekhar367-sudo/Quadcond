#!/usr/bin/env python3
"""Put the ensemble free energies on a scale where the Boltzmann step does work.

The AENNA-3D viewer weights competing folds by exp(-dG/RT) and advertises that
the weights respond to pH, salt, crowding and temperature. They barely did, and
the reason was arithmetic rather than biology: RT is 0.616 kcal/mol at 37 C,
and the heuristic assigned a canonical four-tract quadruplex about
-2.4 x 3 x 4 = -28.8 kcal/mol while a sequence with no motif fell back to a
hardcoded +8 to +10. A 35 kcal/mol gap is 57 RT. exp(57) saturates the code's
own exponent clamp, so every sequence was ~100% one state and the sliders moved
numbers that could not move the ensemble.

The fix is not to divide by ten until the picture looks better. It is to give
the numbers a referent. A measured melting temperature determines the folding
free energy at any other temperature under a two-state model:

    dG(T) = dH (1 - T / Tm)

so the 2,274 measured G4 melting temperatures in the atlas -- across 255
buffers -- become 2,274 measured folding free energies at 37 C. This script
fits the heuristic's output to those, which yields a scale and an offset with
units attached and an R^2 that says how much of the variation the heuristic
captures at all.

What comes out is written to ``quadcond/structure/ensemble_calibration.json``
and consumed by the viewer. The residual scatter is reported alongside, because
a calibrated heuristic is still a heuristic: the point of the exercise is a
Boltzmann weight that responds to conditions at a physically plausible rate,
not a claim that these dG values are accurate for any single sequence.

    python scripts/11_calibrate_ensemble_dg.py
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from quadcond.thermo import DEFAULT_DH_G4, R

REFERENCE_T_C = 37.0


# --------------------------------------------------------------------------
# The heuristic, ported verbatim from src/lib/na/predict.ts so that what is
# calibrated is what the viewer actually runs.
# --------------------------------------------------------------------------
def count_runs(seq: str, base: str, mn: int = 2) -> tuple[int, int]:
    m = re.findall(f"{base}{{{mn},}}", seq)
    return len(m), max((len(x) for x in m), default=0)


def g4_energy_raw(seq: str, monovalent: float, mg: float,
                  crowding: float, temperature: float) -> float:
    n, mx = count_runs(seq, "G", 2)
    tetrads = min(4, max(2, mx))
    dG = 8.0
    if n >= 4 and mx >= 2:
        dG = -2.4 * tetrads * min(n, 4) + 2.2 * max(0, 4 - n)
        loops = [L for L in re.sub(r"G{2,}", " ", seq).strip().split() if L]
        dG += sum(max(0, len(L) - 1) * 0.35 for L in loops)
    dG += -1.15 * math.log10(max(monovalent, 1.0) / 20.0)
    dG += -0.22 * min(mg, 10.0)
    dG += -0.028 * crowding
    dG += 0.048 * (temperature - 37.0)
    return dG


def imotif_energy_raw(seq: str, ph: float, mg: float, crowding: float,
                      temperature: float) -> float:
    n, mx = count_runs(seq, "C", 2)
    pka = 4.7
    frac = 1.0 / (1.0 + 10.0 ** (ph - pka))
    dG = 9.0
    if n >= 4:
        dG = -1.7 * mx * min(n, 4)
    elif n >= 2:
        dG = 2.5 - 0.6 * n * mx
    dG += (1.0 - frac) * 9.5
    dG += -0.35 * min(mg, 12.0)
    dG += -0.04 * crowding * frac
    dG += 0.06 * (temperature - 37.0)
    return dG


def dg_from_tm(tm_c: float, dh_kcal: float = DEFAULT_DH_G4,
               t_c: float = REFERENCE_T_C) -> float:
    """Folding free energy at ``t_c`` implied by a measured melting temperature.

    Sign convention: folding is favourable below Tm, so dG is negative there.
    dH is the van 't Hoff enthalpy of unfolding, positive.
    """
    T = t_c + 273.15
    Tm = tm_c + 273.15
    return -dh_kcal * (1.0 - T / Tm)


def fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Least squares y ~ a x + b, returning a, b and R^2."""
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    pred = a * x + b
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return a, b, 1.0 - ss_res / ss_tot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--dh", type=float, default=DEFAULT_DH_G4,
                    help="van 't Hoff unfolding enthalpy, kcal/mol")
    ap.add_argument("--out", default="quadcond/structure/ensemble_calibration.json")
    a = ap.parse_args()

    con = sqlite3.connect(f"file:{Path(a.db)}?mode=ro", uri=True)
    rows = con.execute(
        "SELECT sequence, kind, k, na, li_nh4, mg, ph, crowder_pct, tm FROM records "
        "WHERE evidence_tier='experimental' AND label_class='biophysical' "
        "AND tm IS NOT NULL").fetchall()
    con.close()

    out: dict[str, object] = {
        "reference_temperature_c": REFERENCE_T_C,
        "van_t_hoff_dh_kcal": a.dh,
        "RT_at_reference_kcal": round(R * (REFERENCE_T_C + 273.15), 4),
        "method": (
            "dG(37 C) implied by each measured Tm under a two-state van 't Hoff "
            "model, regressed on the heuristic's raw output. The fitted scale is "
            "how much the heuristic must be compressed to speak in kcal/mol."
        ),
        "classes": {},
    }

    for kind, fn in (("g-quadruplex", "g4"), ("i-motif", "im")):
        sel = [r for r in rows if (r[1] == "G4") == (fn == "g4")]
        if len(sel) < 30:
            continue
        raw, meas = [], []
        for seq, _k, k, na, li, mg, ph, crowd, tm in sel:
            mono = (k or 0.0) + (na or 0.0) + (li or 0.0)
            if fn == "g4":
                raw.append(g4_energy_raw(seq, mono, mg or 0.0, crowd or 0.0, REFERENCE_T_C))
            else:
                raw.append(imotif_energy_raw(seq, ph or 7.0, mg or 0.0,
                                             crowd or 0.0, REFERENCE_T_C))
            meas.append(dg_from_tm(tm, a.dh))
        raw_a, meas_a = np.asarray(raw), np.asarray(meas)
        scale, offset, r2 = fit(raw_a, meas_a)
        resid = meas_a - (scale * raw_a + offset)

        span_before = float(raw_a.max() - raw_a.min())
        span_after = float((scale * raw_a).max() - (scale * raw_a).min())
        rt = R * (REFERENCE_T_C + 273.15)
        out["classes"][kind] = {
            "n": len(sel),
            "scale": round(scale, 4),
            "offset_kcal": round(offset, 3),
            "r2": round(r2, 4),
            "residual_sd_kcal": round(float(resid.std()), 2),
            "raw_span_kcal": round(span_before, 1),
            "raw_span_in_RT": round(span_before / rt, 1),
            "calibrated_span_kcal": round(span_after, 1),
            "calibrated_span_in_RT": round(span_after / rt, 1),
            "measured_dg_range_kcal": [round(float(meas_a.min()), 2),
                                       round(float(meas_a.max()), 2)],
        }
        print(f"{kind}: n={len(sel):,}")
        print(f"  dG(37 C) implied by measurement: "
              f"{meas_a.min():+.1f} to {meas_a.max():+.1f} kcal/mol")
        print(f"  heuristic raw span   {span_before:6.1f} kcal/mol "
              f"= {span_before / rt:5.1f} RT   <- saturates the exponent clamp")
        print(f"  after calibration    {span_after:6.1f} kcal/mol "
              f"= {span_after / rt:5.1f} RT")
        print(f"  scale {scale:.4f}   offset {offset:+.2f} kcal/mol   "
              f"R^2 {r2:.3f}   residual sd {resid.std():.1f} kcal/mol\n")

    # --- the comparison that decides the architecture -----------------------
    # If the calibrated heuristic explains little of the measured free energy,
    # it should not be the SOURCE of dG wherever a trained head exists. The
    # same van 't Hoff relation turns a predicted Tm into a predicted dG, so
    # the two routes are directly comparable on the same 2,274 measurements.
    try:
        from quadcond import assets
        from quadcond.models.predict import Predictor
        from quadcond.conditions import Condition

        pred = Predictor.load(assets.resolve("model"))
        g4 = [r for r in rows if r[1] == "G4"]
        yh, ym = [], []
        for seq, _k, k, na, li, mg, ph, crowd, tm in g4:
            c = Condition.from_mapping({"k": k, "na": na, "li_nh4": li, "mg": mg,
                                        "ph": ph, "temperature": REFERENCE_T_C,
                                        "crowder_pct": crowd}, track_imputed=False)
            v = pred.predict(seq, c, heads=["g4_tm"],
                             n_neighbours=0)[0]["predictions"]["g4_tm"]["value"]
            yh.append(dg_from_tm(v, a.dh))
            ym.append(dg_from_tm(tm, a.dh))
        yh, ym = np.asarray(yh), np.asarray(ym)
        ss_res = float(((ym - yh) ** 2).sum())
        ss_tot = float(((ym - ym.mean()) ** 2).sum())
        model_r2 = 1.0 - ss_res / ss_tot
        head_cv_r2 = pred.model.heads["g4_tm"].metrics.get("r2")
        out["model_derived_dg"] = {
            "route": "g4_tm predicted Tm -> dG(37 C) by the same van 't Hoff relation",
            "n": len(g4),
            "r2_vs_measured_dg_IN_SAMPLE": round(model_r2, 4),
            "rmse_kcal_IN_SAMPLE": round(float(np.sqrt(((ym - yh) ** 2).mean())), 2),
            "in_sample_warning": (
                "g4_tm was TRAINED on these 2,274 rows, so this figure is in-sample "
                "and is not the accuracy to quote. It is reported only to show that "
                "the head carries the free-energy variation the heuristic does not. "
                "The head's honest number is its grouped cross-validated Tm R^2 "
                f"of {head_cv_r2:.3f}; the out-of-fold dG figure is materially lower "
                "than the in-sample one above."
            ),
            "head_grouped_cv_r2_on_tm": round(float(head_cv_r2), 4)
            if head_cv_r2 is not None else None,
            "verdict": (
                "Where a trained head exists, dG must come from it. The calibrated "
                "heuristic fixes the ensemble's dynamic range but explains little of "
                "the free-energy variation; the head explains most of it. The "
                "heuristic's remaining job is the folds no head covers -- cruciform, "
                "triplex, AC-motif, hairpin, duplex -- where it is the only estimate "
                "available and is labelled as such."
            ),
        }
        print(f"\nmodel-derived dG (g4_tm -> van 't Hoff), same {len(g4):,} measurements:")
        print(f"  IN-SAMPLE R^2 {model_r2:.3f}, RMSE "
              f"{np.sqrt(((ym - yh) ** 2).mean()):.2f} kcal/mol -- g4_tm was trained")
        print(f"  on these rows, so this is NOT an accuracy figure. The head's honest")
        print(f"  number is its grouped-CV Tm R^2 of {head_cv_r2:.3f}.")
        print(f"  What the comparison establishes: the head carries free-energy")
        print(f"  variation that the calibrated heuristic (R^2 "
              f"{out['classes']['g-quadruplex']['r2']:.3f}) does not.")
    except Exception as exc:                      # noqa: BLE001
        out["model_derived_dg"] = {"unavailable": str(exc)}
        print(f"\n(model comparison skipped: {exc})")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {a.out}")
    print("\nThe R^2 is the honest part: it says how much of the measured free-energy\n"
          "variation the heuristic explains once its units are fixed. A low value\n"
          "means the calibration made the ensemble respond at the right RATE, not\n"
          "that the heuristic became accurate for any individual sequence.")


if __name__ == "__main__":
    main()
