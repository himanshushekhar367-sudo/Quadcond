#!/usr/bin/env python3
"""Worked examples on sequences the field knows well.

Not a benchmark -- a sanity surface. If a model's behaviour on hTelo, c-MYC and
the BCL2 promoter does not match what is in the literature, everything else it
says is suspect.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from quadcond.conditions import Condition
from quadcond.evaluate import condition_response, response_sensitivity
from quadcond.models.predict import Predictor

# Well-characterised sequences. Provenance is the point: each one is a construct
# that appears repeatedly in the primary literature, so behaviour on it is checkable.
CASES = [
    ("hTelo (22AG)", "AGGGTTAGGGTTAGGGTTAGGG",
     "Human telomeric repeat. Hybrid in K+, antiparallel (basket) in Na+ -- the "
     "canonical demonstration that topology is cation dependent."),
    ("c-MYC Pu22", "TGAGGGTGGGTAGGGTGGGTAA",
     "c-MYC NHE III1 promoter G4. Parallel, unusually stable."),
    ("c-KIT1", "AGGGAGGGCGCTGGGAGGAGGG",
     "c-KIT promoter G4, parallel with an unusual long-loop architecture."),
    ("BCL2 Pu39 core", "GGGCGCGGGAGGAAGGGGGCGGG",
     "BCL2 P1 promoter G4 region; mixed hybrid populations reported."),
    ("VEGF", "GGGGCGGGCCGGGGGCGGGGT",
     "VEGF promoter G4, parallel."),
    ("hTelo C-rich", "CCCTAACCCTAACCCTAACCCT",
     "Complement of the telomeric repeat -- the classic i-motif, folded below "
     "roughly pH 6."),
    ("BCL2 C-rich", "CCCGCCCCCTTCCTCCCGCGCCC",
     "C-rich BCL2 promoter strand; i-motif formation reported near neutral pH."),
]

CONDITIONS = [
    ("100 mM K+", Condition.preset("k100")),
    ("100 mM Na+", Condition.preset("na100")),
    ("100 mM Li+", Condition.preset("li100")),
    ("physiological", Condition.preset("physiological")),
    ("acidic tumour-like", Condition.preset("acidic_tumour")),
    ("endolysosomal", Condition.preset("endolysosomal")),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--out", default="artifacts/case_studies.json")
    a = ap.parse_args()

    pred = Predictor.load(a.model, a.db if Path(a.db).exists() else None)
    report: dict = {"cases": [], "sweeps": {}, "sensitivity": {}}

    table_rows = []
    for name, seq, note in CASES:
        case = {"name": name, "sequence": seq, "note": note, "by_condition": {}}
        for cname, cond in CONDITIONS:
            r = pred.predict(seq, cond, n_neighbours=2)[0]
            slim = {}
            for hn, hp in r["predictions"].items():
                if "probability" in hp:
                    slim[hn] = round(hp["probability"], 4)
                elif "posterior" in hp:
                    slim[hn] = {"call": hp["argmax"], **{k: round(v, 4)
                                                         for k, v in hp["posterior"].items()}}
                else:
                    slim[hn] = round(hp["value"], 3)
                    if "folded_fraction_at_condition" in hp:
                        slim[hn + ":folded_fraction"] = hp["folded_fraction_at_condition"]
            case["by_condition"][cname] = slim
            row = {"case": name, "condition": cname}
            for k, v in slim.items():
                row[k] = v["call"] if isinstance(v, dict) else v
            table_rows.append(row)
        case["nearest_measurements"] = pred.predict(seq, CONDITIONS[0][1], n_neighbours=3)[0]["evidence"]
        report["cases"].append(case)

    # K+ titration for the telomeric G4; pH titration for the telomeric i-motif
    k_rows = condition_response(pred, CASES[0][1], "k", np.linspace(0, 150, 16),
                               Condition.preset("k100"))
    ph_rows = condition_response(pred, CASES[5][1], "ph", np.linspace(4.0, 8.0, 17),
                                Condition.preset("im_reference"))
    report["sweeps"] = {"hTelo_K_titration": k_rows, "hTelo_C_pH_titration": ph_rows}
    report["sensitivity"] = {
        "hTelo_vs_K": response_sensitivity(k_rows, "k"),
        "hTeloC_vs_pH": response_sensitivity(ph_rows, "ph"),
    }

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(report, indent=2, default=str))

    df = pd.DataFrame(table_rows)
    csv = Path(a.out).with_suffix(".csv")
    df.to_csv(csv, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {a.out} and {csv}")

    print("\nSensitivity to [K+] across 0-150 mM (peak-to-peak):")
    for k, v in report["sensitivity"]["hTelo_vs_K"].items():
        print(f"  {k:<40} {v:.4f}")
    print("\nSensitivity to pH across 4.0-8.0:")
    for k, v in report["sensitivity"]["hTeloC_vs_pH"].items():
        print(f"  {k:<40} {v:.4f}")


if __name__ == "__main__":
    main()
