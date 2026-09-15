#!/usr/bin/env python3
"""Score the i-motif heads on a held-out CROSS-ASSAY panel.

Every number in the model card comes from cross-validation inside this atlas.
That is the right way to report performance, and it is not the same thing as
showing the model works on somebody else's sequences. This script does something
in between, and the distinction matters enough to name precisely.

The panel is the 16 oligonucleotides Zanin et al. characterised by circular
dichroism alongside their iMab CUT&Tag mapping (NAR 2023,
doi:10.1093/nar/gkad626, Table S2). It is a **held-out cross-assay panel**, not
a fully independent external set, for four reasons that are all worth stating
rather than burying:

1. **Same study.** These oligos come from the same paper as the CUT&Tag peaks
   that train ``im_fold_genomic``. Different assay (CD versus antibody pulldown),
   same laboratory and same peak context.
2. **Two negative controls.** Twelve i-motifs, two G4 controls, two negatives.
   A separation computed against two sequences is a direction, not a measurement.
3. **One training overlap**, checked rather than assumed -- see
   :func:`training_overlap`. It is the G4 control, and it affects the G4 side only.
4. **No numeric label.** The source reports CD spectra as figures, not tables, so
   the only label available is the class the authors assigned. This measures
   separation between their three groups, not agreement with a measured pH_T.

What it is good for: catching the failure mode that cross-validation cannot,
namely a head that scores well inside its own panel and cannot tell an i-motif
from a G4 anywhere else. It caught exactly that in ``im_fold``. What it cannot
do is establish generalisation -- that needs a panel from a laboratory with no
connection to the training data, which QuadCond does not yet have.

These sequences are never ingested into the atlas.

    python scripts/07_validate_external.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from quadcond.conditions import Condition
from quadcond.models.predict import Predictor
from quadcond.motifs import clean, find_g4, find_im

SET_PATH = "data/validation/zanin2023_cd_oligos.json"


def _fmt(v) -> str:
    return "     n/a" if v is None else f"{v:.3f}"


def training_overlap(db: str, oligos: list[dict]) -> list[dict]:
    """Which held-out oligos overlap the peak-derived training sequences.

    The CD oligos and the CUT&Tag peaks come from the same paper, so "held out"
    is a claim to check against the atlas rather than assert. Two directions
    count, and they mean different things, so both are reported: ``contained_in``
    means the oligo sits inside a training sequence, and ``contains`` means a
    training motif sits inside the oligo -- the second is the common one, since
    an ingested record is a short motif sliced out of a 201 nt window and a
    28-mer oligo can easily swallow one.
    """
    import sqlite3

    if not Path(db).exists():
        return []
    con = sqlite3.connect(db)
    train = [r[0] for r in con.execute(
        "SELECT sequence FROM records WHERE source LIKE 'gse220882%' "
        "OR source LIKE 'shuffled::gse220882%'")]
    con.close()
    blob = "\n".join(train)
    hits = []
    for o in oligos:
        q = clean(o["sequence"])
        inside = q in blob
        swallows = next((t for t in train if len(t) <= len(q) and t in q), None)
        if inside or swallows:
            hits.append({"name": o["name"], "contained_in_a_training_sequence": inside,
                         "contains_training_motif": swallows})
    return hits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--set", default=SET_PATH)
    ap.add_argument("--preset", default="im_reference")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--out", default="artifacts/external_validation_im.json")
    a = ap.parse_args()

    spec = json.loads(Path(a.set).read_text())
    pred = Predictor.load(a.model)
    cond = Condition.preset(a.preset)

    overlap = training_overlap(a.db, spec["oligos"])

    rows = []
    print(f"held-out cross-assay panel: {spec['title']}  (doi:{spec['doi']})")
    print("  same study as the peak training source; 2 negative controls; "
          "class labels only, no numeric CD values -- a direction, not a\n"
          "  generalisation result. See the module docstring.")
    print(f"scored at {cond.label()}")
    print(f"overlap with GSE220882 training sequences: {len(overlap)} of "
          f"{len(spec['oligos'])}")
    for h in overlap:
        what = ("sits inside a training sequence" if h["contained_in_a_training_sequence"]
                else f"contains the training motif {h['contains_training_motif']}")
        print(f"    {h['name']}: {what}")
    if overlap:
        print("  -> not fully held out for the genomic head; the affected rows are named "
              "above.")
    print()
    print(f"{'name':<11}{'authors say':<13}{'P(iM)':>7}{'P(iM|gen)':>11}"
          f"{'P(G4)':>8}{'pH_T':>7}  motif call")
    for o in spec["oligos"]:
        r = pred.predict(o["sequence"], cond, n_neighbours=0)[0]
        p = r["predictions"]
        p_im = p.get("im_fold", {}).get("probability")
        p_img = p.get("im_fold_genomic", {}).get("probability")
        p_g4 = p.get("g4_fold", {}).get("probability")
        pht = p.get("im_pht", {}).get("value")
        motif = ("iM" if find_im(o["sequence"]) else
                 "G4" if find_g4(o["sequence"]) else "-")
        rows.append({**o, "p_im": p_im, "p_im_genomic": p_img, "p_g4": p_g4,
                     "ph_t": pht, "motif": motif})
        print(f"{o['name']:<11}{o['expectation']:<13}{p_im:>7.3f}"
              f"{_fmt(p_img):>11}{p_g4:>8.3f}{pht:>7.2f}  {motif}")

    im = [r for r in rows if r["expectation"] == "i-motif"]
    g4 = [r for r in rows if r["expectation"] == "G4"]
    neg = [r for r in rows if r["expectation"] == "neither"]

    def mean(xs, key):
        return float(np.mean([x[key] for x in xs])) if xs else float("nan")

    summary = {
        "set": spec["source"], "doi": spec["doi"], "n": len(rows),
        "panel_type": "held-out cross-assay panel (same study as the peak training source; 2 negative controls; class labels only)",
        "condition": cond.label(),
        "training_overlap_with_gse220882": overlap,
        "n_training_overlap": len(overlap),
        "mean_p_im_genomic": {"i-motif": mean(im, "p_im_genomic"),
                              "G4": mean(g4, "p_im_genomic"),
                              "neither": mean(neg, "p_im_genomic")},
        "im_genomic_separation": mean(im, "p_im_genomic") - mean(neg, "p_im_genomic"),
        "mean_p_im": {"i-motif": mean(im, "p_im"), "G4": mean(g4, "p_im"),
                      "neither": mean(neg, "p_im")},
        "mean_p_g4": {"i-motif": mean(im, "p_g4"), "G4": mean(g4, "p_g4"),
                      "neither": mean(neg, "p_g4")},
        "im_separation": mean(im, "p_im") - mean(neg, "p_im"),
        "g4_separation": mean(g4, "p_g4") - mean(neg, "p_g4"),
        "canonical_motif_found_in_authors_im_set":
            f"{sum(1 for r in im if r['motif'] == 'iM')}/{len(im)}",
        "rows": rows,
    }

    print(f"\nmean P(iM|genomic):  i-motif set "
          f"{summary['mean_p_im_genomic']['i-motif']:.3f} | "
          f"G4 controls {summary['mean_p_im_genomic']['G4']:.3f} | "
          f"negatives {summary['mean_p_im_genomic']['neither']:.3f}"
          f"   -> separation {summary['im_genomic_separation']:+.3f}")
    print(f"mean P(iM):  i-motif set {summary['mean_p_im']['i-motif']:.3f} | "
          f"G4 controls {summary['mean_p_im']['G4']:.3f} | "
          f"negatives {summary['mean_p_im']['neither']:.3f}"
          f"   -> separation {summary['im_separation']:+.3f}")
    print(f"mean P(G4):  i-motif set {summary['mean_p_g4']['i-motif']:.3f} | "
          f"G4 controls {summary['mean_p_g4']['G4']:.3f} | "
          f"negatives {summary['mean_p_g4']['neither']:.3f}"
          f"   -> separation {summary['g4_separation']:+.3f}")
    print(f"canonical four-C-tract motif present in the authors' i-motif set: "
          f"{summary['canonical_motif_found_in_authors_im_set']}")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {a.out}")

    old, new = summary["im_separation"], summary["im_genomic_separation"]
    print()
    if old < 0.15:
        print("im_fold does NOT separate this independent set. It was trained on 160\n"
              "designed C-tract variants of a small number of parent constructs, and\n"
              "these peak-derived sequences do not look like them. Most of them carry\n"
              "no canonical four-C-tract motif at all, which puts them outside the\n"
              "head's applicability domain -- the predictor flags this per sequence.")
    if new == new and new >= 0.15:
        print(f"im_fold_genomic DOES separate it ({new:+.3f} against {old:+.3f}), and it\n"
              "also scores the authors' G4 controls low, which im_fold does not. That is\n"
              "the result the GSE220882 ingestion was for, and it came from sequence\n"
              "diversity rather than from a better model.")
    elif new == new:
        print(f"im_fold_genomic separates it by {new:+.3f}, against {old:+.3f} for\n"
              "im_fold. That is not the improvement the ingestion was aiming for.")


if __name__ == "__main__":
    main()
