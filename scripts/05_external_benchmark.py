#!/usr/bin/env python3
"""Score an external G4/i-motif model against the atlas's real experimental records.

Any predictor can be evaluated here as long as you can write a five-line adapter
that turns a sequence plus a condition into a topology call and a folding
probability. The point is to make "does this model work on measurements?" a
routine question with a routine answer, rather than a claim.

Two tasks are run:

**Topology** on every experimental record in the atlas that carries a topology
label, against the always-majority-class and chance baselines. A model can beat
chance on accuracy purely by predicting the plurality class, so balanced
accuracy and the confusion matrix are the numbers that matter.

**Folding** on those same experimental positives against composition-matched
dinucleotide shuffles of them, scored for ranking (AUROC) *and* for calibration
(Brier, ECE), because a score that ranks well and is miscalibrated is still
unusable as a probability.

    python scripts/05_external_benchmark.py --adapter quadcond_external \\
        --model-path path/to/demo_model.joblib --src path/to/their/src

**Leakage is checked, not assumed.** Before reporting anything, the script
reconstructs each head's training set from the metadata the model stores (its
sources, tiers and label classes) and measures how much of the benchmark set the
model has already seen. Any overlap at all and the result is reported as
IN-SAMPLE, with the model's own grouped out-of-fold metrics printed alongside as
the honest number. Scoring a model on its own training data produces a beautiful
table and means nothing; this package should not be able to do it by accident
any more than anyone else's.

Writing an adapter: add an entry to ADAPTERS below returning a callable
``(sequence) -> (topology_label_or_None, p_folded_or_None)``.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.metrics import (balanced_accuracy_score, brier_score_loss,
                             confusion_matrix, roc_auc_score)

TOPOLOGIES = ["parallel", "antiparallel", "hybrid"]


# --------------------------------------------------------------------------- #
# Adapters
# --------------------------------------------------------------------------- #
def adapter_quadcond_self(args):
    """This package's own model, so the comparison has a like-for-like column."""
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor

    pred = Predictor.load(args.model_path or "artifacts/quadcond_model.joblib")
    cond = Condition.preset(args.preset)

    def score(seq: str):
        r = pred.predict(seq, cond, n_neighbours=0)[0]["predictions"]
        topo = r.get("g4_topology", {}).get("argmax")
        p = r.get("g4_fold", {}).get("probability")
        return topo, p

    return score, "quadcond (this package)", pred.model


def adapter_quadcond_external(args):
    """A ModelBundle from an independently authored QuadCond package.

    Requires ``--src`` pointing at that project's ``src`` directory. Its package
    is also called ``quadcond``, so it is imported in a subprocess-free way by
    putting its path first and importing the bundle before anything else.
    """
    if not args.src:
        raise SystemExit("--src is required for the quadcond_external adapter")
    sys.path.insert(0, args.src)
    from quadcond.model import ModelBundle  # noqa: E402  (theirs, not ours)
    from quadcond.schema import PredictionConditions  # noqa: E402

    bundle = ModelBundle.load(args.model_path)
    cond = PredictionConditions(k_mm=100.0, na_mm=0.0, mg_mm=0.0, ph=7.0,
                                temperature_c=25.0)

    def score(seq: str):
        r = bundle.predict(seq, "DNA", cond)
        g4 = r.get("g4") or {}
        return ((g4.get("topology") or {}).get("predicted_class"),
                (g4.get("folding") or {}).get("probability"))

    return score, f"{bundle.name} (provenance={bundle.provenance})", None


ADAPTERS = {"quadcond_self": adapter_quadcond_self, "quadcond_external": adapter_quadcond_external}


def training_overlap(model, atlas, benchmark_sequences) -> dict:
    """How much of the benchmark set did this model already train on?

    Rebuilt from the head metadata the model stores, so it works for any model
    this package produced without keeping a copy of the training sequences in
    the artifact.
    """
    from quadcond.motifs import clean

    if model is None:
        return {"checkable": False,
                "note": "external model: no training-set metadata available, "
                        "assumed held out"}
    bench = {clean(s) for s in benchmark_sequences}
    per_head = {}
    for name, head in model.heads.items():
        meta = head.training_meta
        rows = atlas.query(
            kind=head.kind,
            tiers=tuple(meta.get("tiers") or ()) or None,
            sources=tuple(meta.get("sources") or ()) or None,
            label_classes=tuple(meta.get("label_classes") or ()) or None,
            label=head.target,
        )
        seen = {clean(r["sequence"]) for r in rows}
        hit = bench & seen
        per_head[name] = {
            "training_rows_matched": len(rows),
            "benchmark_sequences_seen": len(hit),
            "fraction_seen": round(len(hit) / max(len(bench), 1), 4),
            "out_of_fold_metric": head.metrics.get(
                {"binary": "auroc", "multiclass": "balanced_accuracy",
                 "regression": "r2"}[head.task]
            ),
        }
    worst = max((h["fraction_seen"] for h in per_head.values()), default=0.0)
    return {"checkable": True, "max_fraction_seen": worst,
            "in_sample": worst > 0.0, "per_head": per_head}


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--adapter", choices=sorted(ADAPTERS), default="quadcond_self")
    ap.add_argument("--model-path", default=None)
    ap.add_argument("--src", default=None, help="external package source dir, if needed")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--preset", default="k100")
    ap.add_argument("--out", default="artifacts/external_benchmark.json")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    # Pull the real experimental records before any external package shadows ours.
    from quadcond.atlas import Atlas
    from quadcond.models.calibration import expected_calibration_error
    from quadcond.negatives import matched_negatives

    atlas = Atlas(a.db)
    rows = atlas.query(kind="G4", tiers=("experimental",),
                       label_classes=("biophysical",), label="topology")
    atlas.close()
    if not rows:
        raise SystemExit(f"no experimental topology records in {a.db}")
    seqs = [r["sequence"] for r in rows][: a.limit]
    truth = [r["topology"] for r in rows][: a.limit]
    negatives = matched_negatives(seqs, kind="G4", n_per_positive=1, seed=7)

    score, label, own_model = ADAPTERS[a.adapter](a)
    print(f"benchmarking: {label}")
    print(f"real experimental records: {len(seqs)}")

    atlas = Atlas(a.db)
    overlap = training_overlap(own_model, atlas, seqs)
    atlas.close()
    if overlap.get("in_sample"):
        print("\n" + "=" * 72)
        print("IN-SAMPLE: this model was trained on these records. The numbers below")
        print("are memorisation, not generalisation, and are NOT comparable to a model")
        print("that has never seen this data. The honest figures are the grouped")
        print("out-of-fold metrics the model reports for itself:")
        for name, h in overlap["per_head"].items():
            if h["fraction_seen"] > 0:
                print(f"    {name:<18} out-of-fold = {h['out_of_fold_metric']:.4f}"
                      f"   ({h['fraction_seen']:.0%} of the benchmark set was in training)")
        print("=" * 72)
    elif not overlap.get("checkable"):
        print(f"note: {overlap['note']}")
    print()

    calls, probs = [], []
    for s in seqs:
        t, p = score(s)
        calls.append(t)
        probs.append(p)

    df = pd.DataFrame({"truth": truth, "call": calls}).dropna(subset=["call"])
    result: dict = {"model": label, "n_records": len(seqs),
                    "training_overlap": overlap,
                    "held_out": not overlap.get("in_sample", False)}
    if len(df):
        acc = float((df["call"] == df["truth"]).mean())
        bal = float(balanced_accuracy_score(df["truth"], df["call"]))
        maj = df["truth"].value_counts().idxmax()
        base = float((df["truth"] == maj).mean())
        print(f"TOPOLOGY   accuracy {acc:.3f} | balanced accuracy {bal:.3f}  (n={len(df)})")
        print(pd.DataFrame(confusion_matrix(df["truth"], df["call"], labels=TOPOLOGIES),
                           index=[f"true {t}" for t in TOPOLOGIES],
                           columns=[f"pred {t}" for t in TOPOLOGIES]).to_string())
        print(f"baselines: always-'{maj}' accuracy {base:.3f} | chance balanced accuracy 0.333\n")
        result["topology"] = {"accuracy": acc, "balanced_accuracy": bal,
                              "n": int(len(df)), "majority_baseline_accuracy": base,
                              "chance_balanced_accuracy": 1 / 3}
    else:
        print("TOPOLOGY   the model returned no topology calls\n")

    y, p = [], []
    for group, lab in ((seqs, 1), (negatives, 0)):
        for s in group:
            _, pr = score(s)
            if pr is not None:
                y.append(lab)
                p.append(float(pr))
    if y:
        y, p = np.array(y), np.array(p)
        auroc = float(roc_auc_score(y, p))
        brier = float(brier_score_loss(y, p))
        ece = float(expected_calibration_error(y, p))
        print(f"FOLDING    real G4s vs composition-matched shuffles (n={len(y)})")
        print(f"  AUROC {auroc:.3f} | Brier {brier:.3f} | ECE {ece:.3f}")
        print(f"  mean p(real G4) {p[y == 1].mean():.3f} vs mean p(shuffle) {p[y == 0].mean():.3f}")
        result["folding"] = {"auroc": auroc, "brier": brier, "ece": ece, "n": int(len(y)),
                             "mean_p_positive": float(p[y == 1].mean()),
                             "mean_p_negative": float(p[y == 0].mean())}

    if overlap.get("in_sample"):
        result["warning"] = (
            "IN-SAMPLE: the model was trained on these records. Use the per-head "
            "out_of_fold_metric values in training_overlap for any comparison."
        )
        result["reportable_metrics"] = {
            name: h["out_of_fold_metric"]
            for name, h in overlap["per_head"].items() if h["fraction_seen"] > 0
        }
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(result, indent=2))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
