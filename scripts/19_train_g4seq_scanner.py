#!/usr/bin/env python3
"""Train the genomic G4 window scanner on G4-seq (K+) and evaluate it honestly.

Input: a table of 124-nt windows with columns id, sequence, y, negtype, split,
built from the G4detector benchmark release of G4-seq (GSE110582):
positives = observed K+ G4-seq windows; negatives = random genomic windows and
canonical-motif windows G4-seq did not observe. Chromosomes 2, 8 and 17 are the
held-out test; the mouse set is a cross-species test.

Output: artifacts/quadcond_g4seq_scanner.joblib and a metrics JSON.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

from quadcond.genome_scan import G4SeqScanner, window_features
from quadcond.models.calibration import BinaryCalibrator, expected_calibration_error
from quadcond.models.train import _estimator


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--out", default="artifacts/quadcond_g4seq_scanner.joblib")
    ap.add_argument("--metrics", default="artifacts/g4seq_scanner_metrics.json")
    ap.add_argument("--predictions", default="artifacts/g4seq_scanner_test_predictions.csv.gz")
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()
    t0 = time.time()
    d = pd.read_csv(a.windows)
    tr = d[d.split == "train"].reset_index(drop=True)
    X = window_features(tr.sequence.tolist())
    y = tr.y.values
    chrom = tr.id.str.split(":").str[0].values
    print(f"train {X.shape} in {time.time()-t0:.0f}s", flush=True)

    # Chromosome-grouped out-of-fold predictions, used only to fit the calibrator
    # and to report internal metrics.
    oof = np.zeros(len(y))
    for k, (i, j) in enumerate(GroupKFold(n_splits=5).split(X, y, chrom)):
        m = _estimator("binary", k, len(i)).fit(X[i], y[i])
        oof[j] = m.predict_proba(X[j])[:, 1]
        print(f"  fold {k}", flush=True)
    cal = BinaryCalibrator("isotonic").fit(oof, y)
    est = [_estimator("binary", 100 + s, len(y)).fit(X, y) for s in range(a.seeds)]
    sc = G4SeqScanner(estimators=est, calibrator=cal)

    metrics = {"internal_chromosome_grouped_cv": {
        "auroc": float(roc_auc_score(y, oof)), "auprc": float(average_precision_score(y, oof))}}
    preds = []
    for split in ["test_heldout_chr", "test_mouse"]:
        te = d[d.split == split].reset_index(drop=True)
        p = sc.score_windows(te.sequence.tolist())
        te = te.assign(score=p)
        preds.append(te)
        res = {}
        for neg in ["random", "pq"]:
            s = te[te.negtype.isin(["pos", neg])]
            res[f"vs_{neg}"] = {"auroc": float(roc_auc_score(s.y, s.score)),
                                "auprc": float(average_precision_score(s.y, s.score)),
                                "n_pos": int(s.y.sum()), "n_neg": int((s.y == 0).sum())}
        res["all"] = {"auroc": float(roc_auc_score(te.y, te.score)),
                      "brier": float(brier_score_loss(te.y, te.score)),
                      "ece": float(expected_calibration_error(te.y.values, te.score.values))}
        metrics[split] = res
        print(split, json.dumps(res), flush=True)
    sc.metrics = metrics
    sc.training_meta = {
        "source": "G4-seq K+ (GSE110582) via the G4detector benchmark release",
        "n_train": int(len(y)), "n_pos": int(y.sum()),
        "negatives": "random genomic windows + canonical-motif windows not observed by G4-seq",
        "held_out_chromosomes": ["chr2", "chr8", "chr17"], "cross_species_test": "mouse G4-seq K+",
        "features": int(X.shape[1]), "seeds": a.seeds, "calibration": "isotonic on chromosome-grouped OOF",
        "seconds": round(time.time() - t0, 1),
    }
    sc.save(a.out)
    Path(a.metrics).write_text(json.dumps(sc.describe(), indent=2))
    pd.concat(preds).to_csv(a.predictions, index=False)
    print("saved", a.out, round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
