"""Benchmarks and diagnostics that go beyond a single accuracy number."""
from __future__ import annotations

import numpy as np

from . import readout
from .conditions import Condition
from .models.calibration import expected_calibration_error, reliability_curve
from .motifs import find_g4, find_im, g4hunter_mean, g4hunter_window_max


def baseline_scores(sequences, kind: str = "G4") -> dict[str, np.ndarray]:
    """Rule-based baselines every new G4/iM model should have to beat."""
    seqs = list(sequences)
    sign = 1.0 if kind == "G4" else -1.0
    return {
        "g4hunter_mean": np.array([sign * g4hunter_mean(s) for s in seqs]),
        "g4hunter_window_max": np.array([g4hunter_window_max(s) for s in seqs]),
        "canonical_motif": np.array(
            [1.0 if (find_g4(s) if kind == "G4" else find_im(s)) else 0.0 for s in seqs]
        ),
        "gc_content": np.array(
            [(s.count("G") + s.count("C")) / max(len(s), 1) for s in seqs]
        ),
    }


def compare_to_baselines(atlas, head, *, kind: str = "G4", target: str = "folded",
                         tiers=("experimental", "derived")) -> dict:
    """Ranking *and* calibration comparison against rule-based scores.

    A rule-based score can rank well (high AUROC) while being useless as a
    probability. Reporting both columns is the point: it shows what calibration
    actually buys.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score

    rows = atlas.query(kind=kind, tiers=tiers, label=target)
    seqs = [r["sequence"] for r in rows]
    y = np.array([int(r[target]) for r in rows])
    if len(np.unique(y)) < 2:
        return {"error": "need both classes to benchmark"}

    out: dict[str, dict] = {}
    for name, score in baseline_scores(seqs, kind).items():
        rng = score.max() - score.min()
        p = (score - score.min()) / rng if rng else np.full_like(score, 0.5)
        out[name] = {
            "auroc": float(roc_auc_score(y, score)),
            "auprc": float(average_precision_score(y, score)),
            "ece_if_used_as_probability": expected_calibration_error(y, p),
            "note": "min-max rescaled to [0,1]; rule-based scores are not probabilities",
        }
    out["quadcond_" + head.name] = {
        "auroc": head.metrics.get("auroc"),
        "auprc": head.metrics.get("auprc"),
        "ece_if_used_as_probability": head.metrics.get("ece"),
        "note": "calibrated, grouped out-of-fold",
    }
    return out


def condition_response(predictor, sequence: str, variable: str,
                       values, base: Condition | None = None) -> list[dict]:
    """How every head responds along one condition axis."""
    base = base or Condition()
    rows = []
    for v in values:
        cond = base.replace(**{variable: float(v)})
        r = predictor.predict(sequence, cond, n_neighbours=0)[0]
        row = {variable: float(v)}
        for name, p in r["predictions"].items():
            if readout.is_refused(p):
                row[f"{name}:refused"] = readout.refusal_reason(p)
                continue
            if "probability" in p:
                row[name] = p["probability"]
            elif "posterior" in p:
                row.update({f"{name}:{k}": val for k, val in p["posterior"].items()})
            else:
                row[name] = p["value"]
                if "folded_fraction_at_condition" in p:
                    row[f"{name}:folded_fraction"] = p["folded_fraction_at_condition"]
        row["in_domain"] = all(p["applicability"]["in_domain"]
                               for p in r["predictions"].values())
        rows.append(row)
    return rows


def response_sensitivity(rows: list[dict], variable: str) -> dict[str, float]:
    """Peak-to-peak change of each output across a condition sweep.

    Near-zero means the head is insensitive to that variable -- either because
    the biology says so, or (far more often) because the training data never
    varied it. The model card's applicability table tells you which.
    """
    keys = [k for k in rows[0] if k not in {variable, "in_domain"}]
    return {k: float(max(r[k] for r in rows) - min(r[k] for r in rows)) for k in keys}
