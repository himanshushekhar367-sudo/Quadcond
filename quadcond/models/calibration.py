"""Probability calibration and calibration diagnostics.

A predictor that says "0.9" should be right about 90% of the time.  Published
G4/iM classifiers report AUROC and accuracy, which are invariant to any
monotone rescaling of the score -- so a tool can rank perfectly and still be
badly miscalibrated, which is precisely what burns you when you use the score
as a probability to prioritise candidates.  QuadCond therefore treats
calibration as a first-class output and reports it for every head.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def expected_calibration_error(y_true, p, n_bins: int = 10, strategy: str = "quantile") -> float:
    """ECE: |accuracy - confidence| averaged over bins, weighted by bin size."""
    y_true = np.asarray(y_true).astype(float)
    p = np.asarray(p, dtype=float)
    if strategy == "quantile":
        edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
    else:
        edges = np.linspace(0, 1, n_bins + 1)
    if len(edges) < 2:
        return 0.0
    idx = np.clip(np.digitize(p, edges[1:-1], right=True), 0, len(edges) - 2)
    ece = 0.0
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        ece += m.mean() * abs(y_true[m].mean() - p[m].mean())
    return float(ece)


def maximum_calibration_error(y_true, p, n_bins: int = 10) -> float:
    y_true = np.asarray(y_true).astype(float)
    p = np.asarray(p, dtype=float)
    edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 2:
        return 0.0
    idx = np.clip(np.digitize(p, edges[1:-1], right=True), 0, len(edges) - 2)
    worst = 0.0
    for b in range(len(edges) - 1):
        m = idx == b
        if m.any():
            worst = max(worst, abs(y_true[m].mean() - p[m].mean()))
    return float(worst)


def reliability_curve(y_true, p, n_bins: int = 10):
    """(bin confidence, bin accuracy, bin count) for a reliability diagram."""
    y_true = np.asarray(y_true).astype(float)
    p = np.asarray(p, dtype=float)
    edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 2:
        return np.array([]), np.array([]), np.array([])
    idx = np.clip(np.digitize(p, edges[1:-1], right=True), 0, len(edges) - 2)
    conf, acc, cnt = [], [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        conf.append(p[m].mean())
        acc.append(y_true[m].mean())
        cnt.append(int(m.sum()))
    return np.array(conf), np.array(acc), np.array(cnt)


def multiclass_ece(y_true, proba, n_bins: int = 10) -> float:
    """Top-label ECE for a multiclass head."""
    proba = np.asarray(proba, dtype=float)
    conf = proba.max(axis=1)
    pred = proba.argmax(axis=1)
    correct = (pred == np.asarray(y_true)).astype(float)
    return expected_calibration_error(correct, conf, n_bins=n_bins)


def brier_multiclass(y_true, proba) -> float:
    proba = np.asarray(proba, dtype=float)
    onehot = np.zeros_like(proba)
    onehot[np.arange(len(y_true)), np.asarray(y_true)] = 1.0
    return float(np.mean(np.sum((proba - onehot) ** 2, axis=1)))


# --------------------------------------------------------------------------- #
# Calibrators
# --------------------------------------------------------------------------- #
@dataclass
class BinaryCalibrator:
    """Isotonic or Platt rescaling for a binary head."""

    method: str = "isotonic"
    _iso: IsotonicRegression | None = None
    _lr: LogisticRegression | None = None

    def fit(self, p, y) -> "BinaryCalibrator":
        p = np.asarray(p, dtype=float).ravel()
        y = np.asarray(y).astype(int).ravel()
        if self.method == "isotonic":
            self._iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
            self._iso.fit(p, y)
        elif self.method == "platt":
            self._lr = LogisticRegression(C=1e6, solver="lbfgs")
            self._lr.fit(_logit(p).reshape(-1, 1), y)
        elif self.method == "none":
            pass
        else:
            raise ValueError(f"unknown calibration method {self.method!r}")
        return self

    def transform(self, p):
        p = np.asarray(p, dtype=float).ravel()
        if self.method == "isotonic" and self._iso is not None:
            return np.clip(self._iso.predict(p), 1e-6, 1 - 1e-6)
        if self.method == "platt" and self._lr is not None:
            return self._lr.predict_proba(_logit(p).reshape(-1, 1))[:, 1]
        return p


@dataclass
class TemperatureScaler:
    """Single-parameter temperature scaling for a multiclass head.

    Divides the logits by T > 0 chosen to minimise NLL on held-out data.  It
    cannot change the argmax, so accuracy is untouched -- it only fixes
    over/under-confidence, which is exactly the failure mode we care about.
    """

    temperature: float = 1.0

    def fit(self, proba, y) -> "TemperatureScaler":
        logits = np.log(np.clip(np.asarray(proba, dtype=float), 1e-12, 1.0))
        y = np.asarray(y).astype(int)

        def nll(logT: float) -> float:
            T = float(np.exp(logT))
            z = logits / T
            z = z - z.max(axis=1, keepdims=True)
            logp = z - np.log(np.exp(z).sum(axis=1, keepdims=True))
            return float(-logp[np.arange(len(y)), y].mean())

        res = minimize_scalar(nll, bounds=(-3.0, 3.0), method="bounded")
        self.temperature = float(np.exp(res.x))
        return self

    def transform(self, proba):
        logits = np.log(np.clip(np.asarray(proba, dtype=float), 1e-12, 1.0)) / self.temperature
        logits = logits - logits.max(axis=1, keepdims=True)
        e = np.exp(logits)
        return e / e.sum(axis=1, keepdims=True)


def _logit(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))
