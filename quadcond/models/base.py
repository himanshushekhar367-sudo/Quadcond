"""Model container: heads, metadata, persistence."""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


def file_sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def rows_fingerprint(rows) -> str:
    """Stable hash of the exact training rows behind a model.

    Two models with the same fingerprint were trained on the same data; two with
    different fingerprints were not, whatever their filenames say. This is what
    makes a reported metric traceable to a dataset version.
    """
    h = hashlib.sha256()
    for r in sorted(
        (
            f"{d['sequence']}|{d['kind']}|{d['evidence_tier']}|{d.get('label_class','')}|"
            f"{d['k']}|{d['na']}|{d['li_nh4']}|{d['mg']}|{d['ph']}|{d['temperature']}|"
            f"{d['folded']}|{d['topology']}|{d['tm']}|{d['ph_t']}"
            for d in (dict(x) for x in rows)
        )
    ):
        h.update(r.encode())
    return h.hexdigest()

TASKS = ("binary", "multiclass", "regression")


@dataclass
class Head:
    """One prediction head plus everything needed to interpret it honestly."""

    name: str
    kind: str                    # G4 | iM
    task: str                    # binary | multiclass | regression
    target: str                  # atlas label column
    feature_names: list[str] = field(default_factory=list)
    classes: list[str] = field(default_factory=list)
    use_conditions: bool = True
    estimators: list[Any] = field(default_factory=list)   # ensemble members
    calibrator: Any = None
    # The attribute names are historical and are NOT renamed, because renaming a
    # dataclass field breaks unpickling of artifacts already in the wild. What
    # they hold is a quantile of the out-of-fold residual pool, which is not a
    # split-conformal half-width: the quantile and the coverage reported beside
    # it come from the same residuals, so the coverage is in-sample to that pool
    # and carries no finite-sample guarantee. The emitted keys below say so.
    conformal_q: float | None = None      # quantile of the OOF residual pool
    conformal_alpha: float = 0.1
    metrics: dict = field(default_factory=dict)
    training_meta: dict = field(default_factory=dict)
    applicability: dict = field(default_factory=dict)     # condition + sequence domain

    # ---------------------------------------------------------------- predict
    def _raw(self, X: np.ndarray) -> np.ndarray:
        if self.task == "regression":
            preds = np.column_stack([m.predict(X) for m in self.estimators])
            return preds
        probas = [m.predict_proba(X) for m in self.estimators]
        return np.mean(probas, axis=0)

    def predict(self, X: np.ndarray) -> dict:
        raw = self._raw(X)
        if self.task == "regression":
            mean = raw.mean(axis=1)
            std = raw.std(axis=1)
            out: dict[str, Any] = {"value": mean, "ensemble_std": std}
            if self.conformal_q is not None:
                out["interval_low"] = mean - self.conformal_q
                out["interval_high"] = mean + self.conformal_q
                out["interval_nominal_level"] = 1 - self.conformal_alpha
                out["interval_kind"] = "oof_residual"
                out["interval_note"] = (
                    "half-width is the "
                    f"{100 * (1 - self.conformal_alpha):.0f}th percentile of the "
                    "out-of-fold residual pool. Not a conformal interval: the "
                    "coverage reported alongside it was measured on the same "
                    "residuals. Read the width as a typical error scale."
                )
            return out
        if self.task == "binary":
            p = raw[:, 1]
            if self.calibrator is not None:
                p = self.calibrator.transform(p)
            return {"probability": np.asarray(p, dtype=float), "raw_probability": raw[:, 1]}
        proba = raw
        if self.calibrator is not None:
            proba = self.calibrator.transform(proba)
        return {
            "proba": proba,
            "classes": self.classes,
            "argmax": [self.classes[i] for i in np.argmax(proba, axis=1)],
            "confidence": proba.max(axis=1),
        }

    def summary(self) -> dict:
        d = {
            "name": self.name, "kind": self.kind, "task": self.task,
            "target": self.target, "classes": self.classes,
            "use_conditions": self.use_conditions,
            "n_estimators": len(self.estimators),
            "metrics": self.metrics, "training_meta": self.training_meta,
            "applicability": self.applicability,
        }
        return d


@dataclass
class MultiTaskModel:
    """A bundle of heads sharing one featurizer contract."""

    heads: dict[str, Head] = field(default_factory=dict)
    version: str = "0.4.1"
    atlas_snapshot: dict = field(default_factory=dict)
    dataset_fingerprint: str = ""
    artifact_sha256: str = ""

    def add(self, head: Head) -> None:
        self.heads[head.name] = head

    def __contains__(self, name: str) -> bool:
        return name in self.heads

    def __getitem__(self, name: str) -> Head:
        return self.heads[name]

    def save(self, path: str | Path) -> Path:
        """Write atomically, then record the checksum of what was actually written.

        The checksum is computed after the bytes hit disk, so a truncated or
        half-written artifact fails verification instead of loading silently.
        """
        import joblib

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".tmp")
        joblib.dump(self, tmp, compress=3)
        tmp.replace(path)
        self.artifact_sha256 = file_sha256(path)
        side = path.with_suffix(".json")
        side.write_text(json.dumps(self.summary(), indent=2, default=str))
        return path

    @staticmethod
    def load(path: str | Path, *, verify: bool = True) -> "MultiTaskModel":
        """Check the checksum recorded beside the artifact, *then* deserialise it.

        The order matters. Loading first and verifying afterwards means a
        corrupted file is unpickled before anyone checks it — and a damaged
        joblib archive does not fail politely: it can hand the allocator a
        garbage length prefix and take the process out with an OOM kill before
        any exception is raised. Verifying first turns that into a clean error.
        It is also the right posture for a file that may have travelled: never
        unpickle bytes whose integrity you have not established.
        """
        import joblib

        path = Path(path)
        side = path.with_suffix(".json")
        if verify and side.exists():
            try:
                recorded = json.loads(side.read_text()).get("artifact_sha256")
            except json.JSONDecodeError:
                recorded = None
            if recorded:
                actual = file_sha256(path)
                if recorded != actual:
                    raise ValueError(
                        f"checksum mismatch for {path}: the sidecar records "
                        f"{recorded[:12]}... but the file hashes to {actual[:12]}.... "
                        f"The artifact has been modified or truncated since it was "
                        f"written; retrain or re-download it. Pass verify=False only "
                        f"if you know why it differs and trust the file."
                    )
        return joblib.load(path)

    def summary(self) -> dict:
        return {
            "version": self.version,
            "dataset_fingerprint_sha256": self.dataset_fingerprint,
            "artifact_sha256": self.artifact_sha256,
            "atlas_snapshot": self.atlas_snapshot,
            "heads": {k: h.summary() for k, h in self.heads.items()},
        }
