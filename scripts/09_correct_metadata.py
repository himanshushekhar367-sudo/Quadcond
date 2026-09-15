#!/usr/bin/env python3
"""Cut a metadata-and-documentation correction without retraining anything.

v0.4.0's fitted estimators were fine. Its *descriptions* of them were not, in
ways that all pointed the same direction -- towards a reader over-reading the
two genomic-proxy heads:

* ``grounded_in_measurements=true`` on both proxy heads, because their positive
  rows are experimental. True, and misleading: what was measured is antibody
  occupancy at a locus, not folding. The model card duly printed
  "measurement-grounded: yes" while the results report called the same heads
  proxies. Three fields replace the one -- see :mod:`quadcond.claims`.
* Claim strings that said the proxy heads predict "P(sequence forms an i-motif
  in cells)". Antibody occupancy cannot license a folding probability.
* Regression metrics named ``conformal_halfwidth`` and ``empirical_coverage``.
  Both come from the same pool of out-of-fold residuals, so the coverage is
  in-sample to that pool and the word "conformal" promises a finite-sample
  guarantee that is not there.

This script rewrites those fields on the existing artifacts and bumps the
version. **The estimators, the calibrators and the atlas are untouched** --
every metric in the frozen report stays valid, and it stays valid *because*
nothing was refitted. The previous artifacts are copied aside first so the
v0.4.0 bytes remain checkable.

    python scripts/09_correct_metadata.py --version 0.4.1
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond import claims
from quadcond.models.base import MultiTaskModel

# metrics key -> corrected key. The values are unchanged; only the names lied.
METRIC_RENAMES = {
    "conformal_halfwidth": "oof_residual_halfwidth",
    "conformal_alpha": "oof_residual_alpha",
    "empirical_coverage": "in_sample_coverage_of_residual_pool",
}

INTERVAL_NOTE = (
    "Half-width is a percentile of the out-of-fold residual pool, and the "
    "coverage beside it was measured on those same residuals -- it is in-sample "
    "to the pool and is not a finite-sample coverage guarantee. Grouped split-"
    "conformal intervals are implemented in scripts/02_train.py for the next "
    "training run; this artifact predates them."
)

CALIBRATION_NOTE_KEY = "calibration_scope"


def correct(model: MultiTaskModel, version: str) -> list[str]:
    changes: list[str] = []
    for name, head in model.heads.items():
        meta, metrics = head.training_meta, head.metrics

        sem = claims.semantics(name, meta, head.task)
        before = meta.pop("grounded_in_measurements", None)
        meta.update({
            "has_experimental_observation": sem["has_experimental_observation"],
            "target_semantics": sem["target_semantics"],
            "target_semantics_label": sem["target_semantics_label"],
            "biophysically_grounded": sem["biophysically_grounded"],
        })
        if before is not None and before != sem["biophysically_grounded"]:
            changes.append(
                f"{name}: grounded_in_measurements={before} -> "
                f"biophysically_grounded={sem['biophysically_grounded']} "
                f"(target_semantics={sem['target_semantics']})"
            )

        metrics[CALIBRATION_NOTE_KEY] = sem["calibration_scope"]

        for old, new in METRIC_RENAMES.items():
            if old in metrics:
                metrics[new] = metrics.pop(old)
                changes.append(f"{name}: metrics.{old} -> metrics.{new}")
        if "oof_residual_halfwidth" in metrics:
            metrics["interval_note"] = INTERVAL_NOTE

        corrected = claims.CORRECTED_CLAIMS.get(name)
        if corrected and meta.get("claim") != corrected:
            meta["claim"] = corrected
            changes.append(f"{name}: claim rewritten (folding probability -> proxy score)")

    model.version = version
    return changes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="0.4.1")
    ap.add_argument("--artifacts", nargs="*", default=[
        "artifacts/quadcond_model.joblib",
        "artifacts/quadcond_model_seqonly.joblib",
    ])
    ap.add_argument("--backup-suffix", default=".v0.4.0")
    ap.add_argument("--no-backup", action="store_true")
    a = ap.parse_args()

    for path in a.artifacts:
        p = Path(path)
        if not p.exists():
            print(f"  ! {p} not found; skipping")
            continue
        if not a.no_backup:
            keep = p.with_suffix(p.suffix + a.backup_suffix)
            if not keep.exists():
                shutil.copy2(p, keep)
                print(f"  kept the previous bytes at {keep}")

        model = MultiTaskModel.load(p, verify=True)
        changes = correct(model, a.version)
        model.save(p)
        print(f"\n{p} -> v{model.version}  sha256 {model.artifact_sha256}")
        print(f"  {claims.headline(model.heads)}")
        for c in changes:
            print(f"    {c}")

    print("\nNothing was refitted: estimators, calibrators and the atlas are byte-identical "
          "to v0.4.0, so every metric in the frozen report still describes these models.")


if __name__ == "__main__":
    main()
