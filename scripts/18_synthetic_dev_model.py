#!/usr/bin/env python3
"""Train a throwaway model on generated data, for software testing only.

Why this file exists at all
---------------------------

A synthetic route is genuinely useful: it lets the service, the viewer and the
whole prediction path be exercised without a 400 MB atlas and an hour of
training. The engineering case for it is fine.

Why it is written the way it is
-------------------------------

A branch of this project carried a script that did the same job and did three
things that, together, made it the most dangerous file in the repository:

1. It wrote to ``artifacts/quadcond_model.joblib`` -- the canonical path the
   service loads. Running it replaced the measurement-trained model in place.
2. It registered its generated rows with ``evidence_tier="experimental"``, so
   ``claims.target_semantics`` classified every head as biophysically anchored.
   The model card, ``/info``, and every provenance badge in the viewer would
   have certified a fit to a generated CSV as measurement-grounded.
3. It assigned buffer conditions by ``random.choice`` from fixed lists, which
   gives each condition axis a wide ``n_unique`` -- so the heads would have been
   placed in the "Condition response" workspace and the interface would have
   drawn response surfaces through noise.

None of that was malicious and none of it was visible from the outside. That is
the point: the safeguards have to be in the code, not in remembering.

So this version:

* refuses to write anywhere that is not marked synthetic (``--out`` must contain
  ``synthetic``), and refuses the canonical filename outright;
* registers its sources at ``evidence_tier="synthetic"``, which
  ``claims.target_semantics`` classifies as ``SYNTHETIC`` -- so the headline
  reads "SYNTHETIC MODEL ... Not for scientific use", every head reports
  ``biophysically_grounded: false`` and ``has_experimental_observation: false``,
  and no amount of downstream formatting can undo it;
* stamps the model's atlas snapshot so ``quadcond.service`` refuses to serve it
  unless started with ``--allow-synthetic-model``, which prints a banner and
  flags every response.

    python3 scripts/18_synthetic_dev_model.py --csv dev/synthetic_dataset.csv \
        --out artifacts/synthetic_dev_model.joblib
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond.atlas import Atlas, Record            # noqa: E402
from quadcond.models.train import TaskSpec, train_all   # noqa: E402

#: The one path this script must never produce.
CANONICAL_NAMES = {"quadcond_model.joblib", "quadcond_model_seqonly.joblib"}

#: Everything trained here is tagged with this tier. `claims.target_semantics`
#: reads it and returns SYNTHETIC; nothing downstream can re-classify it.
SYNTHETIC_TIER = "synthetic"

KIND_MAP = {"g-quadruplex": "G4", "i-motif": "iM"}


def checked_output(raw: str) -> Path:
    out = Path(raw)
    if out.name in CANONICAL_NAMES:
        raise SystemExit(
            f"refusing to write {out.name}: that is the artifact the service "
            f"loads as the measurement-trained model. Choose a name containing "
            f"'synthetic'.")
    if "synthetic" not in out.name.lower():
        raise SystemExit(
            f"refusing to write {out}: a synthetic model's filename must say so, "
            f"because the filename is the only thing visible in a directory "
            f"listing six months from now.")
    return out


def build_atlas(csv_path: Path, db_path: Path, seed: int = 0) -> Atlas:
    import pandas as pd

    df = pd.read_csv(csv_path)
    if db_path.exists():
        db_path.unlink()
    atlas = Atlas(str(db_path))

    for key, note in (
        ("synthetic_dataset",
         "Sequences and labels generated for software testing. No measurement "
         "of any kind stands behind any row."),
        ("synthetic_negative",
         "Generated cross-class negatives. Not measured non-folders."),
    ):
        atlas.register_source(
            key,
            title=f"SYNTHETIC — {key} (software test fixture, not data)",
            evidence_tier=SYNTHETIC_TIER,
            notes=note,
        )

    # Conditions are drawn at random, and that is precisely why they must not be
    # allowed to look like a condition response. Recording them at all is a
    # convenience for exercising the featurizer; `capabilities.responds_to`
    # would otherwise see a wide n_unique and file these heads under "Condition
    # response", where the interface draws a response surface through noise.
    # The single fixed buffer below keeps n_unique at 1, so every axis reports
    # `state: "fixed"` and no surface is drawable.
    rng = random.Random(seed)
    fixed = {"k": 100.0, "na": 0.0, "li_nh4": 0.0, "mg": 0.0,
             "ph": 7.0, "temperature": 25.0, "crowder_pct": 0.0, "strand_conc": 5.0}

    records = []
    for idx, row in df.iterrows():
        kind = KIND_MAP.get(str(row["kind"]))
        if kind is None:
            continue
        records.append(Record(
            sequence=str(row["sequence"]),
            kind=kind,
            evidence_tier=SYNTHETIC_TIER,
            label_class="synthetic",
            folded=int(row.get("folded", 0)),
            topology=(str(row["topology"]) if str(row.get("topology", "")) != "nan" else None),
            tm=(float(row["tm"]) if str(row.get("tm", "")) not in ("", "nan") else None),
            method="generated",
            source="synthetic_dataset" if int(row.get("folded", 0)) else "synthetic_negative",
            organism=None,
            genomic=0,
            qc_flags="synthetic_fixture;no_measurement_upstream",
            **fixed,
        ))
    atlas.add_records(records)
    print(f"  {len(records):,} synthetic rows at one fixed buffer")
    return atlas


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--csv", required=True, help="generated dataset CSV")
    ap.add_argument("--out", default="artifacts/synthetic_dev_model.joblib")
    ap.add_argument("--db", default="/tmp/quadcond_synthetic_atlas.db")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = checked_output(args.out)
    csv = Path(args.csv)
    if not csv.exists():
        raise SystemExit(f"no such CSV: {csv}")

    print("SYNTHETIC MODEL BUILD — software fixture, not science.")
    atlas = build_atlas(csv, Path(args.db), args.seed)

    tasks = [
        TaskSpec(name="g4_fold", kind="G4", task="binary", target="folded",
                 sources=("synthetic_dataset", "synthetic_negative"),
                 claim="SYNTHETIC fixture head. Predicts a generated label."),
        TaskSpec(name="im_fold", kind="iM", task="binary", target="folded",
                 sources=("synthetic_dataset", "synthetic_negative"),
                 claim="SYNTHETIC fixture head. Predicts a generated label."),
        TaskSpec(name="g4_topology", kind="G4", task="multiclass", target="topology",
                 sources=("synthetic_dataset",),
                 claim="SYNTHETIC fixture head. Predicts a generated label."),
    ]
    model = train_all(atlas, tasks=tasks, use_conditions=False, n_seeds=2, n_folds=3)

    # The stamp the service checks. Belt as well as braces: the tier already
    # makes every head report SYNTHETIC, and this makes the whole artifact
    # refuse to load into a service that has not opted in.
    snap = dict(getattr(model, "atlas_snapshot", None) or {})
    snap["provenance"] = SYNTHETIC_TIER
    snap["provenance_note"] = (
        "Trained on generated data by scripts/18_synthetic_dev_model.py. Every "
        "number this model produces is a software test fixture. It is not "
        "evidence about DNA and must not be reported as a prediction.")
    model.atlas_snapshot = snap
    model.version = f"{getattr(model, 'version', '0')}+synthetic"

    out.parent.mkdir(parents=True, exist_ok=True)
    saved = model.save(out)
    print(f"\nSaved SYNTHETIC model to {saved}")
    print(f"  heads: {list(model.heads)}")
    print("  the service will refuse this artifact without --allow-synthetic-model")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
