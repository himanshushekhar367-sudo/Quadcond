#!/usr/bin/env python3
"""Train every head the atlas can support, with and without conditions.

The ablation is not optional: the whole claim of QuadCond is that explicit
condition features matter, so the paired sequence-only run is produced by
default and both sets of metrics go into the model card.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond.atlas import Atlas
from quadcond.models.train import DEFAULT_TASKS, train_all


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--out", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--ablation-out", default="artifacts/quadcond_model_seqonly.joblib")
    ap.add_argument("--no-ablation", action="store_true")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--tasks", nargs="*", default=None)
    args = ap.parse_args()

    tasks = DEFAULT_TASKS
    if args.tasks:
        tasks = [t for t in DEFAULT_TASKS if t.name in set(args.tasks)]

    atlas = Atlas(args.db)
    print(f"Atlas: {atlas.composition()['one_line']}\n")

    print("== full model (sequence + conditions) ==")
    full = train_all(atlas, tasks=tasks, use_conditions=True,
                     n_seeds=args.seeds, n_folds=args.folds)
    p = full.save(args.out)
    print(f"  saved -> {p}")

    if not args.no_ablation:
        print("\n== ablation (sequence only) ==")
        seq = train_all(atlas, tasks=tasks, use_conditions=False,
                        n_seeds=args.seeds, n_folds=args.folds)
        p2 = seq.save(args.ablation_out)
        print(f"  saved -> {p2}")

        rows = []
        for name, h in full.heads.items():
            if name not in seq.heads:
                continue
            key = {"binary": "auroc", "multiclass": "balanced_accuracy",
                   "regression": "r2"}[h.task]
            rows.append({
                "head": name, "metric": key,
                "with_conditions": h.metrics.get(key),
                "sequence_only": seq.heads[name].metrics.get(key),
                "delta": (h.metrics.get(key) or 0) - (seq.heads[name].metrics.get(key) or 0),
                "ece_with_conditions": h.metrics.get("ece"),
                "ece_sequence_only": seq.heads[name].metrics.get("ece"),
            })
        Path("artifacts").mkdir(exist_ok=True)
        Path("artifacts/ablation.json").write_text(json.dumps(rows, indent=2))
        print("\nAblation (condition features on vs off)")
        for r in rows:
            print(f"  {r['head']:<18} {r['metric']:<18} "
                  f"{r['with_conditions']:.4f} vs {r['sequence_only']:.4f} "
                  f"(delta {r['delta']:+.4f})")

    atlas.close()


if __name__ == "__main__":
    main()
