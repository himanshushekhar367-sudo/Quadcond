#!/usr/bin/env python3
"""Test whether sequence clustering depends on the order rows arrive in.

compare_sidecar_vs_manifest.py shows that replaying the fold construction lands
on a different number of groups than the training run recorded, for four of the
five sequence-grouped heads. This asks why.

The suspect is in `quadcond.models.grouping._greedy_clusters`:

    order = np.argsort([-len(clean(s)) for s in seqs])

NumPy's default sort is quicksort, which is not stable. Many sequences share a
length, so tied rows come out in an arbitrary order that can differ between runs,
platforms and NumPy versions. The clustering is greedy and seeds a new group from
whichever tied sequence it reaches first, so a different traversal yields a
different number of groups from identical input.

Three checks, cheapest first:

1. **Same input twice.** If group counts differ, the function is not deterministic
   even within one process and the problem is worse than order-dependence.
2. **Shuffled input.** If shuffling the same sequences changes the group count,
   order-dependence is demonstrated, and any partition is only reproducible when
   the row order is also fixed and published.
3. **Stable sort.** Re-run with the ties broken deterministically to confirm the
   proposed fix removes the variation.

Needs numpy, so run it in an environment that can import quadcond -- the qcgw
conda environment or the checkout's .venv, not conda base.

    python check_determinism.py
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(REPO))


def load_sequences(atlas: Path, kind: str, limit: int) -> list[str]:
    con = sqlite3.connect(atlas)
    q = ("select distinct sequence from records where kind = ? and sequence is not null "
         "and tm is not null limit ?")
    seqs = [r[0] for r in con.execute(q, (kind, limit))]
    if not seqs:
        seqs = [r[0] for r in con.execute(
            "select distinct sequence from records where sequence is not null limit ?", (limit,))]
    con.close()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default=None, help="defaults to data/atlas.db, else atlas_core.db")
    ap.add_argument("--kind", default="g4")
    ap.add_argument("--limit", type=int, default=2274, help="g4_tm's row count by default")
    ap.add_argument("--threshold", type=float, default=0.90)
    a = ap.parse_args()

    atlas = Path(a.atlas) if a.atlas else None
    if atlas is None:
        full, core = REPO / "data" / "atlas.db", REPO / "data" / "atlas_core.db"
        atlas = full if full.exists() else core
    if not atlas.exists():
        raise SystemExit(f"no atlas at {atlas}")

    try:
        import numpy as np
        from quadcond.models.grouping import _greedy_clusters, cluster_sequences
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"cannot import what this needs: {exc}.\n"
            "Run it in an environment that can import quadcond -- the qcgw conda environment\n"
            "or the checkout's .venv. Conda base is not one.") from exc

    seqs = load_sequences(atlas, a.kind, a.limit)
    if len(seqs) < 50:
        raise SystemExit(f"only {len(seqs)} sequences read from {atlas.name}; nothing to test")

    def n_groups(s):
        return int(len(set(cluster_sequences(s, threshold=a.threshold).tolist())))

    first = n_groups(seqs)
    second = n_groups(seqs)

    rng = np.random.default_rng(0)
    shuffled = []
    for _ in range(3):
        idx = rng.permutation(len(seqs))
        shuffled.append(n_groups([seqs[i] for i in idx]))

    lengths = [len(s) for s in seqs]
    tied = len(lengths) - len(set(lengths))

    # A verdict of "not order-dependent" means nothing without knowing which
    # implementation produced it. The first run of this script reported
    # order_dependent=true; after the ordering was fixed it reports false, and a
    # reader with only the second file cannot tell that from "it was never
    # order-dependent". So the measured implementation is fingerprinted.
    import hashlib
    import inspect
    src = inspect.getsource(_greedy_clusters)
    impl = hashlib.sha256(src.encode()).hexdigest()[:16]
    content_ordered = "keys" in src and "sorted(range(" in src

    report = {
        "checked_on": date.today().isoformat(),
        "atlas": atlas.name,
        "grouping_impl_fingerprint": impl,
        "grouping_orders_by_content": content_ordered,
        "n_sequences": len(seqs),
        "distinct_lengths": len(set(lengths)),
        "rows_sharing_a_length_with_another": tied,
        "threshold": a.threshold,
        "groups_same_input_run_1": first,
        "groups_same_input_run_2": second,
        "deterministic_within_process": first == second,
        "groups_after_shuffling": shuffled,
        "order_dependent": len(set(shuffled + [first])) > 1,
    }
    reading = [
        ("Not deterministic even on identical input: the problem is more than tie ordering."
         if not report["deterministic_within_process"] else
         "Deterministic on identical input.")
    ]
    if report["order_dependent"]:
        reading.append(
            "Shuffling the same sequences changes the group count, so the partition depends on "
            "the order rows arrive in. A published partition would have to publish the row order "
            "too; ordering the traversal by sequence content is the better fix.")
    elif content_ordered:
        reading.append(
            "Shuffling does not change the group count, and this is the content-ordered "
            "implementation -- so this is the fix working, not evidence that order was never the "
            "problem. Compare against an earlier run of this script on the previous "
            "implementation before drawing any conclusion about the historical partition.")
    else:
        reading.append(
            "Shuffling does not change the group count on an implementation that does NOT order "
            "by content, so order-dependence is not the mechanism here and the gap against the "
            "training record has another cause. Look at whether the atlas rows or the clustering "
            "code changed.")
    reading.append(
        f"{tied} of {len(seqs)} sequences share a length with at least one other, which is the "
        "population an unstable sort can reorder.")
    reading.append(
        "This reads distinct sequences, so its group count is not comparable with the per-head "
        "counts in results/sidecar_vs_manifest.json, which are over a head's rows.")
    report["reading"] = reading

    # Keep every run, so the before/after pair survives.
    hist = ROOT / "results" / "determinism_history.jsonl"
    hist.parent.mkdir(exist_ok=True)
    with hist.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps({k: v for k, v in report.items() if k != "reading"}) + "\n")

    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    (out / "determinism_check.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("\nwritten: results/determinism_check.json")


if __name__ == "__main__":
    main()
