#!/usr/bin/env python3
"""Is the shuffled-negative generation reproducible?

It matters because `g4_fold` and `im_fold` train on those negatives. If they are
not reproducible, publishing the atlas does not make those heads reproducible
either, and retraining would rebuild the very problem the fold work was meant to
remove.

Two levels, cheap first:

1. **Unit.** Call `matched_negatives` twice on the same positives, then on the
   same positives in a shuffled order. It consumes one RNG stream across the
   positives in sequence, so order matters -- the same question the sequence
   clustering turned out to have.
2. **End to end.** Compare a shuffled source between two independently built
   atlases. This is the one that settles it, and it costs nothing when two
   builds already exist on disk.

    python check_shuffle_determinism.py --a data/atlas_rebuild.db --b data/atlas.db
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import sqlite3
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(REPO))


def signature(db: str, src: str):
    con = sqlite3.connect(db)
    rows = con.execute("select sequence from records where source=? order by record_id",
                       (src,)).fetchall()
    con.close()
    h = hashlib.sha256()
    for (s,) in rows:
        h.update(s.encode())
    return len(rows), h.hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default=str(REPO / "data" / "atlas_rebuild.db"))
    ap.add_argument("--b", default=str(REPO / "data" / "atlas.db"))
    a = ap.parse_args()

    from quadcond.negatives import matched_negatives

    con = sqlite3.connect(a.b)
    pos = [r[0] for r in con.execute(
        "select sequence from records where source='g4stab_experimental_tm' "
        "order by record_id limit 400")]
    srcs = [r[0] for r in con.execute(
        "select distinct source from records where source like 'shuffled::%'")]
    con.close()

    x = matched_negatives(pos, kind="G4", seed=7)
    y = matched_negatives(pos, kind="G4", seed=7)
    order = list(pos); random.Random(1).shuffle(order)
    z = matched_negatives(order, kind="G4", seed=7)

    report = {
        "checked_on": date.today().isoformat(),
        "unit": {
            "n_positives": len(pos),
            "same_input_twice_identical": x == y,
            "order_independent": sorted(x) == sorted(z),
        },
        "end_to_end": {},
    }
    for src in sorted(srcs):
        try:
            sa, sb = signature(a.a, src), signature(a.b, src)
        except sqlite3.Error as exc:
            report["end_to_end"][src] = f"not comparable: {exc}"
            continue
        report["end_to_end"][src] = {"a": sa[1], "b": sb[1], "n": sb[0], "identical": sa == sb}

    comparable = [v for v in report["end_to_end"].values() if isinstance(v, dict)]
    all_same = bool(comparable) and all(v["identical"] for v in comparable)
    report["verdict"] = (
        ["Two independent builds produced identical shuffled sources. The generation is "
         "reproducible, so a published atlas makes the negatives reproducible too."]
        if all_same else
        ["Shuffled sources differ between builds. Do not retrain until this is understood: "
         "g4_fold and im_fold would inherit an unreproducible training set."])
    if not report["unit"]["order_independent"]:
        report["verdict"].append(
            "matched_negatives consumes one RNG stream across the positives in order, so the "
            "negatives depend on the order rows arrive in. That is fine while the row order is "
            "itself deterministic -- insertion order from a deterministic build -- but it means "
            "any change to ingest order silently changes the negatives.")

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "shuffle_determinism.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
