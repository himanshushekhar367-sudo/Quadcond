#!/usr/bin/env python3
"""Compare the grouping recorded at training time against the grouping replayed today.

The shipped model sidecars record summary statistics about the cross-validation
grouping -- how many groups a head's rows fell into -- but not the assignment
itself. `quadcond/models/foldmanifest.py` reconstructs an assignment by replaying
`build_folds` on rows read today, and honestly labels the result
`assignment_provenance: reconstructed` rather than `historical`.

This script asks the question that flag invites: does the replay actually land on
the same grouping the training run used? If it does for every head, the historical
partition is recoverable and the only remaining problem is publishing it. If it
does not, the partition is gone and the manuscript has to say which of its numbers
are records rather than reproducible claims.

It reads two files and computes; it fits nothing and needs no heavy dependencies.

    python compare_sidecar_vs_manifest.py

Writes results/sidecar_vs_manifest.csv and .json.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
DEFAULT_SIDECAR = REPO / "manuscript" / "inputs" / "quadcond_model.json"
DEFAULT_MANIFEST = REPO / "benchmarks" / "published_tools" / "results" / "fold_assignments_meta.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sidecar", default=str(DEFAULT_SIDECAR))
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    a = ap.parse_args()

    for p in (a.sidecar, a.manifest):
        if not Path(p).exists():
            raise SystemExit(f"not found: {p}")

    side = json.loads(Path(a.sidecar).read_text(encoding="utf-8"))["heads"]
    man = json.loads(Path(a.manifest).read_text(encoding="utf-8"))

    rows, mismatched = [], []
    for head, m in man.get("heads", {}).items():
        meta = side.get(head, {}).get("training_meta", {})
        grouping = meta.get("grouping", {})
        rec = {
            "head": head,
            "group_by": m.get("group_by"),
            "rows_sidecar": meta.get("n_rows"),
            "rows_replay": m.get("n_rows"),
            "groups_sidecar": grouping.get("n_groups"),
            "groups_replay": m.get("n_groups"),
            "cluster_threshold": m.get("cluster_threshold"),
            "seed": m.get("seed"),
        }
        rec["rows_match"] = rec["rows_sidecar"] == rec["rows_replay"]
        rec["groups_match"] = rec["groups_sidecar"] == rec["groups_replay"]
        rec["reproduced"] = bool(rec["rows_match"] and rec["groups_match"])
        rows.append(rec)
        if not rec["reproduced"]:
            mismatched.append(head)

    by_grouping = {}
    for r in rows:
        g = by_grouping.setdefault(r["group_by"] or "unknown", {"heads": 0, "reproduced": 0})
        g["heads"] += 1
        g["reproduced"] += int(r["reproduced"])

    report = {
        "compared_on": date.today().isoformat(),
        "sidecar": str(a.sidecar),
        "manifest": str(a.manifest),
        "manifest_assignment_provenance": man.get("assignment_provenance"),
        "manifest_rows_fingerprint": man.get("rows_fingerprint"),
        "heads_compared": len(rows),
        "heads_reproduced": sum(r["reproduced"] for r in rows),
        "heads_not_reproduced": mismatched,
        "by_grouping_variable": by_grouping,
        "rows": rows,
        "reading": reading(rows, by_grouping, man),
    }

    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    cols = list(rows[0]) if rows else []
    with (out / "sidecar_vs_manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    (out / "sidecar_vs_manifest.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    width = max(len(r["head"]) for r in rows) if rows else 8
    print(f"{'head':{width}} {'group_by':>10} {'rows':>14} {'groups':>14}  reproduced")
    for r in rows:
        print(f"{r['head']:{width}} {str(r['group_by']):>10} "
              f"{str(r['rows_sidecar']) + '/' + str(r['rows_replay']):>14} "
              f"{str(r['groups_sidecar']) + '/' + str(r['groups_replay']):>14}  "
              f"{'yes' if r['reproduced'] else 'NO'}")
    print()
    for line in report["reading"]:
        print("- " + line)
    print(f"\nwritten: results/sidecar_vs_manifest.csv and .json")


def reading(rows, by_grouping, man):
    out = []
    bad = [r for r in rows if not r["reproduced"]]
    if not bad:
        out.append("Every head's replayed grouping matches the training record. The historical "
                   "partition is recoverable; publish it and the fold question closes.")
        return out
    if all(r["rows_match"] for r in rows):
        out.append("Row counts match for every head, so the rows being read are the right rows. "
                   "The disagreement is in the grouping, not in the data.")
    seq = by_grouping.get("sequence")
    if seq and seq["reproduced"] == 0:
        others = {k: v for k, v in by_grouping.items() if k != "sequence"}
        if others and all(v["reproduced"] == v["heads"] for v in others.values()):
            out.append("Every sequence-grouped head fails to reproduce and every head grouped on "
                       "anything else reproduces exactly. The non-determinism is isolated to the "
                       "sequence clustering path.")
    out.append("_greedy_clusters orders sequences with np.argsort over negative lengths, and "
               "NumPy's default sort is not stable. Sequences sharing a length are permuted "
               "differently between runs, and a greedy pass then selects different centroids. "
               "check_determinism.py tests that mechanism directly.")
    out.append(f"The manifest reports assignment_provenance={man.get('assignment_provenance')!r}, "
               "which is the correct label: it must not be used to clear a head of leakage, and "
               "metrics computed under it are not the shipped artifact's historical out-of-fold "
               "figures.")
    return out


if __name__ == "__main__":
    main()
