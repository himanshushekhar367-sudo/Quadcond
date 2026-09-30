#!/usr/bin/env python3
"""Check a rebuilt atlas against the composition of the one the models were trained on.

The full atlas does not exist any more (docs/REBUILD_THE_ATLAS.md). What survives
is its per-source composition, read out of the frozen model sidecar into
original_atlas_composition.csv. This compares a rebuild against that, source by
source, because "the build finished" and "the build produced the right atlas" are
different statements and the builder skips a missing source with a warning.

    python docs/atlas/compare_rebuild.py --db data/atlas.db

Writes original_vs_rebuild.csv beside this script and prints the table.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(REPO / "data" / "atlas.db"))
    ap.add_argument("--spec", default=str(ROOT / "original_atlas_composition.csv"))
    a = ap.parse_args()
    if not Path(a.db).exists():
        raise SystemExit(f"no atlas at {a.db}; build one first")

    orig = {r["source"]: r for r in csv.DictReader(open(a.spec, encoding="utf-8"))}
    con = sqlite3.connect(a.db)
    new = {src: {"n": n, "n_tm": tm or 0, "n_pht": pht or 0} for src, n, tm, pht in con.execute(
        "select source, count(*), sum(tm is not null), sum(ph_t is not null) "
        "from records group by source")}

    rows, tot_o, tot_n = [], 0, 0
    for src in sorted(orig, key=lambda s: -int(orig[s]["n"])):
        o = orig[src]
        on, otm, opht = int(o["n"]), int(o.get("n_tm") or 0), int(o.get("n_pht") or 0)
        tot_o += on
        got = new.get(src)
        nn, ntm, npht = (got["n"], got["n_tm"], got["n_pht"]) if got else (0, 0, 0)
        tot_n += nn
        if not got:
            status = "MISSING"
        elif (nn, ntm, npht) == (on, otm, opht):
            status = "exact"
        else:
            # "short by -1834" is how the first version reported an excess.
            def gap(orig, got, unit):
                d = orig - got
                return None if d == 0 else f"{abs(d):,} {unit} {'short' if d > 0 else 'over'}"
            parts = [g for g in (gap(on, nn, "rows"), gap(otm, ntm, "Tm"),
                                 gap(opht, npht, "pH_T")) if g]
            status = ", ".join(parts)
        rows.append({"source": src, "orig_n": on, "rebuild_n": nn, "orig_tm": otm,
                     "rebuild_tm": ntm, "orig_pht": opht, "rebuild_pht": npht,
                     "status": status})
    for src in sorted(set(new) - set(orig)):
        tot_n += new[src]["n"]
        rows.append({"source": src, "orig_n": 0, "rebuild_n": new[src]["n"], "orig_tm": 0,
                     "rebuild_tm": new[src]["n_tm"], "orig_pht": 0,
                     "rebuild_pht": new[src]["n_pht"], "status": "NOT IN ORIGINAL"})

    with (ROOT / "original_vs_rebuild.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print(f"{'source':42} {'orig':>8} {'rebuild':>8}  status")
    for r in rows:
        print(f"{r['source'][:42]:42} {r['orig_n']:>8,} {r['rebuild_n']:>8,}  {r['status']}")
    exact = sum(r["status"] == "exact" for r in rows)
    print(f"\n{tot_n:,} of {tot_o:,} rows ({tot_n / tot_o:.1%}); "
          f"{exact} of {len(orig)} sources exact")
    print(f"written: {ROOT / 'original_vs_rebuild.csv'}")


if __name__ == "__main__":
    main()
