#!/usr/bin/env python3
"""Assay-sensitivity controls for the genome-wide AVI analysis.

gw_04_stats.py returns a null for the motif tests. A null is only worth
reporting if the readout can detect anything at all at this scale, so this
script asks the question a reviewer will ask first: does |AVI| respond to
sequence features already known to carry regulatory weight, measured in
exactly the same windows?

Everything here is computed on the FLANK rows alone -- the within-locus
background -- so no control depends on a motif call being right.

  1. CpG.  A SNV whose reference base sits in a CpG dinucleotide should
     carry more weight than one that does not.  If it does not, the readout
     is flat and the motif null says nothing.
  2. Substitution type.  The six collapsed substitutions should not be
     interchangeable.
  3. Variance decomposition.  How much of the |AVI| rank do CpG and
     substitution type explain, against what motif class explains?

It then restates the motif result as an equivalence statement: the largest
effect the data are compatible with, in percentile points.

Reads gw_out/stats/snv_table.parquet (or .tsv.gz) written by gw_04_stats.py.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gw_common import out_dir  # noqa: E402

SUBSTS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def load(sd: Path) -> pd.DataFrame:
    pq, tsv = sd / "snv_table.parquet", sd / "snv_table.tsv.gz"
    if pq.exists():
        return pd.read_parquet(pq)
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t")
    raise SystemExit(f"no SNV table in {sd}; run gw_04_stats.py first")


def cliffs_delta(a, b, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.choice(a, min(n, len(a)), replace=False)
    b = rng.choice(b, min(n, len(b)), replace=False)
    u = stats.mannwhitneyu(a, b).statistic
    return float(2 * u / (len(a) * len(b)) - 1)


def clustered(formula, data, groups):
    import statsmodels.formula.api as smf
    return smf.ols(formula, data=data).fit(cov_type="cluster", cov_kwds={"groups": groups})


def contrast(model, term):
    ci = model.conf_int().loc[term]
    return {"delta_pctile": float(model.params[term]),
            "ci": [float(ci[0]), float(ci[1])],
            "p": float(model.pvalues[term])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="run directory (default: $GW_OUT or ./gw_out)")
    a = ap.parse_args()
    sd = (Path(a.out) if a.out else out_dir()) / "stats"

    D = load(sd)
    # Percentile rank of |AVI| within chromosome: the same outcome gw_04 models,
    # so coefficients here are on the same scale as the motif coefficients.
    D["pct"] = D.groupby("chrom").abs_avi.rank(pct=True)
    D["cpg"] = D.cpg.astype(int)

    report = {"n_snv": int(len(D)), "kinds": {}}
    for kind in ("G4", "iM"):
        K = D[D.kind == kind]
        F = K[K.cls == "flank"].dropna(subset=["pct", "subst"])
        if len(F) < 1000:
            continue
        rep = {"n_flank": int(len(F))}

        # 1. CpG, inside the flank background only.
        yes = F[F.cpg == 1].abs_avi.values
        no = F[F.cpg == 0].abs_avi.values
        m = clustered("pct ~ cpg", F, F.motif_id)
        rep["cpg_in_flank"] = {
            "n_cpg": int(len(yes)), "n_non_cpg": int(len(no)),
            "median_abs_avi_cpg": float(np.median(yes)),
            "median_abs_avi_non_cpg": float(np.median(no)),
            "cliffs_delta": cliffs_delta(yes, no),
            "p_mwu": float(stats.mannwhitneyu(yes, no).pvalue),
            "adjusted": contrast(m, "cpg"),
        }

        # 2. Substitution type, inside the flank background only.
        present = [s for s in SUBSTS if (F.subst == s).sum() >= 200]
        sub = {}
        if len(present) > 1:
            ref = present[0]
            f2 = F[F.subst.isin(present)].copy()
            f2["subst"] = pd.Categorical(f2.subst, categories=present)
            ms = clustered("pct ~ C(subst)", f2, f2.motif_id)
            for s in present[1:]:
                sub[s] = contrast(ms, f"C(subst)[T.{s}]")
            rep["subst_in_flank"] = {"reference": ref, "terms": sub,
                                     "p_joint": float(ms.f_pvalue)}

        # 3. What explains the rank: sequence composition, or motif class?
        S = K.dropna(subset=["pct", "subst"])
        S = S[S.subst.isin(present)]
        if len(S) > 300000:      # R2 is a ratio; a large sample estimates it fine
            S = S.sample(300000, random_state=0)
        r2 = {}
        for name, formula in (("cpg_plus_subst", "pct ~ cpg + C(subst)"),
                              ("motif_class", "pct ~ C(cls)"),
                              ("both", "pct ~ cpg + C(subst) + C(cls)")):
            import statsmodels.formula.api as smf
            r2[name] = float(smf.ols(formula, data=S).fit().rsquared)
        rep["variance_explained_r2"] = r2

        report["kinds"][kind] = rep

    # 4. Restate the motif result as an equivalence bound.
    summ = sd / "summary.json"
    if summ.exists():
        s = json.loads(summ.read_text())
        eq = {}
        for kind, v in s.get("kinds", {}).items():
            t = v.get("adjusted_rank_model", {}).get("terms", {}).get("motif_lost")
            if not t:
                continue
            hi = max(abs(t["ci"][0]), abs(t["ci"][1]))
            eq[kind] = {
                "motif_lost_vs_flank_delta_pctile": t["delta_pctile"],
                "ci": t["ci"],
                "largest_effect_compatible_pctile_points": round(hi * 100, 3),
                "reading": (
                    "adjusted for CpG and substitution type and clustered by motif, "
                    "destroying the motif shifts the |AVI| percentile by at most "
                    f"{hi * 100:.2f} percentile points in either direction"),
            }
        report["equivalence"] = eq

    path = sd / "sensitivity.json"
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(f"\nwritten: {path}")


if __name__ == "__main__":
    main()
