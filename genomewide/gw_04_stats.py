#!/usr/bin/env python3
"""Step 4: do structure-changing SNVs carry larger AlphaGenome effects?

Classes (per SNV, per motif kind):
  motif_lost        substitution destroys the canonical motif
  destabilising     delta <= -T_strong   (G4: degC; iM: pH units)
  moderate          -T_strong < delta <= -T_weak
  neutral           |delta| < T_weak
  stabilising       delta >= T_weak
  control           SNV in the matched out-of-motif control window

Tests (all written to $GW_OUT/stats/):
  1. |AVI| distribution per class; Mann-Whitney + Cliff's delta vs neutral and vs control
  2. enrichment among the top 1% |AVI| (threshold from control SNVs): Fisher OR
  3. adjusted logistic model (class + substitution type + CpG + motif GC + chromosome),
     standard errors clustered by motif
  4. WITHIN-MOTIF paired test: per motif, mean |AVI| of motif-lost SNVs vs
     motif-retained SNVs (Wilcoxon signed-rank). This removes locus-level
     regulatory activity, which is the main confounder of tests 1-3.
  5. dose-response: Spearman(-delta, |AVI|) among retained SNVs
  6. leave-one-chromosome-out: OR of test 2 per held-out chromosome
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
from gw_common import CHROMS, out_dir, read_chrom  # noqa: E402

THRESH = {"G4": (2.0, 5.0), "iM": (0.1, 0.3)}
KEY = ["chrom", "pos", "ref", "alt"]


def cliffs_delta(a, b, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.choice(a, min(n, len(a)), replace=False)
    b = rng.choice(b, min(n, len(b)), replace=False)
    u = stats.mannwhitneyu(a, b).statistic
    return 2 * u / (len(a) * len(b)) - 1


def classify(df, kind):
    weak, strong = THRESH[kind]
    d = pd.to_numeric(df.delta, errors="coerce")
    cls = np.select(
        [df.state.eq("motif_lost"), d <= -strong, d <= -weak, d.abs() < weak, d >= weak],
        ["motif_lost", "destabilising", "moderate", "neutral", "stabilising"], default="other")
    return pd.Series(cls, index=df.index)


def subst_type(ref, alt):
    pair = ref + alt
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if ref in "GA":        # collapse to pyrimidine reference
        pair = comp[ref] + comp[alt]
    return pair[0] + ">" + pair[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--score", help="AVI column to use (default: first column with 'avi' in its name)")
    a = ap.parse_args()
    base = out_dir()
    sd = base / "stats"
    sd.mkdir(exist_ok=True)
    frames = []
    for chrom in a.chroms:
        sp, vp = base / "structural" / f"{chrom}.tsv.gz", base / "avi" / f"{chrom}.tsv.gz"
        if not (sp.exists() and vp.exists()):
            print(f"{chrom}: missing structural or avi output, skipped")
            continue
        s = pd.read_csv(sp, sep="\t", dtype={"delta": str})
        v = pd.read_csv(vp, sep="\t")
        motifs = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t", usecols=["id", "gc", "ctrl_gc"])
        col = a.score or next((c for c in v.columns[6:] if "avi" in c.lower()), v.columns[6])
        v["abs_avi"] = pd.to_numeric(v[col], errors="coerce").abs()
        vm = v[v.set == "motif"].drop_duplicates(KEY)[KEY + ["abs_avi"]]
        s = s.merge(vm, on=KEY, how="inner")
        s["cls"] = [None] * len(s)
        for kind in ("G4", "iM"):
            k = s.kind == kind
            s.loc[k, "cls"] = classify(s[k], kind)
        s = s.merge(motifs[["id", "gc"]], left_on="motif_id", right_on="id", how="left").drop(columns="id")
        c = v[v.set == "control"].merge(motifs[["id", "ctrl_gc"]], left_on="motif_id", right_on="id", how="left")
        kinds = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t", usecols=["id", "kind"]).set_index("id").kind
        c = pd.DataFrame({"chrom": chrom, "pos": c.pos, "ref": c.ref, "alt": c.alt, "motif_id": c.motif_id,
                          "kind": c.motif_id.map(kinds), "cls": "control", "abs_avi": c.abs_avi,
                          "gc": c.ctrl_gc})
        both = pd.concat([s[["chrom", "pos", "ref", "alt", "motif_id", "kind", "cls", "abs_avi", "gc", "delta"]], c])
        seq = read_chrom(a.fasta, chrom)
        # CpG context: the reference base is the C or the G of a CpG (pos is 1-based).
        both["cpg"] = [int(seq[p - 2:p] == "CG" or seq[p - 1:p + 1] == "CG") for p in both.pos]
        both["subst"] = [subst_type(r, x) for r, x in zip(both.ref, both.alt)]
        frames.append(both.dropna(subset=["abs_avi"]))
        print(f"{chrom}: {len(frames[-1])} SNVs with AVI", flush=True)
    if not frames:
        raise SystemExit(
            "nothing to analyse: no chromosome has both a structural and an AVI table.\n"
            f"  structural: {base / 'structural'}\n"
            f"  avi:        {base / 'avi'}\n"
            "Run gw_03_avi.py (Tabix bundle) or gw_03_avi_live.py (sampled, live API) first.")
    D = pd.concat(frames, ignore_index=True)
    if _has_parquet():
        D.to_parquet(sd / "snv_table.parquet")
    else:
        D.to_csv(sd / "snv_table.tsv.gz", sep="\t", index=False)
    report = {"score_column": col, "n_snv": int(len(D)), "kinds": {}}
    order = ["motif_lost", "destabilising", "moderate", "neutral", "stabilising", "control"]
    for kind in ("G4", "iM"):
        K = D[D.kind == kind]
        ctrl = K[K.cls == "control"].abs_avi.values
        if len(ctrl) < 100:
            continue
        top = float(np.quantile(ctrl, 0.99))
        K = K.assign(top1=(K.abs_avi >= top).astype(int))
        rep = {"top1pct_threshold_from_controls": top, "classes": {}}
        neutral = K[K.cls == "neutral"].abs_avi.values
        for cl in order:
            x = K[K.cls == cl]
            if x.empty:
                continue
            r = {"n": int(len(x)), "median_abs_avi": float(x.abs_avi.median()),
                 "frac_top1pct": float(x.top1.mean())}
            if cl != "control":
                t = [[x.top1.sum(), len(x) - x.top1.sum()],
                     [int((ctrl >= top).sum()), int((ctrl < top).sum())]]
                r["OR_vs_control"], r["p_OR_vs_control"] = stats.fisher_exact(t)
                r["cliffs_delta_vs_control"] = cliffs_delta(x.abs_avi.values, ctrl)
                r["p_mwu_vs_control"] = float(stats.mannwhitneyu(x.abs_avi, ctrl).pvalue)
                if cl != "neutral" and len(neutral):
                    r["cliffs_delta_vs_neutral"] = cliffs_delta(x.abs_avi.values, neutral)
                    r["p_mwu_vs_neutral"] = float(stats.mannwhitneyu(x.abs_avi, neutral).pvalue)
            rep["classes"][cl] = r
        # 3. adjusted model
        try:
            import statsmodels.formula.api as smf
            M = K[K.cls.isin(order)].copy()
            M["cls"] = pd.Categorical(M.cls, categories=["control"] + [o for o in order if o != "control"])
            fit = smf.logit("top1 ~ C(cls) + C(subst) + cpg + gc + C(chrom)", data=M).fit(
                disp=0, cov_type="cluster", cov_kwds={"groups": pd.factorize(M.motif_id)[0]})
            rep["adjusted_logit"] = {
                n.replace("C(cls)[T.", "").rstrip("]"): {"OR": float(np.exp(b)),
                 "ci": [float(np.exp(lo)), float(np.exp(hi))], "p": float(p)}
                for n, b, (lo, hi), p in zip(fit.params.index, fit.params,
                                             fit.conf_int().values, fit.pvalues)
                if n.startswith("C(cls)")}
        except Exception as exc:  # noqa: BLE001
            rep["adjusted_logit"] = f"not fitted: {exc}"
        # 4. within-motif paired test
        W = K[K.cls != "control"].copy()
        W["lost"] = W.cls.eq("motif_lost")
        g = W.groupby(["motif_id", "lost"]).abs_avi.mean().unstack()
        g = g.dropna()
        if len(g) > 20:
            w = stats.wilcoxon(g[True], g[False])
            rep["within_motif_lost_vs_retained"] = {
                "n_motifs": int(len(g)), "median_diff": float((g[True] - g[False]).median()),
                "frac_motifs_lost_higher": float((g[True] > g[False]).mean()), "p_wilcoxon": float(w.pvalue)}
        # 5. dose-response
        R = K[K.cls.isin(["destabilising", "moderate", "neutral", "stabilising"])]
        dd = pd.to_numeric(R.delta, errors="coerce")
        rho = stats.spearmanr(-dd, R.abs_avi, nan_policy="omit")
        rep["dose_response_spearman_minus_delta_vs_abs_avi"] = {"rho": float(rho.statistic), "p": float(rho.pvalue)}
        # 6. leave-one-chromosome-out
        loco = {}
        for ch in sorted(K.chrom.unique()):
            x = K[(K.chrom == ch)]
            lost = x[x.cls == "motif_lost"]
            cc = x[x.cls == "control"]
            if len(lost) > 50 and len(cc) > 50:
                loco[ch] = float(stats.fisher_exact([[lost.top1.sum(), len(lost) - lost.top1.sum()],
                                                     [cc.top1.sum(), len(cc) - cc.top1.sum()]])[0])
        rep["per_chromosome_OR_motif_lost_vs_control"] = loco
        report["kinds"][kind] = rep
        _figure(K, order, kind, sd)
    (sd / "summary.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=float)[:4000])


def _has_parquet():
    try:
        import pyarrow  # noqa: F401
        return True
    except ImportError:
        return False


def _figure(K, order, kind, sd):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    for cl, col in zip(order, ["#2a78d6", "#4a3aa7", "#1baf7a", "#9d9c97", "#eda100", "#52514e"]):
        x = np.sort(K[K.cls == cl].abs_avi.values)
        if len(x):
            ax.plot(x, 1 - np.arange(len(x)) / len(x), label=f"{cl} (n={len(x):,})", color=col, lw=1.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("|AVI score|")
    ax.set_ylabel("fraction of SNVs >= x")
    ax.set_title(f"{kind}: AlphaGenome effect by QuadCond structural class", loc="left")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(sd / f"{kind}_avi_by_class.png", dpi=200)


if __name__ == "__main__":
    main()
