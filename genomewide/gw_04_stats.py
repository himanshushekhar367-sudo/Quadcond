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
        try:
            v = pd.read_csv(vp, sep="\t")
        except Exception as exc:                              # noqa: BLE001
            print(f"{chrom}: AVI table unreadable ({exc}); skipped")
            continue
        # The structural table is every SNV in every motif on the chromosome --
        # millions of rows -- while the AVI table is the sampled or bundled
        # subset actually joined. Reading the whole structural file to throw
        # nearly all of it away is what made this step run out of patience (and
        # memory) on chr1; it is streamed and filtered instead.
        wanted = set(v.pos.astype(int))
        parts, kept, seen = [], 0, 0
        try:
            for chunk in pd.read_csv(sp, sep="\t", dtype={"delta": str}, chunksize=2_000_000):
                seen += len(chunk)
                part = chunk[chunk.pos.isin(wanted)]
                if len(part):
                    parts.append(part)
                    kept += len(part)
        except Exception as exc:                              # noqa: BLE001
            print(f"{chrom}: structural table unreadable after {seen:,} rows ({exc}); skipped")
            continue
        if not parts:
            print(f"{chrom}: no structural row matches an AVI position; skipped")
            continue
        s = pd.concat(parts, ignore_index=True)
        print(f"{chrom}: {kept:,} of {seen:,} structural rows at AVI positions", flush=True)
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
        f = v[v.set == "flank"].merge(motifs[["id", "gc"]], left_on="motif_id", right_on="id", how="left")
        kinds0 = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t",
                             usecols=["id", "kind"]).set_index("id").kind
        f = pd.DataFrame({"chrom": chrom, "pos": f.pos, "ref": f.ref, "alt": f.alt,
                          "motif_id": f.motif_id, "kind": f.motif_id.map(kinds0),
                          "cls": "flank", "abs_avi": f.abs_avi, "gc": f.gc})
        c = v[v.set == "control"].merge(motifs[["id", "ctrl_gc"]], left_on="motif_id", right_on="id", how="left")
        kinds = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t", usecols=["id", "kind"]).set_index("id").kind
        c = pd.DataFrame({"chrom": chrom, "pos": c.pos, "ref": c.ref, "alt": c.alt, "motif_id": c.motif_id,
                          "kind": c.motif_id.map(kinds), "cls": "control", "abs_avi": c.abs_avi,
                          "gc": c.ctrl_gc})
        both = pd.concat([s[["chrom", "pos", "ref", "alt", "motif_id", "kind", "cls", "abs_avi", "gc", "delta"]], f, c])
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
    order = ["motif_lost", "destabilising", "moderate", "neutral", "stabilising", "flank", "control"]
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
        # 3. adjusted models.
        #
        # The rank model is the primary one. Top-1% membership is rare among motif
        # SNVs (tens of events), and a logistic fit on that separates: the first
        # run of this returned an odds ratio of 1e-8 with a convergence warning,
        # which is a fitting artefact and not a finding. A rank-transformed
        # outcome uses the whole distribution and cannot separate.
        try:
            import statsmodels.formula.api as smf
            M = K[K.cls.isin(order)].copy()
            present = [o for o in order if (M.cls == o).any()]
            base_cls = "flank" if (M.cls == "flank").sum() > 100 else "control"
            if base_cls not in present:
                base_cls = present[-1]
            # Only classes that HAVE rows. Declaring a category with no rows gives
            # patsy an all-zero dummy column, and that -- not the data -- is what
            # made the design matrix rank-deficient.
            M["cls"] = pd.Categorical(M.cls, categories=[base_cls] + [o for o in present if o != base_cls])
            M["avi_rank"] = M.groupby("chrom").abs_avi.rank(pct=True)
            fit_r = smf.ols("avi_rank ~ C(cls) + C(subst) + cpg + gc + C(chrom)", data=M).fit(
                cov_type="cluster", cov_kwds={"groups": pd.factorize(M.motif_id)[0]})
            rep["adjusted_rank_model"] = {
                "reference_class": base_cls,
                "note": "outcome is |AVI| percentile rank within chromosome; a coefficient is a "
                        "shift in mean percentile against the reference class, clustered by motif",
                "terms": {n.replace("C(cls)[T.", "").rstrip("]"):
                          {"delta_pctile": float(b), "ci": [float(lo), float(hi)], "p": float(p)}
                          for n, b, (lo, hi), p in zip(fit_r.params.index, fit_r.params,
                                                       fit_r.conf_int().values, fit_r.pvalues)
                          if n.startswith("C(cls)")}}
        except Exception as exc:                              # noqa: BLE001
            rep["adjusted_rank_model"] = f"not fitted: {exc}"
        try:
            events = int(K.top1.sum())
            if events < 50:
                raise ValueError(f"only {events} top-1% events in this kind; a logistic fit on "
                                 "that many separates rather than estimates")
            import statsmodels.formula.api as smf
            M = K[K.cls.isin(order)].copy()
            present = [o for o in order if (M.cls == o).any()]
            # Same reference class as the rank model above. The two were reported
            # against different baselines in the first version, which made the
            # odds ratios and the percentile shifts look like they disagreed.
            logit_base = "flank" if (M.cls == "flank").sum() > 100 else "control"
            if logit_base not in present:
                logit_base = present[-1]
            M["cls"] = pd.Categorical(M.cls, categories=[logit_base] + [o for o in present if o != logit_base])
            fit = smf.logit("top1 ~ C(cls) + C(subst) + cpg + gc + C(chrom)", data=M).fit(
                disp=0, cov_type="cluster", cov_kwds={"groups": pd.factorize(M.motif_id)[0]})
            rep["adjusted_logit"] = {
                "reference_class": logit_base,
                "note": "outcome is membership of the top 1% of |AVI| as defined by the distant "
                        "control distribution; rare among motif SNVs, so read the rank model first",
                "terms": {
                    n.replace("C(cls)[T.", "").rstrip("]"): {"OR": float(np.exp(b)),
                     "ci": [float(np.exp(lo)), float(np.exp(hi))], "p": float(p)}
                    for n, b, (lo, hi), p in zip(fit.params.index, fit.params,
                                                 fit.conf_int().values, fit.pvalues)
                    if n.startswith("C(cls)")}}
        except Exception as exc:  # noqa: BLE001
            rep["adjusted_logit"] = f"not fitted: {exc}"
        # What the classes differ on besides structure. The distant control is
        # matched on GC per motif, but a motif whose GC cannot be matched within
        # 2-20 kb gets no control at all, so the control set is drawn from the
        # easier, lower-GC loci -- a selection effect that no adjustment in the
        # model can see. Reported so it is not invisible.
        rep["covariate_balance"] = {
            cl: {"n": int((K.cls == cl).sum()),
                 "mean_gc": round(float(K.gc[K.cls == cl].mean()), 4),
                 "frac_cpg": round(float(K.cpg[K.cls == cl].mean()), 4)}
            for cl in order if (K.cls == cl).any()}
        rep["top1_note"] = ("`frac_top1pct` for the control class is 1% by construction: the "
                            "threshold IS its 99th percentile. Read the other classes against "
                            "that, and prefer the rank model and the within-locus tests.")
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
        # 4b. motif-destroying SNVs vs the flanks of the same window. The distant
        # control window is matched on GC but not on what the locus does, and its
        # |AVI| tail runs well above the motif's; the flank holds promoter,
        # chromatin and conservation fixed by construction.
        FL = K[K.cls == "flank"]
        if len(FL) > 100:
            lost_m = K[K.cls == "motif_lost"].groupby("motif_id").abs_avi.mean()
            fl_m = FL.groupby("motif_id").abs_avi.mean()
            pair = pd.concat([lost_m.rename("lost"), fl_m.rename("flank")], axis=1).dropna()
            if len(pair) > 20:
                w = stats.wilcoxon(pair.lost, pair.flank)
                rep["within_locus_lost_vs_flank"] = {
                    "n_motifs": int(len(pair)),
                    "median_diff": float((pair.lost - pair.flank).median()),
                    "frac_motifs_lost_higher": float((pair.lost > pair.flank).mean()),
                    "p_wilcoxon": float(w.pvalue)}
            ret_m = K[K.cls.isin(["destabilising", "moderate", "neutral", "stabilising"])] \
                .groupby("motif_id").abs_avi.mean()
            pair2 = pd.concat([ret_m.rename("motif"), fl_m.rename("flank")], axis=1).dropna()
            if len(pair2) > 20:
                w2 = stats.wilcoxon(pair2.motif, pair2.flank)
                rep["within_locus_motif_vs_flank"] = {
                    "n_motifs": int(len(pair2)),
                    "median_diff": float((pair2.motif - pair2.flank).median()),
                    "frac_motifs_higher": float((pair2.motif > pair2.flank).mean()),
                    "p_wilcoxon": float(w2.pvalue)}

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
