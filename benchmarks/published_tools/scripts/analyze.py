import pandas as pd, numpy as np, json, warnings
warnings.filterwarnings("ignore")
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, r2_score, mean_absolute_error, brier_score_loss
from scipy.stats import spearmanr, pearsonr
rng = np.random.default_rng(0)
T = pd.read_csv("tool_scores_all.csv")
TS = T.set_index("sequence")
NB = 500
def boot(fn, y, s, groups=None):
    y = np.asarray(y); s = np.asarray(s, float)
    ok = ~np.isnan(s); y, s = y[ok], s[ok]
    g = np.asarray(groups)[ok] if groups is not None else np.arange(len(y))
    est = fn(y, s)
    ug = np.unique(g); idx = {k: np.where(g == k)[0] for k in ug}
    vals = []
    for _ in range(NB):
        pick = rng.choice(ug, len(ug), replace=True)
        ii = np.concatenate([idx[k] for k in pick])
        try: vals.append(fn(y[ii], s[ii]))
        except Exception: pass
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return est, lo, hi, int(ok.sum())
auroc = lambda y, s: roc_auc_score(y, s)
auprc = lambda y, s: average_precision_score(y, s)
sp = lambda y, s: spearmanr(y, s)[0]
rows = []
def add(task, tool, metric, res, note=""):
    est, lo, hi, n = res
    rows.append(dict(task=task, tool=tool, metric=metric, value=round(est, 4), ci_low=round(lo, 4), ci_high=round(hi, 4), n=n, note=note))

G4_TOOLS = [("G4Hunter (whole-seq mean)", "g4h_mean"), ("G4Hunter (max 25-nt window)", "g4h_max"),
            ("pqsfinder", "pqsfinder"), ("QGRS Mapper G-score", "qgrs"), ("G4Catchall", "g4catchall_default"),
            ("G4Catchall (+G2)", "g4catchall_g2"), ("G4Boost (default motif search)", "g4boost_prob"), ("G4Boost (canonical G3 motifs)", "g4boost_g3_prob"), ("DeepG4", "deepg4"),
            ("G4detector (K, PQ-neg model)", "g4detector_K_pq"), ("G4detector (K, random-neg model)", "g4detector_K_random"),
            ("G4detector (K, dishuffle model)", "g4detector_K_dishuffle"),
            ("G4mismatch (K)", "g4mismatch"), ("Canonical G3+ regex count", "g4_regex"), ("GC content", "gc")]
def tcol(seqs, col):
    return TS.reindex(seqs)[col].values if col in TS.columns else np.full(len(seqs), np.nan)

# ---------------- T1 G4 folding vs matched shuffles (unique sequences)
d = pd.read_csv("oof_g4_fold_full.csv"); d["sequence"] = d.sequence.str.upper()
u = d.groupby("sequence").agg(y=("y", "first"), pred=("pred", "mean"), group=("group", "first")).reset_index()
add("G4 folding vs dinucleotide shuffles", "QuadCond g4_fold (grouped OOF)", "AUROC", boot(auroc, u.y, u.pred, u.group))
add("G4 folding vs dinucleotide shuffles", "QuadCond g4_fold (grouped OOF)", "AUPRC", boot(auprc, u.y, u.pred, u.group))
for name, col in G4_TOOLS:
    s = tcol(u.sequence, col)
    if np.all(np.isnan(s)): continue
    add("G4 folding vs dinucleotide shuffles", name, "AUROC", boot(auroc, u.y, s, u.group))
    add("G4 folding vs dinucleotide shuffles", name, "AUPRC", boot(auprc, u.y, s, u.group))

# ---------------- T2 external G4-seq
e = pd.read_csv("ext_g4seq.csv")
try:
    q = pd.read_csv("scores_qc_ext.csv"); e = e.merge(q, on="sequence", how="left")
except Exception as ex: print(ex)
for neg in ["pq", "random"]:
    sub = e[(e.negtype == "pos") | (e.negtype == neg)]
    task = f"G4-seq K+ (human) positives vs {'PQS-matched unobserved' if neg=='pq' else 'random genomic'} negatives"
    for name, col in [("QuadCond g4_fold (whole window)", "qc_whole"), ("QuadCond g4_fold (motif scan, max)", "qc_scan"),
                      ("QuadCond g4_fold_genomic (whole window)", "qc_genomic")]:
        if col in sub: add(task, name, "AUROC", boot(auroc, sub.y, sub[col]))
    for name, col in G4_TOOLS:
        s = tcol(sub.sequence, col)
        if np.all(np.isnan(s)): continue
        note = "trained on this benchmark (G4detector held-out split unknown)" if "G4detector" in name else ""
        add(task, name, "AUROC", boot(auroc, sub.y, s), note)

# ---------------- T3 topology
d = pd.read_csv("oof_g4_topology_full.csv")
P = d[["p_parallel", "p_antiparallel", "p_hybrid"]].values
ba = lambda y, s: balanced_accuracy_score(y, s)
add("G4 topology (K+, 3-class)", "QuadCond g4_topology (grouped OOF)", "balanced accuracy", boot(ba, d.y, P.argmax(1), d.group))
add("G4 topology (K+, 3-class)", "Majority class baseline", "balanced accuracy", (1/3, 1/3, 1/3, len(d)))
try:
    g = pd.read_csv("scores_g4sp.csv")
    for m in ["LightGBMClassifier", "CatBoostClassifier", "XGBClassifier"]:
        for kind, note in [("pretrained", "released model; trained on this same dataset (in-sample, optimistic)"),
                           ("grouped", "same estimator + hyper-parameters retrained on QuadCond grouped folds")]:
            cols = [f"{m}_{kind}_p_{c}" for c in ["parallel", "antiparallel", "hybrid"]]
            if all(c in g for c in cols):
                add("G4 topology (K+, 3-class)", f"G4ShapePredictor {m.replace('Classifier','')} ({kind})", "balanced accuracy",
                    boot(ba, g.y, g[cols].values.argmax(1), d.group), note)
except Exception as ex: print("g4sp", ex)

# ---------------- T4 G4 Tm
d = pd.read_csv("oof_g4_tm_full.csv"); ds = pd.read_csv("oof_g4_tm_seqonly.csv")
r2 = lambda y, s: r2_score(y, s); mae = lambda y, s: mean_absolute_error(y, s)
TT = "G4 melting temperature (buffer-resolved)"
for nm, x in [("QuadCond g4_tm (grouped OOF)", d.pred), ("QuadCond g4_tm sequence-only (grouped OOF)", ds.pred)]:
    for met, fn in [("R2", r2), ("MAE (°C)", mae), ("Spearman", sp)]:
        add(TT, nm, met, boot(fn, d.y, x, d.group))
try:
    g = pd.read_csv("scores_g4stab.csv")
    for met, fn in [("R2", r2), ("MAE (°C)", mae), ("Spearman", sp)]:
        add(TT, "G4STAB released ensemble", met, boot(fn, d.y, g.ensemble_mean, d.group), "trained on this same dataset (in-sample, optimistic)")
except Exception as ex: print(ex)
try:
    g = pd.read_csv("scores_g4stab_retrained.csv")
    for sch in ["grouped", "random"]:
        for met, fn in [("R2", r2), ("MAE (°C)", mae), ("Spearman", sp)]:
            add(TT, f"G4STAB architecture retrained ({sch} 5-fold)", met, boot(fn, d.y, g[f"g4stab_retrained_{sch}"], d.group))
except Exception as ex: print("g4stab retrain", ex)
try:
    q = pd.read_csv("oof_g4_tm_random.csv")
    for met, fn in [("R2", r2), ("MAE (°C)", mae), ("Spearman", sp)]:
        add(TT, "QuadCond g4_tm (random 5-fold, for comparison with published G4STAB protocol)", met, boot(fn, q.y, q.pred, d.group))
except Exception as ex: print("qc random", ex)
for name, col, sign in [("G4Hunter (whole-seq mean)", "g4h_mean", 1), ("pqsfinder", "pqsfinder", 1), ("QGRS Mapper G-score", "qgrs", 1),
                        ("G4Boost (−predicted MFE)", "g4boost_mfe", -1), ("G4Boost G3 (−predicted MFE)", "g4boost_g3_mfe", -1), ("G4mismatch (K)", "g4mismatch", 1), ("DeepG4", "deepg4", 1)]:
    s = tcol(d.sequence.str.upper(), col)
    if np.all(np.isnan(s)): continue
    add(TT, name, "Spearman", boot(sp, d.y, sign * s, d.group), "buffer-blind score")

# ---------------- T4b within-sequence buffer response (same sequence, different buffers)
WT = "Within-sequence buffer response of G4 Tm"
dd0 = pd.read_csv("oof_g4_tm_full.csv"); dd0["sequence"] = dd0.sequence.str.upper()
def within(pred):
    x = dd0.assign(_p=pred)
    rs, ws = [], []
    for sq, g in x.groupby("sequence"):
        if g.y.nunique() >= 3 and len(g) >= 4:
            dy = g.y - g.y.mean(); dp = g._p - g._p.mean()
            rs.append(spearmanr(g.y, g._p)[0] if g._p.nunique() > 1 else 0.0)
            ws.append((dy, dp))
    rs = np.nan_to_num(np.array(rs))
    alldy = np.concatenate([a for a, b in ws]); alldp = np.concatenate([b for a, b in ws])
    return rs, alldy, alldp
def add_within(name, pred, note=""):
    rs, dy, dp = within(pred)
    vals = [np.median(rng.choice(rs, len(rs))) for _ in range(NB)]
    rows.append(dict(task=WT, tool=name, metric="median per-sequence Spearman", value=round(float(np.median(rs)), 4),
                     ci_low=round(float(np.percentile(vals, 2.5)), 4), ci_high=round(float(np.percentile(vals, 97.5)), 4), n=len(rs), note=note))
    rows.append(dict(task=WT, tool=name, metric="R2 of centred Tm (buffer effect only)", value=round(float(1 - np.sum((dy - dp) ** 2) / np.sum(dy ** 2)), 4),
                     ci_low=np.nan, ci_high=np.nan, n=len(dy), note=note))
add_within("QuadCond g4_tm (grouped OOF)", dd0.pred.values)
add_within("QuadCond g4_tm sequence-only (grouped OOF)", pd.read_csv("oof_g4_tm_seqonly.csv").pred.values, "no buffer input: cannot respond")
try:
    add_within("G4STAB released ensemble", pd.read_csv("scores_g4stab.csv").ensemble_mean.values, "in-sample")
except Exception as ex: print(ex)
try:
    add_within("G4STAB architecture retrained (grouped 5-fold)", pd.read_csv("scores_g4stab_retrained.csv").g4stab_retrained_grouped.values)
except Exception as ex: print(ex)
add_within("Any buffer-blind tool (G4Hunter, pqsfinder, QGRS, G4Catchall, G4Boost, DeepG4, G4detector, G4mismatch)", np.zeros(len(dd0)), "constant across buffers by construction")

# ---------------- T5 single-substitution ΔTm
p = pd.read_csv("pairs_g4_tm.csv")
VT = "Single-substitution effect on G4 Tm (ΔTm, same buffer)"
pg = p.ref  # group by reference
def sign_acc(y, s):
    # directional accuracy on |dTm| >= 2 C; a predicted change of exactly 0 is a coin flip (0.5)
    m = np.abs(y) >= 2
    ys, ss = np.sign(y[m]), np.sign(s[m])
    return float(np.where(ss == 0, 0.5, (ys == ss).astype(float)).mean())
add(VT, "QuadCond g4_tm (Δ of grouped OOF predictions)", "Spearman", boot(sp, p.d_true, p.d_qc, pg))
add(VT, "QuadCond g4_tm (Δ of grouped OOF predictions)", "sign accuracy |ΔTm|≥2°C", boot(sign_acc, p.d_true, p.d_qc, pg))
dd = pd.read_csv("oof_g4_tm_full.csv"); dd["sequence"] = dd.sequence.str.upper()
keycols = ["k", "na", "li_nh4", "mg", "ph", "strand_conc"]
def delta_from(values):
    dd["_v"] = values
    lut = dd.groupby(["sequence"] + keycols, dropna=False)._v.mean()
    def get(s, row): return lut.get(tuple([s] + [row[c] for c in keycols]), np.nan)
    return np.array([get(r.alt, r) - get(r.ref, r) for _, r in p.iterrows()])
try:
    g = pd.read_csv("scores_g4stab.csv")
    x = delta_from(g.ensemble_mean.values)
    add(VT, "G4STAB released ensemble (Δ)", "Spearman", boot(sp, p.d_true, x, pg), "in-sample")
    add(VT, "G4STAB released ensemble (Δ)", "sign accuracy |ΔTm|≥2°C", boot(sign_acc, p.d_true, x, pg), "in-sample")
except Exception as ex: print(ex)
try:
    g = pd.read_csv("scores_g4stab_retrained.csv")
    x = delta_from(g.g4stab_retrained_grouped.values)
    add(VT, "G4STAB architecture retrained, grouped (Δ)", "Spearman", boot(sp, p.d_true, x, pg))
    add(VT, "G4STAB architecture retrained, grouped (Δ)", "sign accuracy |ΔTm|≥2°C", boot(sign_acc, p.d_true, x, pg))
except Exception as ex: print(ex)
for name, col in [("ΔG4Hunter (G4SNVHunter core score)", "g4h_mean"), ("Δpqsfinder", "pqsfinder"), ("ΔQGRS G-score", "qgrs"),
                  ("ΔG4mismatch", "g4mismatch"), ("ΔDeepG4", "deepg4")]:
    if col not in TS.columns: continue
    x = TS.reindex(p.alt)[col].values - TS.reindex(p.ref)[col].values
    add(VT, name, "Spearman", boot(sp, p.d_true, x, pg), "buffer-blind")
    add(VT, name, "sign accuracy |ΔTm|≥2°C", boot(sign_acc, p.d_true, x, pg), "buffer-blind")
try:
    g = pd.read_csv("scores_g4snvhunter.csv")
    x = g.set_index("pair_id").reindex(range(len(p)))["delta"].values
    add(VT, "G4SNVHunter G4VarImpact (ΔG4Hunter window score)", "Spearman", boot(sp, p.d_true, x, pg), "buffer-blind; package run")
    add(VT, "G4SNVHunter G4VarImpact (ΔG4Hunter window score)", "sign accuracy |ΔTm|≥2°C", boot(sign_acc, p.d_true, x, pg), "buffer-blind; package run")
except Exception as ex: print("g4snv", ex)

# ---------------- T6 iM folding
d = pd.read_csv("oof_im_fold_full.csv"); d["sequence"] = d.sequence.str.upper()
u = d.groupby("sequence").agg(y=("y", "first"), pred=("pred", "mean"), group=("group", "first")).reset_index()
IT = "i-motif folding vs dinucleotide shuffles"
add(IT, "QuadCond im_fold (grouped OOF)", "AUROC", boot(auroc, u.y, u.pred, u.group))
for name, col, sign, fill in [("iM-Seeker (folding probability)", "imseeker_prob", 1, 0.0),
                              ("G4Hunter on C-strand (−min window)", "g4h_min", -1, None),
                              ("Canonical C3+ i-motif regex count", "im_regex", 1, None), ("GC content", "gc", 1, None)]:
    s = tcol(u.sequence, col).astype(float)
    if fill is not None: s = np.where(np.isnan(s), fill, s)
    add(IT, name, "AUROC", boot(auroc, u.y, sign * s, u.group),
        "iM-Seeker returns no call when no C3+ motif is found; scored 0" if "iM-Seeker" in name else "")
# ---------------- T7 iM pHT
d = pd.read_csv("oof_im_pht_full.csv"); d["sequence"] = d.sequence.str.upper()
PT = "i-motif transitional pH (pH_T)"
for met, fn in [("R2", r2), ("Spearman", sp)]:
    add(PT, "QuadCond im_pht (grouped OOF)", met, boot(fn, d.y, d.pred, d.group))
s = tcol(d.sequence, "imseeker_strength")
add(PT, "iM-Seeker folding strength", "Spearman", boot(sp, d.y, s, d.group), "trained on these pH_T values (in-sample); min-max scaled so R2 not comparable")
try:
    gi = pd.read_csv("scores_imseeker_grouped.csv")
    okm = gi.imseeker_grouped.notna()
    add(PT, "iM-Seeker architecture retrained (grouped folds)", "Spearman", boot(sp, gi.y[okm], gi.imseeker_grouped[okm], gi.group[okm]), "137/160 sequences with an iM-Seeker motif call")
    add(PT, "iM-Seeker architecture retrained (grouped folds)", "R2", boot(r2, gi.y[okm], gi.imseeker_grouped[okm], gi.group[okm]), "137/160")
    dq = d.set_index(d.index)
    add(PT, "QuadCond im_pht on the same 137 sequences", "Spearman", boot(sp, d.y[okm.values], d.pred[okm.values], d.group[okm.values]))
    add(PT, "QuadCond im_pht on the same 137 sequences", "R2", boot(r2, d.y[okm.values], d.pred[okm.values], d.group[okm.values]))
except Exception as ex: print("ims grouped", ex)
s = tcol(d.sequence, "g4h_min")
add(PT, "G4Hunter on C-strand (−min window)", "Spearman", boot(sp, d.y, -s, d.group), "")

R = pd.DataFrame(rows)
R.to_csv("benchmark_results.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70); pd.set_option("display.max_rows", 300)
print(R.to_string(index=False))
