import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np
from scipy.stats import spearmanr
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
QC, OTHER, INS = "#2a78d6", "#9d9c97", "#d6d5d0"
R = pd.read_csv("benchmark_results.csv")
def panel(ax, task, metric, title, xlabel, ref=None, drop=()):
    s = R[(R.task == task) & (R.metric == metric)]
    s = s[~s.tool.isin(drop)].sort_values("value")
    y = np.arange(len(s))
    cols = [QC if t.startswith("QuadCond") else (INS if "in-sample" in str(n) or "trained on this" in str(n) else OTHER)
            for t, n in zip(s.tool, s.note.fillna(""))]
    ax.barh(y, s.value, color=cols, height=0.62, edgecolor="#fcfcfb", linewidth=1)
    err = np.array([s.value - s.ci_low, s.ci_high - s.value]).clip(0)
    ax.errorbar(s.value, y, xerr=err, fmt="none", ecolor="#52514e", elinewidth=0.8, capsize=2)
    ax.set_yticks(y, [t if len(t) < 48 else t[:46] + "…" for t in s.tool], fontsize=7.5)
    if ref is not None: ax.axvline(ref, color="#52514e", lw=0.8, ls="--")
    ax.set_title(title, loc="left", fontsize=10, color="#0b0b0b"); ax.set_xlabel(xlabel)
    ax.grid(axis="y", visible=False)
T1 = "G4 folding vs dinucleotide shuffles"
T2p = "G4-seq K+ (human) positives vs PQS-matched unobserved negatives"
T2r = "G4-seq K+ (human) positives vs random genomic negatives"
fig, axs = plt.subplots(1, 3, figsize=(16, 6.2), constrained_layout=True)
panel(axs[0], T1, "AUROC", "A  G4 folding vs matched shuffles (atlas, grouped CV)", "AUROC", 0.5)
panel(axs[1], T2r, "AUROC", "B  External: G4-seq K+ vs random genome", "AUROC", 0.5)
panel(axs[2], T2p, "AUROC", "C  External: G4-seq K+ vs unobserved PQS", "AUROC", 0.5)
fig.text(0.01, -0.02, "Blue = QuadCond; grey = published tool; light grey = tool trained on the evaluation data (optimistic). Bars: 95% bootstrap CI.", fontsize=8, color="#52514e")
fig.savefig("fig1_g4_folding.png", dpi=180, bbox_inches="tight")

fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
panel(axs[0], "G4 topology (K+, 3-class)", "balanced accuracy", "A  G4 topology (grouped CV)", "balanced accuracy", 1/3)
panel(axs[1], "G4 melting temperature (buffer-resolved)", "Spearman", "B  G4 Tm ranking", "Spearman ρ", 0)
panel(axs[2], "i-motif folding vs dinucleotide shuffles", "AUROC", "C  i-motif folding vs matched shuffles", "AUROC", 0.5)
fig.savefig("fig2_topology_tm_im.png", dpi=180, bbox_inches="tight")

fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
VT = "Single-substitution effect on G4 Tm (ΔTm, same buffer)"
panel(axs[0], VT, "Spearman", "A  Single-substitution ΔTm (1,020 measured pairs)", "Spearman ρ (predicted Δ vs measured ΔTm)", 0)
panel(axs[1], VT, "sign accuracy |ΔTm|≥2°C", "B  Direction of effect, |ΔTm| ≥ 2 °C", "directional accuracy (ties = 0.5)", 0.5)
p = pd.read_csv("pairs_g4_tm.csv")
ax = axs[2]
ax.scatter(p.d_true, p.d_qc, s=9, color=QC, alpha=0.45, edgecolors="none")
lim = [-35, 50]; ax.plot(lim, lim, color="#52514e", lw=0.8, ls="--")
ax.set_xlim(lim); ax.set_ylim(-25, 25); ax.axhline(0, color="#8a8984", lw=0.6); ax.axvline(0, color="#8a8984", lw=0.6)
ax.set_xlabel("measured ΔTm (°C)"); ax.set_ylabel("QuadCond predicted ΔTm (°C)")
ax.set_title(f"C  QuadCond ΔTm, ρ = {spearmanr(p.d_true, p.d_qc)[0]:.2f}", loc="left", fontsize=10)
fig.savefig("fig3_variant_dtm.png", dpi=180, bbox_inches="tight")

WT = "Within-sequence buffer response of G4 Tm"
fig, axs = plt.subplots(1, 2, figsize=(12, 3.8), constrained_layout=True)
panel(axs[0], WT, "median per-sequence Spearman", "A  Same sequence, different buffers", "median per-sequence Spearman ρ", 0)
panel(axs[1], "i-motif transitional pH (pH_T)", "Spearman", "B  i-motif pH_T", "Spearman ρ", 0)
fig.savefig("fig4_condition_pht.png", dpi=180, bbox_inches="tight")
print("ok")
