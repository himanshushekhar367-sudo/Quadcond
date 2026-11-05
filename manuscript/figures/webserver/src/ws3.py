"""Figure WS3: what /scan/conditions returns for one sequence. Predicted Tm of
the 21-nt telomeric repeat as K+ or Na+ is swept 0-150 mM, the other cation
held at zero. Data are the server's JSON responses, saved in raw/."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K_COL, NA_COL = "#0A8F9C", "#B35F02"          # validated categorical pair (light)
INK, MUTED, GRID = "#1b2328", "#5b6870", "#e3e8eb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42,
                     "axes.edgecolor": "#9aa7b0", "axes.linewidth": 0.6,
                     "xtick.color": MUTED, "ytick.color": MUTED, "xtick.major.width": 0.6, "ytick.major.width": 0.6})

def series(name):
    d = json.load(open(Path("raw") / name))
    s = d["series"]["g4_tm"]
    prov = d["run_record"]["provenance"]["model"]
    assert prov["matches_manifest"] is True, "model does not match the release manifest"
    assert all(p["applicability"]["state"] == "in_domain" for p in s["points"])
    xs = [p["value_of_axis"] for p in s["points"]]
    ys = [p["value"] for p in s["points"]]
    lo = [p["interval"][0] for p in s["points"]]; hi = [p["interval"][1] for p in s["points"]]
    return xs, ys, lo, hi, d, prov

kx, ky, klo, khi, kd, prov = series("sweep_k.json")
nx, ny, nlo, nhi, nd, _ = series("sweep_na.json")

fig, ax = plt.subplots(figsize=(3.4, 2.75), dpi=600)
fig.subplots_adjust(left=0.15, right=0.80, bottom=0.16, top=0.95)
ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)

for xs, ys, lo, hi, col, lab in [(kx, ky, klo, khi, K_COL, "K⁺ (no Na⁺)"), (nx, ny, nlo, nhi, NA_COL, "Na⁺ (no K⁺)")]:
    ax.fill_between(xs, lo, hi, color=col, alpha=0.10, lw=0)
    ax.plot(xs, ys, color=col, lw=1.3, label=lab, zorder=3)
    ax.plot(xs, ys, "o", ms=3.2, color=col, mec="white", mew=0.6, zorder=4)

# direct labels at the line ends (ink text, the mark carries the colour)
ax.annotate(f"K⁺  {ky[-1]:.1f} °C", (kx[-1], ky[-1]), xytext=(6, 0), textcoords="offset points",
            va="center", fontsize=6.8, color=INK, annotation_clip=False)
ax.annotate(f"Na⁺  {ny[-1]:.1f} °C", (nx[-1], ny[-1]), xytext=(6, 0), textcoords="offset points",
            va="center", fontsize=6.8, color=INK, annotation_clip=False)

# the 100 mM comparison
i = kx.index(100.0)
ax.annotate("", xy=(100, ny[i]), xytext=(100, ky[i]),
            arrowprops={"arrowstyle": "<->", "color": INK, "lw": 0.7, "shrinkA": 2, "shrinkB": 2})
ax.text(97.5, (ky[i] + ny[i]) / 2, f"ΔTm {ky[i] - ny[i]:.1f} °C\nat 100 mM", ha="right", va="center", fontsize=6.3, color=INK, linespacing=1.15)
ax.annotate(f"no added cation\n{ky[0]:.1f} °C", (0, ky[0]), xytext=(9, -2), textcoords="offset points",
            va="center", fontsize=6.3, color=MUTED)

ax.set_xlim(-4, 154); ax.set_ylim(20, 86)
ax.set_xticks([0, 25, 50, 75, 100, 125, 150])
ax.set_xlabel("Monovalent cation (mM)", color=INK)
ax.set_ylabel("Predicted Tm, g4_tm (°C)", color=INK)
leg = ax.legend(loc="upper left", frameon=False, fontsize=6.6, handlelength=1.6)

out = Path("out")
for ext in ("png", "pdf"):
    fig.savefig(out / f"Figure_WS3_condition_sweep.{ext}", dpi=600, facecolor="white")
fig.savefig(out / "Figure_WS3_condition_sweep.tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("model sha", prov["model_file_sha256"][:16], "| K+", ky[0], "->", ky[-1], "| Na+", ny[0], "->", ny[-1], "| 100 mM", ky[i], ny[i])
print("base", kd["base_condition"], "| imputed", kd["condition_imputed_fields"])
