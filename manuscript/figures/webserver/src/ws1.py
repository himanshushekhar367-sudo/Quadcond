"""Figure WS1: how the QuadCond web server is put together. Same palette, font
and box style as manuscript Figure 1 (build_manuscript.py)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

NAVY, GREEN, OCHRE, PURPLE = "#184e72", "#226d59", "#9a6020", "#68527b"
FACE, ARROW, INK, MUTED = "#f7f9fb", "#607080", "#1b2328", "#4a5965"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "pdf.fonttype": 42})

fig = plt.figure(figsize=(7.09, 6.1), dpi=600)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 10); ax.set_ylim(0, 8.6); ax.set_axis_off()

def box(x, y, w, h, title=None, body=None, color=NAVY, face=FACE, lw=0.9, tsize=8, bsize=7, dash=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                lw=lw, ec=color, fc=face, ls=(0, (3, 2)) if dash else "-"))
    if title and body:
        ax.text(x + w / 2, y + h - 0.13, title, ha="center", va="top", color=color, weight="bold", fontsize=tsize)
        ax.text(x + w / 2, y + (h - 0.3) / 2, body, ha="center", va="center", color=INK, fontsize=bsize, linespacing=1.35)
    elif title:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", color=color, weight="bold", fontsize=tsize)

def arrow(x1, y1, x2, y2, label=None, both=False):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops={"arrowstyle": "<|-|>" if both else "-|>", "color": ARROW, "lw": 1.1})
    if label:
        ax.text((x1 + x2) / 2 + 0.12, (y1 + y2) / 2, label, ha="left", va="center", fontsize=6.8, color=MUTED, style="italic")

def band(y, h, label):
    ax.text(0.05, y + h - 0.02, label, ha="left", va="top", fontsize=7.2, color=MUTED, weight="bold")

# ---- browser ----------------------------------------------------------------
band(7.0, 1.55, "BROWSER  ·  AENNA-3D workbench (React, TypeScript, Three.js)")
box(0.05, 7.0, 3.05, 1.2, "Structure view", "Schematic 3D archetype per strand\nat published helical parameters")
box(3.475, 7.0, 3.05, 1.2, "Evidence cards", "Value, interval, folded fraction,\nclaim limit and source per head")
box(6.9, 7.0, 3.05, 1.2, "Compare tools", "Variant join · Mutations · Conditions\nBatch · Model guide")

arrow(5.0, 6.98, 5.0, 6.42, "same origin · HTTP + JSON · no credentials", both=True)

# ---- the one process -------------------------------------------------------
ax.add_patch(FancyBboxPatch((0.05, 1.72), 9.9, 4.66, boxstyle="round,pad=0.02,rounding_size=0.1",
                            lw=1.3, ec=NAVY, fc="white"))
ax.text(0.22, 6.24, "quadcond serve", ha="left", va="top", fontsize=9, weight="bold", color=NAVY, family="DejaVu Sans Mono")
ax.text(2.05, 6.235, "one process, one port: the workbench at / and the API beside it", ha="left", va="top", fontsize=7.4, color=MUTED)

# endpoints row
eps = [("Health", "/health\n/ready"), ("Describe", "/info\n/schema/prediction"),
       ("Score", "/predict\n/batch"), ("Scan", "/scan/mutations\n/scan/conditions\n/scan/genome"),
       ("Join", "/scan/variant")]
x = 0.25; w = (9.55 - 4 * 0.12) / 5
for title, routes in eps:
    ax.add_patch(FancyBboxPatch((x, 4.72), w, 1.18, boxstyle="round,pad=0.02,rounding_size=0.06", lw=0.7, ec="#9fb3c4", fc=FACE))
    ax.text(x + w / 2, 5.80, title, ha="center", va="top", fontsize=7.4, weight="bold", color=NAVY)
    ax.text(x + w / 2, 5.18, routes, ha="center", va="center", fontsize=6.6, family="DejaVu Sans Mono", color=INK, linespacing=1.3)
    x += w + 0.12

ax.text(5.0, 4.52, "Request checks:  DNA only (U is refused as RNA)  ·  /predict ≤ 64 sequences, each ≤ 5,000 nt  ·  /scan/genome ≤ 20 kb",
        ha="center", va="center", fontsize=6.7, color=MUTED)
arrow(5.0, 4.38, 5.0, 4.12)

# heads
ax.text(5.0, 4.05, "12 prediction heads", ha="center", va="top", fontsize=7.8, weight="bold", color=INK)
hw = 3.0
box(0.35, 2.95, hw, 0.82, "Biophysically anchored (7)", "Tm, topology, transitional pH,\nfolding vs matched shuffles", GREEN, bsize=6.6, tsize=7.4)
box(3.5, 2.95, hw, 0.82, "Genomic proxies (3)", "CUT&Tag peak association;\nnot a folding probability", OCHRE, bsize=6.6, tsize=7.4)
box(6.65, 2.95, hw, 0.82, "Derived or predicted (2)", "Labels from another model\nor from the sequence itself", PURPLE, bsize=6.6, tsize=7.4)
arrow(5.0, 2.93, 5.0, 2.70)

# applicability + provenance
box(0.35, 1.92, 4.55, 0.74, "Applicability gate, per head",
    "Supported estimate  |  extrapolation hidden  |  refused, no number", NAVY, bsize=6.6, tsize=7.4)
box(5.1, 1.92, 4.55, 0.74, "Provenance block on every response",
    "model SHA-256 · matches_manifest (true / false / unknown)\nbuild versus runtime environment", NAVY, bsize=6.4, tsize=7.4)

# ---- what feeds the process --------------------------------------------------
band(0.05, 1.45, "")
arrow(2.55, 1.30, 2.55, 1.70)
arrow(7.45, 1.30, 7.45, 1.70)
box(0.05, 0.05, 4.95, 1.22, "Release assets, verified before use",
    "model · atlas_core · G4-seq scanner from GitHub release v0.5.1\nSHA-256 checked against assets_manifest.json when the\nimage is built; a mismatch stops the build", NAVY, bsize=6.5, tsize=7.4)
box(5.2, 0.05, 4.75, 1.22, "Optional: AlphaGenome variant scores",
    "AVI Tabix bundle · table export · live API\nexact join on chromosome, position, ref, alt;\nnever by proximity", NAVY, bsize=6.5, tsize=7.4, dash=True)

ax.text(9.93, 6.24, "no accounts · no sessions · no stored requests", ha="right", va="top", fontsize=6.6, color=MUTED, style="italic")

out = Path("out")
for ext in ("png", "pdf"):
    fig.savefig(out / f"Figure_WS1_architecture.{ext}", dpi=600, facecolor="white")
fig.savefig(out / "Figure_WS1_architecture.tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("WS1 written")
