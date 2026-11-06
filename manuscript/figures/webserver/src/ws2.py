"""Figure WS2: the AENNA-3D web interface, composed from 3x screenshots of the
running server (fresh clone of the GitHub repository, v0.5.1 release assets)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

RAW = Path("raw"); OUT = Path("out")
NAVY, INK, MUTED, RULE = "#1a4d6e", "#1b2328", "#5b6870", "#c9d3d8"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42})

def img(name, crop=None):
    im = Image.open(RAW / name).convert("RGB")
    return im.crop(crop) if crop else im

A = img("A_overview.png", (0, 0, 5040, 2150))
B = img("C_card_K140.png"); C = img("D_card_Na150.png"); D = img("E_card_imotif.png")
E = img("F_refusal.png")
F = img("F_compare_conditions.png", (0, 470, 2988, 2070))

W = 7.09                       # 180 mm, double column
gap, lab = 0.12, 0.18          # inter-row gap, room for panel letter/title (in)
hA = W * A.height / A.width
colw = (W - 2 * 0.12) / 3
hBD = max(colw * im.height / im.width for im in (B, C, D))
wE = 2.0; wF = W - wE - 0.14
hE = wE * E.height / E.width; hF = wF * F.height / F.width
h3 = max(hE, hF)
H = 0.36 + hA + gap + lab + hBD + gap + lab + h3 + 0.05

fig = plt.figure(figsize=(W, H), dpi=600)
def put(im, x, ytop, w, letter=None, title=None):
    h = w * im.height / im.width
    ax = fig.add_axes([x / W, (H - ytop - h) / H, w / W, h / H])
    ax.imshow(im, interpolation="lanczos"); ax.set_axis_off()
    for s in ("left", "right", "top", "bottom"):
        ax.spines[s].set_visible(False)
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, lw=0.5, ec=RULE))
    if letter:
        fig.text(x / W, (H - ytop + 0.06) / H, letter, fontsize=10, fontweight="bold", color=INK, va="bottom")
    if title:
        fig.text((x + 0.2) / W, (H - ytop + 0.075) / H, title, fontsize=6.6, color=MUTED, va="bottom")
    return ax, h

# --- A: whole interface, with the three working areas named above it ----------
y = 0.36
axA, _ = put(A, 0, y, W)
fig.text(0, (H - 0.02) / H, "A", fontsize=10, fontweight="bold", color=INK, va="top")
for (x0, x1), name in [((12, 308), "1  Sequence and buffer"),
                       ((320, 1316), "2  Schematic 3D archetype / Compare tools"),
                       ((1328, 1668), "3  Evidence cards")]:
    fx0, fx1 = x0 / 1680 * W, x1 / 1680 * W
    yb = H - y + 0.05
    fig.add_artist(plt.Line2D([fx0 / W, fx1 / W], [yb / H] * 2, color=NAVY, lw=0.8))
    for fx in (fx0, fx1):
        fig.add_artist(plt.Line2D([fx / W] * 2, [yb / H, (yb - 0.035) / H], color=NAVY, lw=0.8))
    fig.text(((fx0 + fx1) / 2 + (0.25 if name.startswith("1") else 0)) / W, (yb + 0.035) / H, name,
             ha="center", va="bottom", fontsize=6.8, color=NAVY, fontweight="bold")

# --- B-D: evidence cards --------------------------------------------------------
y += hA + gap + lab
put(B, 0, y, colw, "B", "Telomeric G4 · 140 mM K⁺, 10 mM Na⁺")
put(C, colw + 0.12, y, colw, "C", "Same strand · no K⁺, 150 mM Na⁺")
put(D, 2 * (colw + 0.12), y, colw, "D", "Telomeric C-strand · reference buffer")

# --- E-F: refusal and the applicability-aware condition sweep ----------------------
y += hBD + gap + lab
put(E, 0, y, wE, "E", "C-strand in cytosol-like buffer: declined")
put(F, wE + 0.14, y, wF, "F", "Compare › Conditions: K⁺ swept 0–150 mM, a verdict per head")

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure_WS2_interface.{ext}", dpi=600, facecolor="white")
fig.savefig(OUT / "Figure_WS2_interface.tiff", dpi=600, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print(f"WS2: {W:.2f} x {H:.2f} in")
