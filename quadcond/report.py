"""Self-contained HTML report: atlas composition, calibration, ablation, case studies.

No external JS. Fonts come from Google Fonts (the one host an Artifact may load
from); everything else -- charts included -- is inline SVG generated here, so the
file opens from disk, survives being emailed, and can be published as-is.
"""
from __future__ import annotations

import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from . import claims
from .models.base import MultiTaskModel

# Validated categorical slots (dataviz reference palette, all-pairs clean at 3).
SERIES = ["var(--s1)", "var(--s2)", "var(--s3)"]
TIER_SLOT = {"experimental": "var(--s1)", "derived": "var(--s2)", "predicted": "var(--s3)"}
KEY_METRIC = {"binary": "auroc", "multiclass": "balanced_accuracy", "regression": "r2"}
METRIC_LABEL = {"auroc": "AUROC", "balanced_accuracy": "balanced accuracy", "r2": "R²"}


def _e(x) -> str:
    return html.escape(str(x))


def _num(v, nd=4) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return "—" if v != v else f"{v:.{nd}f}"
    if isinstance(v, int):
        return f"{v:,}"
    return _e(v)


# --------------------------------------------------------------------------- #
# Charts
# --------------------------------------------------------------------------- #
def _stacked_composition(items, total, title) -> str:
    """One 100% bar: how the atlas is composed, with counts labelled outside."""
    if not total:
        return ""
    w, barh = 620, 26
    x = 0.0
    segs, labels = [], []
    for i, (name, n, colour) in enumerate(items):
        frac = n / total
        seg_w = frac * w
        gap = 2 if i else 0
        segs.append(
            f'<rect x="{x + gap:.2f}" y="0" width="{max(seg_w - gap, 3):.2f}" height="{barh}" '
            f'fill="{colour}" rx="3"><title>{_e(name)}: {n:,} records ({frac:.1%})</title></rect>'
        )
        labels.append((name, n, frac, colour))
        x += seg_w
    legend = "".join(
        f'<span class="key"><i style="background:{c}"></i>{_e(nm)}'
        f'<b>{n:,}</b><em>{f:.1%}</em></span>'
        for nm, n, f, c in labels
    )
    return (
        f'<figure class="fig"><figcaption>{_e(title)}</figcaption>'
        f'<svg viewBox="0 0 {w} {barh}" preserveAspectRatio="none" class="compbar" '
        f'role="img" aria-label="{_e(title)}">{"".join(segs)}</svg>'
        f'<div class="keys">{legend}</div></figure>'
    )


def _log_bars(rows, title, note="") -> str:
    """Horizontal bars on a log10 scale -- counts here span three orders of magnitude."""
    if not rows:
        return ""
    labels = [r[0] for r in rows]
    vals = [max(float(r[1]), 1.0) for r in rows]
    colours = [r[2] for r in rows]
    padl, padr, rowh, padt = 244, 92, 26, 24
    w = 700
    iw = w - padl - padr
    vmax = math.log10(max(vals)) or 1.0
    h = padt + len(rows) * rowh + 26
    parts = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{_e(title)}">']
    for d in range(0, int(math.ceil(vmax)) + 1):
        gx = padl + iw * (d / vmax)
        parts.append(f'<line x1="{gx:.1f}" y1="{padt - 6}" x2="{gx:.1f}" y2="{padt + len(rows)*rowh}" '
                     f'stroke="var(--grid)" />')
        parts.append(f'<text x="{gx:.1f}" y="{padt + len(rows)*rowh + 16}" text-anchor="middle" '
                     f'class="tick">{10**d:,}</text>')
    for i, (lab, v, col) in enumerate(zip(labels, vals, colours)):
        y = padt + i * rowh
        bw = max(2.0, iw * (math.log10(v) / vmax))
        short = lab if len(lab) <= 34 else lab[:32] + "\u2026"
        parts.append(f'<text x="{padl-10}" y="{y+14}" text-anchor="end" class="blabel">{_e(short)}'
                     f'<title>{_e(lab)}</title></text>')
        parts.append(f'<rect x="{padl}" y="{y+3}" width="{bw:.1f}" height="{rowh-9}" fill="{col}" rx="3">'
                     f'<title>{_e(lab)}: {int(v):,}</title></rect>')
        parts.append(f'<text x="{padl+bw+7:.1f}" y="{y+14}" class="bval">{int(v):,}</text>')
    parts.append("</svg>")
    n = f'<p class="figNote">{_e(note)}</p>' if note else ""
    return (f'<figure class="fig"><figcaption>{_e(title)}</figcaption>'
            f'<div class="scroll">{"".join(parts)}</div>{n}</figure>')


def _paired_bars(rows, title, name_a, name_b, note="") -> str:
    """Two series per row, direct-labelled -- the condition ablation."""
    if not rows:
        return ""
    padl, padr, padt, rowh = 200, 66, 44, 46
    w = 640
    iw = w - padl - padr
    vmax = max(max(r[1], r[2]) for r in rows) or 1.0
    h = padt + len(rows) * rowh + 12
    parts = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{_e(title)}">']
    parts.append(f'<rect x="{padl}" y="16" width="10" height="10" rx="2" fill="{SERIES[0]}"/>'
                 f'<text x="{padl+16}" y="25" class="blabel">{_e(name_a)}</text>'
                 f'<rect x="{padl+168}" y="16" width="10" height="10" rx="2" fill="{SERIES[1]}"/>'
                 f'<text x="{padl+184}" y="25" class="blabel">{_e(name_b)}</text>')
    for i, (lab, a, b) in enumerate(rows):
        y = padt + i * rowh
        parts.append(f'<text x="{padl-10}" y="{y+22}" text-anchor="end" class="blabel">{_e(lab)}</text>')
        for j, (v, col) in enumerate(((a, SERIES[0]), (b, SERIES[1]))):
            bw = max(2.0, iw * (v / vmax))
            yy = y + 2 + j * 17
            parts.append(f'<rect x="{padl}" y="{yy}" width="{bw:.1f}" height="13" fill="{col}" rx="3">'
                         f'<title>{_e(lab)} — {name_a if j == 0 else name_b}: {v:.4f}</title></rect>')
            parts.append(f'<text x="{padl+bw+7:.1f}" y="{yy+11}" class="bval">{v:.3f}</text>')
    parts.append("</svg>")
    n = f'<p class="figNote">{_e(note)}</p>' if note else ""
    return (f'<figure class="fig"><figcaption>{_e(title)}</figcaption>'
            f'<div class="scroll">{"".join(parts)}</div>{n}</figure>')


def _reliability(rel, title) -> str:
    conf, acc, cnt = rel.get("confidence", []), rel.get("accuracy", []), rel.get("count", [])
    if not conf:
        return '<p class="figNote">No reliability data for this head.</p>'
    w = h = 268
    pad = 42
    iw = ih = w - 2 * pad

    def X(v):
        return pad + v * iw

    def Y(v):
        return pad + (1 - v) * ih

    p = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{_e(title)}">']
    p.append(f'<rect x="{pad}" y="{pad}" width="{iw}" height="{ih}" fill="var(--plot)" '
             f'stroke="var(--grid)" rx="2"/>')
    for g in (0.0, 0.25, 0.5, 0.75, 1.0):
        p.append(f'<line x1="{X(g):.1f}" y1="{pad}" x2="{X(g):.1f}" y2="{pad+ih}" stroke="var(--grid)"/>')
        p.append(f'<line x1="{pad}" y1="{Y(g):.1f}" x2="{pad+iw}" y2="{Y(g):.1f}" stroke="var(--grid)"/>')
        p.append(f'<text x="{X(g):.1f}" y="{pad+ih+15}" text-anchor="middle" class="tick">{g:g}</text>')
        p.append(f'<text x="{pad-7}" y="{Y(g)+3:.1f}" text-anchor="end" class="tick">{g:g}</text>')
    p.append(f'<line x1="{X(0)}" y1="{Y(0)}" x2="{X(1)}" y2="{Y(1)}" stroke="var(--muted)" '
             f'stroke-dasharray="3 4"/>')
    p.append(f'<polyline points="{" ".join(f"{X(c):.1f},{Y(a):.1f}" for c, a in zip(conf, acc))}" '
             f'fill="none" stroke="{SERIES[0]}" stroke-width="2"/>')
    mx = max(cnt) if cnt else 1
    for c, a, n in zip(conf, acc, cnt):
        r = 3.5 + 4.0 * (n / mx) ** 0.5
        p.append(f'<circle cx="{X(c):.1f}" cy="{Y(a):.1f}" r="{r:.1f}" fill="{SERIES[0]}" '
                 f'stroke="var(--plot)" stroke-width="2">'
                 f'<title>predicted {c:.3f} · observed {a:.3f} · n={n}</title></circle>')
    p.append(f'<text x="{pad+iw/2:.0f}" y="{h-6}" text-anchor="middle" class="axis">predicted probability</text>')
    p.append(f'<text x="11" y="{pad+ih/2:.0f}" text-anchor="middle" class="axis" '
             f'transform="rotate(-90 11 {pad+ih/2:.0f})">observed frequency</text>')
    p.append("</svg>")
    return (f'<figure class="fig"><figcaption>{_e(title)}</figcaption>'
            f'<div class="scroll">{"".join(p)}</div>'
            f'<p class="figNote">Dashed diagonal is perfect calibration; marker area is bin count.</p>'
            f'</figure>')


def _line_chart(xs, series, xlabel, ylabel, title, note="") -> str:
    """Small multiples-free line chart for a condition sweep."""
    if not xs or not series:
        return ""
    w, h = 620, 280
    padl, padr, padt, padb = 62, 108, 22, 44
    iw, ih = w - padl - padr, h - padt - padb
    xmin, xmax = min(xs), max(xs)
    allv = [v for _, vals in series for v in vals]
    ymin, ymax = min(allv), max(allv)
    if ymax - ymin < 1e-9:
        ymin, ymax = ymin - 0.5, ymax + 0.5
    padv = (ymax - ymin) * 0.08
    prob_axis = min(allv) >= 0.0 and max(allv) <= 1.0
    ymin, ymax = ymin - padv, ymax + padv
    if prob_axis:  # probabilities never leave [0, 1]; a negative tick is a chart bug
        ymin, ymax = max(ymin, 0.0), min(ymax, 1.0)

    def X(v):
        return padl + (v - xmin) / (xmax - xmin) * iw if xmax > xmin else padl

    def Y(v):
        return padt + (1 - (v - ymin) / (ymax - ymin)) * ih

    p = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{_e(title)}">']
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4
        p.append(f'<line x1="{padl}" y1="{Y(yv):.1f}" x2="{padl+iw}" y2="{Y(yv):.1f}" stroke="var(--grid)"/>')
        p.append(f'<text x="{padl-8}" y="{Y(yv)+3:.1f}" text-anchor="end" class="tick">{yv:.3g}</text>')
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4
        p.append(f'<text x="{X(xv):.1f}" y="{padt+ih+18}" text-anchor="middle" class="tick">{xv:.3g}</text>')
    for si, (name, vals) in enumerate(series):
        col = SERIES[si % len(SERIES)]
        pts = " ".join(f"{X(x):.1f},{Y(v):.1f}" for x, v in zip(xs, vals))
        p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2" '
                 f'stroke-linejoin="round"/>')
        for x, v in zip(xs, vals):
            p.append(f'<circle cx="{X(x):.1f}" cy="{Y(v):.1f}" r="3" fill="{col}" '
                     f'stroke="var(--plot)" stroke-width="1.5">'
                     f'<title>{name}: {xlabel} {x:g} → {v:.3f}</title></circle>')
        p.append(f'<text x="{padl+iw+8}" y="{Y(vals[-1])+4:.1f}" class="blabel" fill="var(--ink)">'
                 f'{_e(name)}</text>')
    p.append(f'<text x="{padl+iw/2:.0f}" y="{h-6}" text-anchor="middle" class="axis">{_e(xlabel)}</text>')
    p.append(f'<text x="12" y="{padt+ih/2:.0f}" text-anchor="middle" class="axis" '
             f'transform="rotate(-90 12 {padt+ih/2:.0f})">{_e(ylabel)}</text>')
    p.append("</svg>")
    n = f'<p class="figNote">{_e(note)}</p>' if note else ""
    return (f'<figure class="fig"><figcaption>{_e(title)}</figcaption>'
            f'<div class="scroll">{"".join(p)}</div>{n}</figure>')


# --------------------------------------------------------------------------- #
CSS = """
:root{
  color-scheme: light;
  --bg:#f5f7f8; --surface:#ffffff; --plot:#fbfcfc; --raised:#eef2f4;
  --ink:#12171b; --ink2:#39424a; --muted:#68727b; --line:#dde3e7; --grid:#e7ecef;
  --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a;
  --good:#1a6b47; --goodbg:#e6f3ed; --warn:#9a4a1e; --warnbg:#fbeee6;
  --shadow:0 1px 2px rgba(18,23,27,.05), 0 8px 24px -18px rgba(18,23,27,.25);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    color-scheme: dark;
    --bg:#0f1316; --surface:#171c20; --plot:#141a1e; --raised:#1e252a;
    --ink:#e7ecef; --ink2:#c0c8ce; --muted:#8b959d; --line:#283036; --grid:#232b31;
    --s1:#3987e5; --s2:#d95926; --s3:#199e70;
    --good:#7cc7a4; --goodbg:#152a22; --warn:#e0a179; --warnbg:#2c1d14;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -20px rgba(0,0,0,.8);
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --bg:#0f1316; --surface:#171c20; --plot:#141a1e; --raised:#1e252a;
  --ink:#e7ecef; --ink2:#c0c8ce; --muted:#8b959d; --line:#283036; --grid:#232b31;
  --s1:#3987e5; --s2:#d95926; --s3:#199e70;
  --good:#7cc7a4; --goodbg:#152a22; --warn:#e0a179; --warnbg:#2c1d14;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -20px rgba(0,0,0,.8);
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font-family:"IBM Plex Sans",ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  font-size:15px;line-height:1.6;-webkit-font-smoothing:antialiased}
.shell{display:grid;grid-template-columns:216px minmax(0,1fr);gap:40px;
  max-width:1180px;margin:0 auto;padding:36px 26px 96px}
@media (max-width:900px){.shell{grid-template-columns:1fr;gap:20px;padding:24px 18px 64px}
  nav.rail{position:static!important;max-height:none!important}}
nav.rail{position:sticky;top:28px;align-self:start;max-height:calc(100vh - 56px);overflow:auto}
nav.rail .brand{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:17px;
  letter-spacing:-.01em;margin-bottom:2px}
nav.rail .tag{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:10.5px;
  letter-spacing:.09em;text-transform:uppercase;color:var(--muted);margin-bottom:18px}
nav.rail a{display:block;color:var(--ink2);text-decoration:none;font-size:13.5px;
  padding:4px 0 4px 11px;border-left:2px solid var(--line)}
nav.rail a:hover{color:var(--s1);border-left-color:var(--s1)}
nav.rail a.sub{padding-left:22px;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px}
nav.rail .grp{font-size:10.5px;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);
  margin:18px 0 6px}
h1{font-family:"IBM Plex Serif",Georgia,serif;font-size:31px;line-height:1.2;margin:0 0 6px;
  letter-spacing:-.018em;text-wrap:balance}
h2{font-family:"IBM Plex Serif",Georgia,serif;font-size:21px;margin:0 0 4px;letter-spacing:-.012em}
h3{font-size:15px;margin:0}
.lede{color:var(--ink2);max-width:66ch;margin:0 0 4px}
.stamp{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:11.5px;color:var(--muted);
  letter-spacing:.02em}
section{margin-top:44px;scroll-margin-top:24px}
section > .head{border-bottom:1px solid var(--line);padding-bottom:9px;margin-bottom:18px}
section > .head p{margin:5px 0 0;color:var(--muted);font-size:13.5px;max-width:70ch}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:1px;
  background:var(--line);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:22px 0 0}
.tile{background:var(--surface);padding:14px 16px}
.tile .k{font-size:10.5px;letter-spacing:.085em;text-transform:uppercase;color:var(--muted)}
.tile .v{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:25px;line-height:1.25;
  font-variant-numeric:tabular-nums;letter-spacing:-.02em;margin-top:3px}
.tile .s{font-size:12px;color:var(--muted)}
.fig{margin:18px 0 0}
.fig figcaption{font-size:11px;letter-spacing:.085em;text-transform:uppercase;color:var(--muted);
  margin-bottom:9px}
.figNote{font-size:12.5px;color:var(--muted);margin:7px 0 0;max-width:70ch}
.scroll{overflow-x:auto}
svg{max-width:100%;height:auto;display:block}
.compbar{width:100%;height:26px}
.keys{display:flex;flex-wrap:wrap;gap:14px;margin-top:10px}
.key{display:inline-flex;align-items:baseline;gap:6px;font-size:12.5px;color:var(--ink2)}
.key i{width:9px;height:9px;border-radius:2px;display:inline-block;align-self:center}
.key b{font-family:"IBM Plex Mono",ui-monospace,monospace;font-variant-numeric:tabular-nums;
  color:var(--ink)}
.key em{font-style:normal;color:var(--muted)}
.tick,.axis{font-size:10.5px;fill:var(--muted);font-family:"IBM Plex Mono",ui-monospace,monospace}
.blabel{font-size:11.5px;fill:var(--ink2);font-family:"IBM Plex Sans",sans-serif}
.bval{font-size:11px;fill:var(--muted);font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-variant-numeric:tabular-nums}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:16px 0 0}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--grid);vertical-align:top}
thead th{font-size:10.5px;letter-spacing:.075em;text-transform:uppercase;color:var(--muted);
  border-bottom:1px solid var(--line);font-weight:600;white-space:nowrap}
td.n,th.n{text-align:right;font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-variant-numeric:tabular-nums}
code,.mono{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px}
code{background:var(--raised);border-radius:3px;padding:1px 5px}
.seq{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;letter-spacing:.03em;
  word-break:break-all}
.pill{display:inline-flex;align-items:center;gap:5px;font-size:11px;padding:2px 8px;
  border-radius:999px;border:1px solid var(--line);color:var(--muted);white-space:nowrap}
.pill i{width:7px;height:7px;border-radius:50%;display:inline-block}
.pill.good{color:var(--good);background:var(--goodbg);border-color:transparent}
.pill.warn{color:var(--warn);background:var(--warnbg);border-color:transparent}
.headcard{background:var(--surface);border:1px solid var(--line);border-radius:9px;
  box-shadow:var(--shadow);margin:16px 0 0;overflow:hidden}
.headcard.ungrounded{border-left:3px solid var(--s2)}
.headcard > .hc{display:flex;justify-content:space-between;align-items:center;gap:12px;
  flex-wrap:wrap;padding:14px 18px 0}
.headcard .claim{padding:8px 18px 0;color:var(--ink2);font-size:13.5px;max-width:74ch}
.hcbody{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,290px);gap:22px;padding:14px 18px 18px}
@media (max-width:820px){.hcbody{grid-template-columns:1fr}}
.kv{width:100%;font-size:13px}
.kv th{color:var(--muted);font-weight:400;width:52%;text-transform:none;letter-spacing:0;font-size:13px;
  border-bottom:1px solid var(--grid)}
details{border-top:1px solid var(--grid);padding:10px 18px 14px}
summary{cursor:pointer;color:var(--s1);font-size:13px}
summary:focus-visible,a:focus-visible{outline:2px solid var(--s1);outline-offset:2px}
.note{border-left:3px solid var(--s2);background:var(--surface);padding:12px 16px;
  border-radius:0 8px 8px 0;margin:20px 0 0;font-size:14px;box-shadow:var(--shadow)}
.note b{color:var(--warn)}
.src{border-bottom:1px solid var(--grid);padding:11px 0;display:grid;
  grid-template-columns:110px minmax(0,1fr);gap:14px}
@media (max-width:600px){.src{grid-template-columns:1fr;gap:4px}}
.src .meta{font-size:12px;color:var(--muted)}
pre{background:var(--raised);border:1px solid var(--line);border-radius:7px;padding:12px;
  overflow-x:auto;font-size:11.5px;line-height:1.5;max-height:460px}
a{color:var(--s1)}
.delta{font-family:"IBM Plex Mono",ui-monospace,monospace;font-variant-numeric:tabular-nums}
.delta.up{color:var(--good)} .delta.flat{color:var(--muted)}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
"""


def build_report(model_path, db_path, out_path, *, ablation: str | None = None,
                 case_studies: str | None = "artifacts/case_studies.json") -> Path:
    model = MultiTaskModel.load(model_path)
    seq_only = None
    if ablation and Path(ablation).exists():
        try:
            seq_only = MultiTaskModel.load(ablation)
        except Exception:
            seq_only = None
    cases = None
    if case_studies and Path(case_studies).exists():
        try:
            cases = json.loads(Path(case_studies).read_text())
        except Exception:
            cases = None

    from .atlas import Atlas

    atlas = Atlas(db_path)
    summary = atlas.summary()
    sources = [dict(r) for r in atlas.conn.execute(
        "SELECT * FROM sources ORDER BY evidence_tier, source")]
    total = int(summary["n"].sum())
    atlas.close()

    stamp = datetime.now(timezone.utc).strftime("%d %b %Y, %H:%M UTC")

    bench_html = ""
    bench_files = sorted(Path("artifacts").glob("external_benchmark*.json")) \
        if Path("artifacts").exists() else []
    benches = []
    for bf in bench_files:
        try:
            benches.append(json.loads(bf.read_text()))
        except Exception:
            pass
    if benches:
        rows_b, chart_rows = [], []
        for b in benches:
            t = b.get("topology") or {}
            f = b.get("folding") or {}
            held_out = b.get("held_out", True)
            reportable = b.get("reportable_metrics") or {}
            if held_out:
                topo_val = t.get("balanced_accuracy")
                fold_val = f.get("auroc")
                basis = "<span class='pill good'>held out</span>"
                extra = ""
            else:
                # The model trained on these records, so the in-sample numbers are
                # memorisation. Report its own grouped out-of-fold metrics instead
                # and say what the in-sample figures were, so the gap is visible.
                topo_val = reportable.get("g4_topology")
                fold_val = reportable.get("g4_fold")
                basis = "<span class='pill warn'>in-sample &rarr; out-of-fold shown</span>"
                extra = (f"<div class='meta'>in-sample was "
                         f"{_num(t.get('balanced_accuracy'), 3)} / "
                         f"{_num(f.get('auroc'), 3)} &mdash; the gap is what a naive "
                         f"self-evaluation would have added</div>")
            rows_b.append(
                f"<tr><td>{_e(b.get('model', '?'))}{extra}</td><td>{basis}</td>"
                f"<td class='n'>{_num(topo_val, 3)}</td>"
                f"<td class='n'>{_num(fold_val, 3)}</td>"
                f"<td class='n'>{_num(f.get('brier') if held_out else None, 3)}</td>"
                f"<td class='n'>{_num(f.get('ece') if held_out else None, 3)}</td></tr>"
            )
            chart_rows.append((b.get("model", "?")[:26],
                               float(topo_val or 0.0), float(fold_val or 0.0)))
        bench_chart = _paired_bars(
            chart_rows, "Scored on the same real measurements",
            "topology balanced accuracy", "folding AUROC",
            "Chance balanced accuracy is 0.333; chance AUROC is 0.5. Where a model "
            "trained on these records, its grouped out-of-fold value is plotted, "
            "never its in-sample one.",
        )
        bench_html = f"""
  <section id="benchmark">
    <div class="head"><h2>External benchmark</h2>
      <p>Every model below is scored on the same held-out set: the atlas's experimental,
      biophysical G4 topology measurements, and those same sequences against
      composition-matched dinucleotide shuffles. The script reconstructs each model's
      training set first and refuses to present an in-sample result as a held-out one.
      Reproduce with <code>scripts/05_external_benchmark.py</code>.</p></div>
    {bench_chart}
    <table><thead><tr><th>model</th><th>basis</th><th class="n">topology bal. acc.</th>
    <th class="n">folding AUROC</th><th class="n">folding Brier</th>
    <th class="n">folding ECE</th></tr></thead>
    <tbody>{''.join(rows_b)}</tbody></table>
    <p class="figNote">A model trained on synthetic data can report ordinary-looking
    metrics on its own held-out synthetic split and still land at chance here. That is the
    reason evidence tier, label class and dataset fingerprint are stored per row rather
    than described in a README.</p>
  </section>"""

    grounded = [n for n, h in model.heads.items()
                if claims.target_semantics(n, h.training_meta) == claims.BIOPHYSICAL]
    ungrounded = [n for n, h in model.heads.items()
                  if claims.target_semantics(n, h.training_meta) != claims.BIOPHYSICAL]

    # ---- ablation rows
    ab_rows = []
    for name, h in model.heads.items():
        if seq_only is None or name not in seq_only.heads:
            continue
        key = KEY_METRIC[h.task]
        va, vb = h.metrics.get(key), seq_only.heads[name].metrics.get(key)
        if va is None or vb is None:
            continue
        ab_rows.append((f"{name} \u00b7 {METRIC_LABEL[key]}", float(va), float(vb), name, key))
    best_delta = max((a - b for _, a, b, _, _ in ab_rows), default=0.0)

    # ---- headline tiles
    n_exp = int(summary.loc[summary["evidence_tier"] == "experimental", "n"].sum())
    tiles = [
        ("Atlas records", f"{total:,}", f"{n_exp:,} experimental"),
        ("Trained heads", f"{len(model.heads)}", f"{len(grounded)} measurement-grounded"),
        ("Distinct buffers", f"{int(summary['n_conditions'].max())}",
         "ionic conditions in the largest source"),
        ("Best ablation gain", f"{best_delta:+.3f}", "from adding condition features"),
    ]
    tiles_html = "".join(
        f'<div class="tile"><div class="k">{_e(k)}</div><div class="v">{_e(v)}</div>'
        f'<div class="s">{_e(sub)}</div></div>' for k, v, sub in tiles
    )

    # ---- atlas charts
    tier_totals = summary.groupby("evidence_tier")["n"].sum()
    order = [t for t in ("experimental", "derived", "predicted") if t in tier_totals.index]
    comp = _stacked_composition(
        [(t, int(tier_totals[t]), TIER_SLOT[t]) for t in order],
        total, "Records by evidence tier",
    )
    src_bars = _log_bars(
        [(f"{r['source']} · {r['kind']}", int(r["n"]), TIER_SLOT.get(r["evidence_tier"], SERIES[0]))
         for _, r in summary.iterrows()],
        "Records per source",
        "Log scale — the predicted tier outnumbers the experimental tier by two orders of "
        "magnitude, which is the whole problem this project is about.",
    )
    atlas_tbl = "".join(
        f"<tr><td><code>{_e(r['source'])}</code></td><td>{_e(r['kind'])}</td>"
        f"<td><span class='pill'><i style='background:{TIER_SLOT.get(r['evidence_tier'], SERIES[0])}'></i>"
        f"{_e(r['evidence_tier'])}</span></td>"
        f"<td><span class='pill'>{_e(r.get('label_class', '—'))}</span></td>"
        f"<td class='n'>{int(r['n']):,}</td><td class='n'>{int(r['n_unique_seq']):,}</td>"
        f"<td class='n'>{int(r['n_conditions'])}</td><td class='n'>{int(r['n_topology']):,}</td>"
        f"<td class='n'>{int(r['n_tm']):,}</td><td class='n'>{int(r['n_pht']):,}</td></tr>"
        for _, r in summary.iterrows()
    )

    # ---- ablation section
    if ab_rows:
        ab_chart = _paired_bars(
            [(r[0], r[1], r[2]) for r in ab_rows],
            "Head performance with and without condition features",
            "sequence + conditions", "sequence only",
            "Only the head whose training rows span more than one buffer moves. "
            "The flat rows are a measurement about the data, not about the method.",
        )
        ab_tbl = "".join(
            f"<tr><td><code>{_e(nm)}</code></td><td>{_e(METRIC_LABEL[k])}</td>"
            f"<td class='n'>{a:.4f}</td><td class='n'>{b:.4f}</td>"
            f"<td class='n'><span class='delta {'up' if a - b > 0.02 else 'flat'}'>"
            f"{a - b:+.4f}</span></td></tr>"
            for _, a, b, nm, k in ab_rows
        )
        ablation_html = f"""{ab_chart}
        <table><thead><tr><th>head</th><th>metric</th><th class="n">with conditions</th>
        <th class="n">sequence only</th><th class="n">Δ</th></tr></thead><tbody>{ab_tbl}</tbody></table>"""
    else:
        ablation_html = '<p class="figNote">No sequence-only model found; run <code>scripts/02_train.py</code>.</p>'

    # ---- case studies
    case_html = ""
    if cases:
        sweep = cases.get("sweeps", {}).get("hTelo_K_titration", [])
        if sweep:
            xs = [r["k"] for r in sweep]
            keys = [k for k in sweep[0] if k not in ("k", "in_domain")
                    and isinstance(sweep[0][k], (int, float))]
            tm_keys = [k for k in keys if "tm" in k and "folded" not in k]
            prob_keys = [k for k in keys if k.startswith("g4_topology")]
            if tm_keys:
                case_html += _line_chart(
                    xs, [(k, [r[k] for r in sweep]) for k in tm_keys],
                    "[K⁺] (mM)", "predicted Tm (°C)",
                    "Human telomeric G4 (AGGGTTAGGGTTAGGGTTAGGG): Tm across a K⁺ titration",
                    "The response is real because the underlying rows span five ionic "
                    "conditions. Heads trained on a single buffer produce a flat line here.",
                )
            if prob_keys:
                case_html += _line_chart(
                    xs, [(k.split(":")[-1], [r[k] for r in sweep]) for k in prob_keys],
                    "[K⁺] (mM)", "posterior probability",
                    "Topology posterior across the same K⁺ titration",
                    "Near-flat, and correctly so: the topology head's training rows carry no "
                    "reported buffer at all, so it has no basis to move.",
                )
        rows = []
        for c in cases.get("cases", [])[:8]:
            k100 = c["by_condition"].get("100 mM K+", {})
            na100 = c["by_condition"].get("100 mM Na+", {})
            topo = k100.get("g4_topology")
            rows.append(
                f"<tr><td>{_e(c['name'])}<div class='seq'>{_e(c['sequence'])}</div></td>"
                f"<td>{_e(topo['call']) if isinstance(topo, dict) else '—'}</td>"
                f"<td class='n'>{_num(k100.get('g4_fold'), 3)}</td>"
                f"<td class='n'>{_num(k100.get('g4_tm_distilled'), 1)}</td>"
                f"<td class='n'>{_num(na100.get('g4_tm_distilled'), 1)}</td>"
                f"<td class='n'>{_num((k100.get('g4_tm_distilled') or 0) - (na100.get('g4_tm_distilled') or 0), 1)}</td>"
                f"</tr>"
            )
        if rows:
            case_html += f"""<table><thead><tr><th>construct</th><th>topology call</th>
            <th class="n">P(G4)</th><th class="n">Tm in K⁺</th><th class="n">Tm in Na⁺</th>
            <th class="n">ΔTm</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
            <p class="figNote">Tm columns come from the distilled head and are not measurements.
            The K⁺ − Na⁺ difference is the quantity worth reading: it is the only column here
            that could not be produced by a condition-blind model.</p>"""

    # ---- head cards
    cards, nav_heads = [], []
    for name, h in model.heads.items():
        tm_, mt = h.training_meta, h.metrics
        _sem = claims.target_semantics(name, tm_)
        is_grounded = _sem == claims.BIOPHYSICAL
        nav_heads.append(f'<a class="sub" href="#h-{_e(name)}">{_e(name)}</a>')
        badge = (f'<span class="pill good">biophysically anchored</span>' if is_grounded
                 else f'<span class="pill warn">{_e(claims.SEMANTICS_LABEL[_sem])}</span>')

        if h.task == "binary":
            kv = [("AUROC", mt.get("auroc")), ("AUPRC", mt.get("auprc")),
                  ("Brier", mt.get("brier")), ("ECE (calibrated)", mt.get("ece")),
                  ("ECE (uncalibrated)", mt.get("ece_uncalibrated")),
                  ("calibrator", mt.get("calibration_method")),
                  ("accuracy @ 0.5", mt.get("accuracy_at_0.5"))]
            chart = _reliability(mt.get("reliability", {}), f"{name} — reliability")
        elif h.task == "multiclass":
            kv = [("balanced accuracy", mt.get("balanced_accuracy")),
                  ("accuracy", mt.get("accuracy")),
                  ("AUROC (OvR macro)", mt.get("auroc_ovr_macro")),
                  ("Brier", mt.get("brier")), ("ECE (calibrated)", mt.get("ece")),
                  ("ECE (uncalibrated)", mt.get("ece_uncalibrated")),
                  ("temperature", mt.get("temperature"))]
            chart = _reliability(mt.get("reliability", {}), f"{name} — reliability (top label)")
        else:
            lo, hi = mt.get("target_range", [None, None])
            kv = [("R²", mt.get("r2")), ("RMSE", mt.get("rmse")), ("MAE", mt.get("mae")),
                  ("Spearman ρ", mt.get("spearman")),
                  (f"conformal ± @{1-mt.get('conformal_alpha', .1):.0%}",
                   mt.get("conformal_halfwidth")),
                  ("empirical coverage", mt.get("empirical_coverage"))]
            chart = (f'<figure class="fig"><figcaption>target range in training</figcaption>'
                     f'<p class="figNote mono">{_num(lo, 2)} to {_num(hi, 2)}</p></figure>'
                     if lo is not None else "")

        kv_html = "".join(f"<tr><th>{_e(k)}</th><td class='n'>{_num(v)}</td></tr>" for k, v in kv)
        if mt.get("per_class"):
            kv_html += "<tr><th colspan='2' style='padding-top:12px'>per class</th></tr>"
            for c, d in mt["per_class"].items():
                kv_html += (f"<tr><th>{_e(c)} <span class='mono' style='color:var(--muted)'>"
                            f"n={d['n']}</span></th><td class='n'>{_num(d.get('recall'), 3)}</td></tr>")

        g = tm_.get("grouping", {})
        applic = "".join(
            f"<tr><td>{_e(k)}</td><td class='n'>{v['min']:g} – {v['max']:g}</td>"
            f"<td class='n'>{v['n_unique']}</td></tr>"
            for k, v in (h.applicability or {}).items()
        )
        frozen = [k for k, v in (h.applicability or {}).items() if v["n_unique"] <= 1]
        frozen_note = (f"<p class='figNote'><b>Single value seen for:</b> "
                       f"{', '.join(f'<code>{_e(k)}</code>' for k in frozen)}. "
                       f"Varying those at prediction time is extrapolation, and the predictor "
                       f"flags it.</p>" if frozen else "")

        cards.append(f"""
<article class="headcard{'' if is_grounded else ' ungrounded'}" id="h-{_e(name)}">
  <div class="hc"><h3><code>{_e(name)}</code></h3>
    <div>{badge}<span class="pill">{_e(h.task)}</span><span class="pill">{_e(h.kind)}</span>
      <span class="pill">{tm_.get('n_rows', 0):,} rows</span></div></div>
  <p class="claim">{_e(tm_.get('claim', ''))}</p>
  <div class="hcbody">
    <div><table class="kv">{kv_html}</table></div>
    <div>{chart}</div>
  </div>
  <details><summary>Training provenance and applicability domain</summary>
    <table class="kv">
      <tr><th>evidence tiers</th><td>{_e(', '.join(tm_.get('tiers', [])))}</td></tr>
      <tr><th>sources</th><td class="mono">{_e(', '.join(tm_.get('sources', [])))}</td></tr>
      <tr><th>sequence groups</th><td class="n">{g.get('n_groups', 0):,} from
        {g.get('n_items', 0):,} ({g.get('largest_group', 0)} in the largest)</td></tr>
      <tr><th>cross-validation</th><td class="n">{tm_.get('n_folds', 0)} grouped folds ×
        {tm_.get('n_seeds', 0)} seeds</td></tr>
      <tr><th>features</th><td class="n">{tm_.get('n_features', 0)}</td></tr>
    </table>
    <table><thead><tr><th>condition variable</th><th class="n">range seen</th>
      <th class="n">distinct values</th></tr></thead><tbody>{applic}</tbody></table>
    {frozen_note}
  </details>
</article>""")

    src_html = "".join(
        f'<div class="src"><div><span class="pill"><i style="background:'
        f'{TIER_SLOT.get(s["evidence_tier"], SERIES[0])}"></i>{_e(s["evidence_tier"])}</span></div>'
        f'<div><code>{_e(s["source"])}</code><div>{_e(s["title"] or "")}</div>'
        + (f'<div class="meta">doi:{_e(s["doi"])}</div>' if s.get("doi") else "")
        + (f'<div class="meta">{_e(s["notes"])}</div>' if s.get("notes") else "")
        + "</div></div>"
        for s in sources
    )

    ungrounded_note = ""
    if ungrounded:
        ungrounded_note = f"""<div class="note"><b>Read the tiers before the numbers.</b>
        {', '.join(f'<code>{_e(n)}</code>' for n in ungrounded)}
        {'were' if len(ungrounded) > 1 else 'was'} trained on derived or model-predicted labels,
        because the corresponding experimental tables are not in this atlas. They are
        architecture and distillation heads, not folding predictors. Ingesting G4STAB
        Supplementary Table 1 and the iM-Seeker pH_T tables and retraining is what changes
        that.</div>"""

    doc = f"""<title>QuadCond Model Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>{CSS}</style>
<div class="shell">
<nav class="rail">
  <div class="brand">QuadCond</div>
  <div class="tag">model report</div>
  <a href="#atlas">Atlas</a>
  <a href="#ablation">Condition ablation</a>
  <a href="#cases">Case studies</a>
  <a href="#heads">Heads</a>
  <a href="#sources">Sources</a>
  <a href="#raw">Raw metadata</a>
  <div class="grp">heads</div>
  {''.join(nav_heads)}
</nav>
<main>
  <h1>Condition-aware G-quadruplex and i-motif prediction</h1>
  <p class="lede">Every score on this page refers to a buffer. Cation identity and
  concentration, pH and temperature are model inputs, probabilities are calibrated on
  out-of-fold predictions, and cross-validation is grouped by sequence similarity so
  near-duplicates cannot straddle a split.</p>
  <p class="stamp">Generated {stamp} · atlas {total:,} records · {len(model.heads)} heads</p>
  <div class="tiles">{tiles_html}</div>
  {ungrounded_note}

  <section id="atlas">
    <div class="head"><h2>Atlas</h2>
      <p>One row is one measurement of one sequence under one buffer. Every row declares
      two independent provenance axes — <em>evidence tier</em> (was it measured, derived or
      predicted) and <em>label class</em> (does it measure a structure in a buffer, occupancy
      in a nucleus, or just a motif call) — plus its method, its DOI, and which condition
      fields were reported rather than imputed. A ChIP peak and a melting curve are both
      experimental; they are not the same evidence, and the atlas will not pool them.</p></div>
    {comp}
    {src_bars}
    <table><thead><tr><th>source</th><th>kind</th><th>tier</th><th>label class</th><th class="n">records</th>
    <th class="n">unique seq</th><th class="n">buffers</th><th class="n">topology</th>
    <th class="n">Tm</th><th class="n">pH_T</th></tr></thead><tbody>{atlas_tbl}</tbody></table>
  </section>

  <section id="ablation">
    <div class="head"><h2>Does explicit conditioning help?</h2>
      <p>Same data, same protocol, condition and interaction features removed. This is the
      claim of the project, tested.</p></div>
    {ablation_html}
  </section>

  <section id="cases">
    <div class="head"><h2>Case studies</h2>
      <p>Constructs the field knows well, scored across buffers. Not a benchmark — a sanity
      surface. Behaviour that contradicts the primary literature here would invalidate
      everything else on this page.</p></div>
    {case_html or '<p class="figNote">Run <code>scripts/04_case_studies.py</code> to populate this section.</p>'}
  </section>

  <section id="heads">
    <div class="head"><h2>Heads</h2>
      <p>Each head states the claim it is entitled to make, the tier of data behind it, and the
      region of condition space it actually saw.</p></div>
    {''.join(cards)}
  </section>

  {bench_html}

  <section id="sources">
    <div class="head"><h2>Sources</h2>
      <p>Cite these, not the framework. Full provenance, including the tables that require a
      manual download, is in <code>docs/DATA_SOURCES.md</code>.</p></div>
    {src_html}
  </section>

  <section id="raw">
    <div class="head"><h2>Raw model metadata</h2></div>
    <details><summary>Show the full JSON summary</summary>
    <pre>{_e(json.dumps(model.summary(), indent=2, default=str))}</pre></details>
  </section>
</main>
</div>"""

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding="utf-8")
    return out
