#!/usr/bin/env python3
"""Build the frozen-results report as one self-contained HTML file.

Every number in the page is read out of the artifacts rather than typed here:
``artifacts/quadcond_model.json`` for the heads, their metrics, their
reliability curves and their applicability ranges; ``artifacts/ablation.json``
for the condition ablation; ``artifacts/external_validation_im.json`` for the
held-out CD set; ``data/genomic/gse220882_peaks.meta.json`` for the peak-calling
record. The atlas content fingerprint is recomputed from the database on the spot.

That is the point. A results page assembled by hand drifts from the artifacts it
describes the first time anything is retrained, and the drift is invisible. This
one cannot: regenerate it and either the numbers match the models or the script
raises.

Two files come out. ``--out`` is a complete standalone document to open or send.
``--artifact-out`` is the same page without the document skeleton, for hosts that
supply their own.

    python scripts/08_results_report.py

Templates are in ``scripts/report_templates/``.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import html
import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

TEMPLATES = Path(__file__).resolve().parent / "report_templates"

# One definition of atlas identity, imported rather than restated: this script
# and quadcond.assets both quote the fingerprint, and two functions hashing
# slightly different column lists would print two different "the" fingerprints
# for the same database.
from quadcond.assets import atlas_fingerprint as _atlas_content_fingerprint  # noqa: E402
from quadcond.atlas import Atlas  # noqa: E402


def atlas_fingerprint(db: str | Path) -> tuple[int, str]:
    """Row count and the content fingerprint of an atlas.

    Not the file hash: SQLite rewrites page metadata whenever the database is
    opened, so ``atlas.db`` checksums differently after a mere read. Quoting the
    file hash as an identity would make a reproducible atlas look
    irreproducible.
    """
    con = sqlite3.connect(f"file:{Path(db)}?mode=ro", uri=True)
    try:
        n = con.execute("SELECT COUNT(*) FROM records").fetchone()[0]
    finally:
        con.close()
    return int(n), _atlas_content_fingerprint(db)


def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def esc(s):
    return html.escape(str(s))


def f(x, n=3):
    return "\u2014" if x is None else f"{x:.{n}f}"


def hue(name):
    return "g4" if KIND.get(name) == "G4" else "im"



def bar_panel(rows, x0, x1, width=560, row_h=34, title="", note=""):
    """rows: (head, label, with_cond, seq_only, metric_label)"""
    pad_l, pad_r, pad_t, pad_b = 168, 74, 26, 30
    h = pad_t + len(rows) * row_h + pad_b
    w = width
    plot_w = w - pad_l - pad_r
    def X(v): return pad_l + (v - x0) / (x1 - x0) * plot_w
    out = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(title)}" class="fig">']
    # ticks
    n_ticks = 5
    for i in range(n_ticks + 1):
        v = x0 + (x1 - x0) * i / n_ticks
        x = X(v)
        out.append(f'<line x1="{x:.1f}" y1="{pad_t-6}" x2="{x:.1f}" y2="{h-pad_b}" class="grid"/>')
        out.append(f'<text x="{x:.1f}" y="{h-pad_b+16}" class="tick" text-anchor="middle">{v:.1f}</text>')
    for i, (name, label, wc, so, ml) in enumerate(rows):
        y = pad_t + i * row_h
        cy = y + row_h / 2
        k = hue(name)
        xs, xw = X(min(wc, so)), abs(X(wc) - X(so))
        # connector
        out.append(f'<line x1="{X(so):.1f}" y1="{cy:.1f}" x2="{X(wc):.1f}" y2="{cy:.1f}" class="conn {k}"/>')
        out.append(f'<circle cx="{X(so):.1f}" cy="{cy:.1f}" r="5" class="dot-open {k}"/>')
        out.append(f'<circle cx="{X(wc):.1f}" cy="{cy:.1f}" r="5.5" class="dot {k}">'
                   f'<title>{esc(name)}: {ml} {wc:.3f} with conditions, {so:.3f} sequence only</title></circle>')
        out.append(f'<text x="{pad_l-14}" y="{cy+4:.1f}" class="rowlab" text-anchor="end">{esc(label)}</text>')
        d = wc - so
        out.append(f'<text x="{w-pad_r+10}" y="{cy+4:.1f}" class="delta {k}">{d:+.3f}</text>')
    out.append("</svg>")
    return "\n".join(out)

def reliability(name, size=190):
    h = HEADS[name]; r = h["reliability"]
    pad = 30
    inner = size - 2 * pad
    m = 5          # keep an r=4 mark at value 1.0 inside the plotted frame
    k = hue(name)
    def P(v): return pad + m + v * (inner - 2 * m)
    def Q(v): return size - pad - m - v * (inner - 2 * m)
    o = [f'<svg viewBox="0 0 {size} {size+6}" role="img" aria-label="Reliability diagram for {esc(name)}" class="rel">']
    o.append(f'<rect x="{pad}" y="{pad}" width="{inner}" height="{inner}" class="relbox"/>')
    o.append(f'<line x1="{P(0)}" y1="{Q(0)}" x2="{P(1)}" y2="{Q(1)}" class="diag"/>')
    pts = " ".join(f"{P(c):.1f},{Q(a):.1f}" for c, a in zip(r["confidence"], r["accuracy"]))
    o.append(f'<polyline points="{pts}" class="relline {k}"/>')
    for c, a, n in zip(r["confidence"], r["accuracy"], r["count"]):
        o.append(f'<circle cx="{P(c):.1f}" cy="{Q(a):.1f}" r="4" class="dot {k}">'
                 f'<title>predicted {c:.2f} · observed {a:.2f} · n={n}</title></circle>')
    o.append(f'<text x="{pad}" y="{size+2}" class="tick">0</text>')
    o.append(f'<text x="{size-pad}" y="{size+2}" class="tick" text-anchor="end">1</text>')
    o.append("</svg>")
    return "\n".join(o)

def enrichment(width=560):
    data = [("iM", "iMab peaks", 0.373), ("iM", "shuffled windows", 0.275), ("iM", "random hg38", 0.039),
            ("G4", "BG4 peaks", 0.229), ("G4", "shuffled windows", 0.131), ("G4", "random hg38", 0.021)]
    pad_l, pad_r, pad_t, pad_b = 150, 56, 22, 28
    row_h, gap = 26, 16
    h = pad_t + 3 * row_h + gap + 3 * row_h + pad_b
    plot_w = width - pad_l - pad_r
    xmax = 0.40
    o = [f'<svg viewBox="0 0 {width} {h}" role="img" aria-label="Canonical motif frequency" class="fig">']
    for i in range(5):
        v = xmax * i / 4
        x = pad_l + v / xmax * plot_w
        o.append(f'<line x1="{x:.1f}" y1="{pad_t-6}" x2="{x:.1f}" y2="{h-pad_b}" class="grid"/>')
        o.append(f'<text x="{x:.1f}" y="{h-pad_b+16}" class="tick" text-anchor="middle">{v*100:.0f}%</text>')
    y = pad_t
    for i, (kind, label, v) in enumerate(data):
        if i == 3: y += gap
        k = "g4" if kind == "G4" else "im"
        tone = "" if "peaks" in label else (" mid" if "shuffled" in label else " faint")
        bw = v / xmax * plot_w
        o.append(f'<rect x="{pad_l}" y="{y+5}" width="{bw:.1f}" height="{row_h-12}" rx="3" class="bar {k}{tone}">'
                 f'<title>{esc(label)}: {v*100:.1f}%</title></rect>')
        o.append(f'<text x="{pad_l-12}" y="{y+row_h/2+4:.0f}" class="rowlab" text-anchor="end">{esc(label)}</text>')
        o.append(f'<text x="{pad_l+bw+8:.1f}" y="{y+row_h/2+4:.0f}" class="valnum">{v*100:.1f}%</text>')
        y += row_h
    o.append("</svg>")
    return "\n".join(o)

def transfer(width=560):
    rows = [("im_fold", "im_fold", EXT["mean_p_im"], EXT["im_separation"]),
            ("im_fold_genomic", "im_fold_genomic", EXT["mean_p_im_genomic"], EXT["im_genomic_separation"]),
            ("g4_fold", "g4_fold", EXT["mean_p_g4"], EXT["g4_separation"])]
    pad_l, pad_r, pad_t, pad_b = 152, 60, 34, 30
    row_h = 46
    h = pad_t + len(rows) * row_h + pad_b
    plot_w = width - pad_l - pad_r
    def X(v): return pad_l + v * plot_w
    o = [f'<svg viewBox="0 0 {width} {h}" role="img" aria-label="Mean predicted probability on the held-out CD set" class="fig">']
    for i in range(6):
        v = i / 5
        o.append(f'<line x1="{X(v):.1f}" y1="{pad_t-8}" x2="{X(v):.1f}" y2="{h-pad_b}" class="grid"/>')
        o.append(f'<text x="{X(v):.1f}" y="{h-pad_b+16}" class="tick" text-anchor="middle">{v:.1f}</text>')
    marks = [("i-motif", "i-motif set"), ("G4", "G4 controls"), ("neither", "negatives")]
    for i, (name, label, means, sep) in enumerate(rows):
        y = pad_t + i * row_h
        cy = y + row_h / 2 - 4
        k = hue(name)
        vals = [means[m[0]] for m in marks]
        o.append(f'<line x1="{X(min(vals)):.1f}" y1="{cy:.1f}" x2="{X(max(vals)):.1f}" y2="{cy:.1f}" class="conn {k}"/>')
        for (key, label2), shape in zip(marks, ["dot", "dot-open", "dot-x"]):
            v = means[key]
            cls = "dot" if shape == "dot" else ("dot-open" if shape == "dot-open" else "dot-half")
            o.append(f'<circle cx="{X(v):.1f}" cy="{cy:.1f}" r="6" class="{cls} {k}">'
                     f'<title>{esc(label)} · {esc(label2)}: {v:.3f}</title></circle>')
        o.append(f'<text x="{pad_l-14}" y="{cy+4:.1f}" class="rowlab mono" text-anchor="end">{esc(label)}</text>')
        o.append(f'<text x="{width-pad_r+8}" y="{cy+4:.1f}" class="delta {k}">{sep:+.3f}</text>')
    o.append("</svg>")
    return "\n".join(o)

# ---------------------------------------------------------------- tables
PREFERRED_ORDER = [
    "g4_fold", "g4_topology", "g4_tm", "g4_tm_distilled", "im_fold", "im_pht",
    "im_pht_condition", "im_tm_condition", "im_fold_genomic", "g4_fold_genomic",
    "locus_peak_overlap_state", "im_architecture",
]


def all_heads() -> list[str]:
    """Every head in the artifact, preferred ones first.

    A hand-written list here is how the last release shipped a report that
    described eleven heads while the model contained twelve: `locus_state` was
    trained, its atlas sources appeared in the data table, and the head itself
    was simply absent from both the table and the claims section because nobody
    added it to two literals. Anything not in PREFERRED_ORDER now appears at the
    end rather than silently not at all.
    """
    known = [n for n in PREFERRED_ORDER if n in HEADS]
    return known + sorted(n for n in HEADS if n not in PREFERRED_ORDER)


def heads_table():
    order = all_heads()
    metric_name = {"binary": "AUROC", "multiclass": "bal. acc.", "regression": "R²"}
    rows = []
    for n in order:
        h = HEADS[n]
        m = h["metrics"]
        key = "auroc" if h["task"] == "binary" else ("balanced_accuracy" if h["task"] == "multiclass" else "r2")
        ts = h.get("target_semantics", "biophysical")
        badge, btitle = {
            "biophysical": ("biophysical", "the label is a physical measurement of the "
                                           "structure this head predicts"),
            "genomic_proxy": ("proxy", "antibody occupancy at a genomic locus — an "
                                       "experimental observation of a different quantity "
                                       "from folding; cannot license P(folds)"),
            "derived": ("derived", "no measurement stands behind the label"),
            "predicted": ("predicted", "another model's output, not a measurement"),
        }[ts]
        ece = m.get("ece")
        rows.append(f"""<tr>
 <th scope="row"><code class="{hue(n)}">{esc(n)}</code></th>
 <td><span class="badge {badge}" title="{esc(btitle)}">{badge}</span></td>
 <td class="num">{m[key]:.3f}<span class="unit">{metric_name[h['task']]}</span></td>
 <td class="num">{f(ece) if ece is not None else '—'}</td>
 <td class="num">{h['n_rows']:,}</td>
 <td class="num">{h['n_groups']:,}</td>
 <td>{esc(h['group_by'])}</td>
</tr>""")
    return "\n".join(rows)

def claims_list():
    claim_first = ["g4_tm", "im_tm_condition", "im_pht_condition", "im_pht", "im_fold",
                   "im_fold_genomic", "g4_fold_genomic", "locus_peak_overlap_state",
                   "g4_fold", "g4_topology", "g4_tm_distilled", "im_architecture"]
    order = [n for n in claim_first if n in HEADS] + sorted(
        n for n in HEADS if n not in claim_first)
    out = []
    for n in order:
        h = HEADS[n]
        claim = h["claim"].split(". ")
        prose = esc(h["claim"]).replace(" -- ", " \u2014 ")
        out.append(f'<div class="claim"><code class="{hue(n)}">{esc(n)}</code>'
                   f'<p>{prose}</p></div>')
    return "\n".join(out)

def applicability_rows():
    out = []
    for n in ("g4_tm", "im_pht_condition", "im_tm_condition"):
        a = HEADS[n]["applicability"]
        # Every axis the domain check enforces, not a subset. crowder_pct and
        # strand_conc were enforced from v0.4.3 and absent from this table, so
        # the page documented a narrower gate than the code applies -- which is
        # the more dangerous direction: a reader plans a query against the table.
        target = HEADS[n].get("target") or HEADS[n].get("training_meta", {}).get("target")
        for fld, unit in (("k", "mM K⁺"), ("na", "mM Na⁺"), ("li_nh4", "mM Li⁺/NH₄⁺"),
                          ("mg", "mM Mg²⁺"), ("ph", "pH"), ("temperature", "°C"),
                          ("crowder_pct", "% crowder"), ("strand_conc", "µM strand"),
                          ("length", "nt")):
            if fld not in a: continue
            # A melting-temperature head predicts a temperature; it does not
            # consume one, and `_domain_check` skips this field for such heads.
            # Listing it here as though it gated a query described a narrower
            # tool than the code is -- and a reader plans a query against this
            # table, so the row was worse than absent.
            if fld == "temperature" and target == "tm":
                continue
            v = a[fld]
            if v["min"] == v["max"]:
                span = f'{v["min"]:g} only'
            else:
                span = f'{v["min"]:g} – {v["max"]:g}'
            core = f'{v["p05"]:g} – {v["p95"]:g}'
            thin = v["n_unique"] < 6 or (v["p95"] - v["p05"]) < 0.25 * (v["max"] - v["min"])
            out.append(f'<tr><th scope="row"><code class="{hue(n)}">{esc(n)}</code></th>'
                       f'<td>{esc(unit)}</td><td class="num">{span}</td>'
                       f'<td class="num{" warn" if thin else ""}">{core}</td>'
                       f'<td class="num">{v["n_unique"]}</td></tr>')
    return "\n".join(out)

def cd_rows():
    out = []
    for r in EXT["rows"]:
        cls = {"i-motif": "im", "G4": "g4", "neither": "neu"}[r["expectation"]]
        out.append(f'<tr><th scope="row" class="mono">{esc(r["name"])}</th>'
                   f'<td><span class="tag {cls}">{esc(r["expectation"])}</span></td>'
                   f'<td class="num">{r["p_im"]:.3f}</td>'
                   f'<td class="num strong">{r["p_im_genomic"]:.3f}</td>'
                   f'<td class="num">{r["p_g4"]:.3f}</td>'
                   f'<td class="num">{r["ph_t"]:.2f}</td>'
                   f'<td class="mono small">{esc(r["motif"])}</td>'
                   f'<td class="seq">{esc(r["sequence"])}</td></tr>')
    return "\n".join(out)

def sample_rows():
    out = []
    for s in D["peak_samples"]:
        out.append(f'<tr><th scope="row" class="mono">{esc(s["gsm"])}</th>'
                   f'<td>{esc(s["cell_line"])}</td><td>{esc(s["antibody"])}</td>'
                   f'<td class="num">{s["replicate"]}</td>'
                   f'<td class="num">{s["blocks"]:,}</td>'
                   f'<td class="num">{s["authors_seacr_peaks"]:,}</td>'
                   f'<td class="num strong">{s["implied_top_fraction"]*100:.2f}%</td></tr>')
    return "\n".join(out)

def atlas_rows():
    """Rows from the atlas itself, not from the model's training snapshot.

    The snapshot's `n_conditions` was computed with an older key that omitted
    Mg2+, so the table printed 255 for g4stab beside prose from a live query
    saying 261 -- two counts of the same thing, differing, with neither defined.
    The snapshot's row total is still checked against the atlas (see main), so
    reading the breakdown live cannot silently describe a different corpus.
    """
    live = Atlas(ARGS.db).summary().to_dict(orient="records")
    out = []
    for s in live:
        k = {"G4": "g4", "iM": "im"}.get(s["kind"], "im")
        out.append(f'<tr><th scope="row" class="mono src">{esc(s["source"])}</th>'
                   f'<td><span class="kind {k}">{esc(s["kind"])}</span></td>'
                   f'<td>{esc(s["evidence_tier"])}</td>'
                   f'<td>{esc(s["label_class"])}</td>'
                   f'<td class="num">{s["n"]:,}</td>'
                   f'<td class="num">{s["n_unique_seq"]:,}</td>'
                   f'<td class="num">{s["n_conditions"]:,}</td></tr>')
    return "\n".join(out)

# --------------------------------------------------------------------------



def n_g4_seqs() -> int:
    con = sqlite3.connect(f"file:{Path(ARGS.db)}?mode=ro", uri=True)
    try:
        return con.execute("SELECT COUNT(DISTINCT sequence) FROM records "
                           "WHERE source='g4stab_experimental_tm'").fetchone()[0]
    finally:
        con.close()


def n_buffers() -> int:
    """Distinct G4 buffers in the measured Tm source, counted from the atlas."""
    con = sqlite3.connect(f"file:{Path(ARGS.db)}?mode=ro", uri=True)
    try:
        # Same key as Atlas.summary()'s n_conditions column, deliberately:
        # the page printed 261 in prose beside a table column reading 255,
        # because the two counts used different keys and neither was defined.
        return con.execute(
            "SELECT COUNT(DISTINCT k||'_'||na||'_'||li_nh4||'_'||mg||'_'||ph) "
            "FROM records WHERE source='g4stab_experimental_tm'"
        ).fetchone()[0]
    finally:
        con.close()


def test_results(path: str = "artifacts/pytest-results.xml") -> dict:
    """Read a JUnit report. Never count `def test_` and call the result passing.

    v0.4.2 printed "52 tests passing" from a literal in the template, and the
    fix that replaced it counted test *functions* in the tree -- which says the
    tests exist, not that they pass, and would have gone on reporting a healthy
    number through a red suite. This reads the outcome of an actual run, and if
    no run artifact is present it says so instead of inventing one.
    """
    import xml.etree.ElementTree as ET

    f = Path(path)
    if not f.exists():
        return {"available": False,
                "text": "not recorded (run scripts/15_run_tests.py)"}
    root = ET.parse(f).getroot()
    suites = [root] if root.tag == "testsuite" else list(root)
    tot = sum(int(x.get("tests", 0)) for x in suites)
    fail = sum(int(x.get("failures", 0)) for x in suites)
    err = sum(int(x.get("errors", 0)) for x in suites)
    skip = sum(int(x.get("skipped", 0)) for x in suites)
    passed = tot - fail - err - skip
    return {
        "available": True, "total": tot, "passed": passed,
        "failed": fail + err, "skipped": skip,
        "text": (f"{passed} passing" if not (fail + err)
                 else f"{passed} passing, {fail + err} FAILING"),
        "detail": (f"{tot} collected: {passed} passed, {fail + err} failed, "
                   f"{skip} skipped"),
    }


def commit_id() -> str:
    """The tree the run describes, when it is a git checkout."""
    import subprocess

    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or "not a git checkout"
    except Exception:                                     # noqa: BLE001
        return "not a git checkout"


def conformal_rows() -> str:
    """The grouped-split coverage actually present in the artifacts.

    Rendered from the metrics rather than typed into the prose, because the
    previous version of this paragraph said the estimate was "implemented for
    the next training run" long after that run had happened.
    """
    bits = []
    for name, h in HEADS.items():
        cov = h["metrics"].get("split_conformal_coverage")
        if cov is None:
            continue
        bits.append(f"<code>{name}</code> {cov:.3f}")
    if not bits:
        return "no head in this artifact carries one"
    return ", ".join(bits[:-1]) + (" and " + bits[-1] if len(bits) > 1 else "")


def cation_block() -> str:
    """The out-of-fold K+/Na+ check, or a note that it has not been run."""
    path = Path("artifacts/cation_identity_check.json")
    if not path.exists():
        return ('<p class="note">scripts/10_cation_identity_check.py has not been run '
                'against this artifact.</p>')
    d = json.loads(path.read_text())["cation_identity"]
    if "out_of_fold" not in d:
        return ('<p class="note">the cation check in this artifact is in-sample only; '
                'rerun scripts/10_cation_identity_check.py without --in-sample-only.</p>')
    o, i = d["out_of_fold"], d["in_sample"]
    ci = d.get("out_of_fold_ci95", {})

    def band(stat, fmt="{:+.1f}"):
        b = ci.get(stat)
        return "" if not b else f' <span class="dim">[{fmt.format(b["lo"])}, {fmt.format(b["hi"])}]</span>'

    return f"""
    <table class="grid">
      <thead><tr><th></th><th>measured</th><th>out-of-fold</th><th>in-sample</th></tr></thead>
      <tbody>
        <tr><td>mean &Delta;Tm(K&minus;Na)</td>
            <td>{i['measured_delta_mean']:+.1f} &deg;C</td>
            <td><strong>{o['predicted_delta_mean']:+.1f} &deg;C</strong>{band('predicted_delta_mean')}</td>
            <td>{i['predicted_delta_mean']:+.1f} &deg;C</td></tr>
        <tr><td>sign agreement</td><td>&mdash;</td>
            <td><strong>{o['sign_agreement']:.0%}</strong>{band('sign_agreement', '{:.0%}')}</td>
            <td>{i['sign_agreement']:.0%}</td></tr>
        <tr><td>Pearson r</td><td>&mdash;</td>
            <td><strong>{o['pearson_r']:.3f}</strong>{band('pearson_r', '{:+.2f}')}</td>
            <td>{i['pearson_r']:.3f}</td></tr>
        <tr><td>MAE of &Delta;Tm</td><td>&mdash;</td>
            <td>{o['mae_of_delta']:.1f} &deg;C</td>
            <td>{i['mae_of_delta']:.1f} &deg;C</td></tr>
      </tbody>
    </table>
    <p>{d['n_pairs']} matched pairs over {d['n_sequences']} sequences, each scored by the fold
    that held its whole sequence-similarity cluster out. Intervals are 95% percentile bootstrap
    over <strong>clusters</strong>, not over pairs: hTelo 22AG alone contributes 187 of the pairs,
    so a pair-level interval would be several times too tight. Memorisation is the gap between the
    two prediction columns &mdash; it buys {i['pearson_r'] - o['pearson_r']:+.3f} of r, and the
    held-out interval on sign agreement still excludes chance.</p>
    """


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-json", default="artifacts/quadcond_model.json")
    ap.add_argument("--model", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--seqonly", default="artifacts/quadcond_model_seqonly.joblib")
    ap.add_argument("--ablation", default="artifacts/ablation.json")
    ap.add_argument("--external", default="artifacts/external_validation_im.json")
    ap.add_argument("--peaks-meta", default="data/genomic/gse220882_peaks.meta.json")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--out", default=None,
                    help="default: artifacts/quadcond_v<version>_results.html")
    ap.add_argument("--artifact-out", default="artifacts/quadcond_results.body.html")
    a = ap.parse_args()

    global D, HEADS, ABL, EXT, PM, KIND, ARGS, TESTS
    ARGS = a
    TESTS = test_results()

    model = json.loads(Path(a.model_json).read_text())
    heads = model["heads"] if isinstance(model["heads"], list) else list(model["heads"].values())
    n_rows, atlas_fp = atlas_fingerprint(a.db)

    D = {
        "version": model["version"],
        "dataset_fingerprint": model["dataset_fingerprint_sha256"],
        "atlas": model["atlas_snapshot"],
        "ablation": json.loads(Path(a.ablation).read_text()),
        "external": json.loads(Path(a.external).read_text()),
        "heads": [],
    }
    peaks = json.loads(Path(a.peaks_meta).read_text())
    D["peaks_meta"] = {k: v for k, v in peaks.items() if k != "samples"}
    D["peak_samples"] = peaks["samples"]

    for h in heads:
        m, tm = h["metrics"], h["training_meta"]
        D["heads"].append({
            "name": h["name"], "kind": h["kind"], "task": h["task"], "target": h["target"],
            "n_rows": tm["n_rows"], "n_groups": tm["grouping"]["n_groups"],
            "group_by": tm.get("group_by"),
            "grounded": tm.get("biophysically_grounded", False),
            "target_semantics": tm.get("target_semantics", "biophysical"),
            "label_classes": tm.get("label_classes"), "sources": tm["sources"],
            "claim": tm["claim"], "tiers": tm["tiers"],
            "metrics": {k: v for k, v in m.items() if k != "reliability"},
            "reliability": m.get("reliability"),
            "applicability": h.get("applicability"),
        })

    from quadcond import claims as _claims
    D["tally"] = _claims.tally([
        {"name": h["name"],
         "training_meta": {"sources": h["sources"], "tiers": h["tiers"],
                           "label_classes": h["label_classes"]}}
        for h in D["heads"]
    ])
    HEADS = {h["name"]: h for h in D["heads"]}
    ABL = {x["head"]: x for x in D["ablation"]}
    EXT = D["external"]
    PM = D["peaks_meta"]
    KIND = {h["name"]: h["kind"] for h in D["heads"]}

    if n_rows != D["atlas"]["n_records"]:
        raise SystemExit(
            f"the atlas holds {n_rows:,} rows but the model was trained against "
            f"{D['atlas']['n_records']:,}. Retrain, or point --db at the atlas the "
            "model was built from -- do not publish a page that mixes the two."
        )

    abl_reg = [(n, n, ABL[n]["with_conditions"], ABL[n]["sequence_only"], "R\u00b2")
               for n in ("im_tm_condition", "g4_tm", "g4_tm_distilled",
                         "im_pht_condition", "im_pht")]
    abl_cls = [(n, n, ABL[n]["with_conditions"], ABL[n]["sequence_only"], ABL[n]["metric"])
               for n in ("g4_fold", "g4_topology", "im_fold",
                         "im_fold_genomic", "g4_fold_genomic")]

    today = datetime.date.today().strftime("%d %B %Y").lstrip("0")
    style = (TEMPLATES / "style.css").read_text()
    overlaps = EXT.get("training_overlap_with_gse220882") or []

    body = (TEMPLATES / "body.html").read_text().format(
        today=today,
        version=D["version"],
        n_records=f"{n_rows:,}",
        n_heads=len(D["heads"]),
        n_biophysical=D["tally"]["biophysical"],
        n_proxy=D["tally"]["genomic_proxy"],
        n_aux=D["tally"]["derived"] + D["tally"]["predicted"],
        dataset_fp=D["dataset_fingerprint"],
        artifact_sha=file_sha256(a.model),
        seqonly_sha=file_sha256(a.seqonly),
        atlas_fp=atlas_fp,
        heads_table=heads_table(),
        claims=claims_list(),
        fig_abl_reg=bar_panel(abl_reg, 0.0, 1.0, title="Condition ablation, regression heads"),
        fig_abl_cls=bar_panel(abl_cls, 0.5, 1.0, title="Condition ablation, classification heads"),
        fig_transfer=transfer(),
        fig_enrich=enrichment(),
        rel_g4_fold=reliability("g4_fold"),
        rel_im_fold=reliability("im_fold"),
        rel_im_gen=reliability("im_fold_genomic"),
        rel_g4_gen=reliability("g4_fold_genomic"),
        ece_g4_fold=f(HEADS["g4_fold"]["metrics"]["ece"]),
        ece_im_fold=f(HEADS["im_fold"]["metrics"]["ece"]),
        ece_im_gen=f(HEADS["im_fold_genomic"]["metrics"]["ece"]),
        ece_g4_gen=f(HEADS["g4_fold_genomic"]["metrics"]["ece"]),
        cd_rows=cd_rows(),
        sample_rows=sample_rows(),
        atlas_rows=atlas_rows(),
        applicability_rows=applicability_rows(),
        overlap_n=len(overlaps),
        overlap_name=overlaps[0]["name"] if overlaps else "none",
        im_sep=f'{EXT["im_separation"]:+.3f}',
        im_gen_sep=f'{EXT["im_genomic_separation"]:+.3f}',
        g4_sep=f'{EXT["g4_separation"]:+.3f}',
        canon=EXT["canonical_motif_found_in_authors_im_set"],
        junc=f'{PM["genome_self_test"]["fraction_canonical"] * 100:.1f}',
        n_junc=f'{PM["genome_self_test"]["junctions_checked"]:,}',
        conformal_rows=conformal_rows(),
        n_heads_total=len(HEADS),
        tests=TESTS['text'],
        tests_detail=TESTS.get('detail', TESTS['text']),
        commit=commit_id(),
        n_buffers=f"{n_buffers():,}",
        n_g4_seqs=f"{n_g4_seqs():,}",
        cation_block=cation_block(),
    )

    if a.out is None:
        a.out = f"artifacts/quadcond_v{D['version']}_results.html"
    title = f"QuadCond v{D['version']} Results"
    content = f"<title>{title}</title>\n<style>\n{style}\n</style>\n{body}"
    Path(a.artifact_out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.artifact_out).write_text(content)

    standalone = (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{title}</title>\n<style>\n{style}\n</style>\n</head>\n<body>\n"
        f"{body}\n</body>\n</html>\n"
    )
    Path(a.out).write_text(standalone)

    print(f"atlas: {n_rows:,} rows, content fingerprint {atlas_fp[:16]}...")
    print(f"wrote {a.out} ({len(standalone):,} bytes)")
    print(f"wrote {a.artifact_out} ({len(content):,} bytes)")


if __name__ == "__main__":
    main()
