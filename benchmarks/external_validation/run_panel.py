#!/usr/bin/env python3
"""Score a held-out panel through the frozen model, without touching it.

The benchmark in ../published_tools measures QuadCond against other tools on
collections QuadCond was trained on, using grouped out-of-fold predictions. That
is a fair internal protocol, but it is not external validation: every sequence
in it was available when the model was built.

A panel here is different. It is a set of sequences whose measurements were
published *after* the model artifact was frozen, scored by the frozen artifact
with no refitting of any kind. That is the only comparison that answers "does
this model work on data nobody could have tuned it to", and it is the claim the
manuscript currently cannot make.

The rules a panel has to satisfy, and why each one matters:

1. Its measurements postdate the model artifact. Record the artifact SHA-256 and
   the panel's publication date in the panel file and check them before use --
   a panel whose data predate the freeze is an ordinary test set, not an
   external one, and calling it external is the failure this directory exists
   to prevent.
2. Its conditions are recorded per row, not assumed. The point of this model is
   that buffer is an input; scoring a physiological mixed-cation panel at a
   defaulted 100 mM K+ measures nothing.
3. Nothing here is ever used for training or for selecting a model. If a panel
   is consulted while making a modelling decision it stops being held out, and
   it should be moved to ../published_tools and relabelled.

Usage:
    python run_panel.py panels/<panel>.tsv --model ../../artifacts/quadcond_model.joblib

Writes results/<panel>.csv (one row per sequence, with the applicability verdict
and refusals kept) and results/<panel>_summary.json (metrics, plus the model
identity the numbers came from).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(REPO))

REQUIRED = ["sequence_id", "sequence"]
CONDITION_COLS = {"k": "buffer_k_mm", "na": "buffer_na_mm", "li_nh4": "buffer_li_nh4_mm",
                  "mg": "buffer_mg_mm", "ph": "ph", "temperature": "temperature_c",
                  "crowder_pct": "crowder_pct", "strand_conc": "strand_conc_um"}
# Observations a panel may carry. A panel supplies whichever it has; metrics are
# computed only for those, so a folding-only panel does not silently report a
# melting-temperature correlation over an empty column.
OBSERVED = ["observed_folded", "observed_topology", "observed_tm_c", "observed_pht"]


def read_panel(path: Path):
    meta, rows = {}, []
    with path.open(encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    body = []
    for ln in lines:
        if ln.startswith("#"):
            if ":" in ln:
                k, v = ln[1:].split(":", 1)
                meta[k.strip()] = v.strip()
            continue
        if ln.strip():
            body.append(ln)
    if not body:
        return meta, []
    rdr = csv.DictReader(body, delimiter="\t")
    for r in rdr:
        rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()})
    for c in REQUIRED:
        if rows and c not in rows[0]:
            raise SystemExit(f"{path.name}: missing required column {c!r}")
    return meta, rows


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("--model", default=str(REPO / "artifacts" / "quadcond_model.joblib"))
    ap.add_argument("--atlas", default=str(REPO / "data" / "atlas_core.db"))
    ap.add_argument("--allow-undated", action="store_true",
                    help="score a panel whose measurements are not dated after the model freeze")
    a = ap.parse_args()

    path = Path(a.panel)
    meta, rows = read_panel(path)
    if not rows:
        raise SystemExit(
            f"{path.name} has a header and provenance block but no data rows yet.\n"
            "Extract the sequences and measurements from the source listed in the file,\n"
            "fill the rows, and run this again. Nothing is scored from an empty panel.")

    from quadcond import Condition
    from quadcond.models.predict import Predictor

    pred = Predictor.load(a.model, a.atlas if Path(a.atlas).exists() else None)
    frozen = meta.get("model_artifact_sha256_at_freeze", "")
    if frozen and pred.model_sha256 and frozen != pred.model_sha256:
        raise SystemExit(
            f"{path.name} was declared held out against model {frozen[:12]}, but the loaded\n"
            f"artifact is {pred.model_sha256[:12]}. Either score it with the declared model or\n"
            "update the panel file deliberately; a silent substitution makes the result\n"
            "unauditable.")
    if not meta.get("measurements_published") and not a.allow_undated:
        raise SystemExit(f"{path.name}: no 'measurements_published' date in the provenance "
                         "block, so it cannot be shown to postdate the model freeze. "
                         "Add it, or pass --allow-undated and say so in the write-up.")

    out_rows = []
    for r in rows:
        cond_map = {field: num(r.get(col)) for field, col in CONDITION_COLS.items()}
        cond = Condition.from_mapping({k: v for k, v in cond_map.items() if v is not None})
        res = pred.predict(r["sequence"], cond)[0]
        rec = {"sequence_id": r["sequence_id"], "sequence": r["sequence"],
               "condition": res["condition"],
               "condition_imputed": ";".join(res["condition_imputed_fields"])}
        for obs in OBSERVED:
            if obs in r:
                rec[obs] = r[obs]
        for head, p in res["predictions"].items():
            if p.get("refused"):
                rec[f"{head}__refused"] = p["refusal_reason"]
                continue
            for key in ("probability", "value", "label"):
                if key in p:
                    rec[f"{head}__{key}"] = p[key]
            rec[f"{head}__domain"] = p["applicability"].get("status", "")
        out_rows.append(rec)

    ROOT.joinpath("results").mkdir(exist_ok=True)
    stem = path.stem
    cols = sorted({k for r in out_rows for k in r},
                  key=lambda c: (c not in ("sequence_id", "sequence", "condition"), c))
    with (ROOT / "results" / f"{stem}.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)

    summary = {"panel": stem, "n": len(out_rows), "scored_on": date.today().isoformat(),
               "model_path": pred.model_path, "model_sha256": pred.model_sha256,
               "model_version": getattr(pred.model, "version", None),
               "provenance": meta, "metrics": metrics(out_rows)}
    (ROOT / "results" / f"{stem}_summary.json").write_text(json.dumps(summary, indent=2),
                                                           encoding="utf-8")
    print(json.dumps(summary["metrics"], indent=2))
    print(f"\nwritten: results/{stem}.csv and results/{stem}_summary.json")


def metrics(rows):
    """Only the metrics the panel's own observations support."""
    import numpy as np
    out = {}

    def pairs(obs_col, pred_col):
        xs, ys = [], []
        for r in rows:
            o, p = num(r.get(obs_col)), num(r.get(pred_col))
            if o is not None and p is not None:
                xs.append(o); ys.append(p)
        return np.array(xs), np.array(ys)

    y, p = pairs("observed_folded", "g4_fold__probability")
    if len(y) >= 10 and len(set(y.tolist())) == 2:
        out["g4_fold"] = {"n": int(len(y)), "auroc": _auroc(y, p)}
    for obs, head, label in (("observed_tm_c", "g4_tm__value", "g4_tm"),
                             ("observed_pht", "im_pht__value", "im_pht")):
        y, p = pairs(obs, head)
        if len(y) >= 5:
            from scipy import stats
            out[label] = {"n": int(len(y)),
                          "spearman": float(stats.spearmanr(y, p).statistic),
                          "pearson": float(stats.pearsonr(y, p).statistic),
                          "mae": float(np.mean(np.abs(y - p))),
                          "bias": float(np.mean(p - y))}
    topo = [(r.get("observed_topology"), r.get("g4_topology__label")) for r in rows
            if r.get("observed_topology") and r.get("g4_topology__label")]
    if len(topo) >= 10:
        agree = sum(1 for o, q in topo if str(o).strip().lower() == str(q).strip().lower())
        out["g4_topology"] = {"n": len(topo), "exact_agreement": agree / len(topo)}
    refused = sum(1 for r in rows if any(k.endswith("__refused") for k in r))
    out["rows_with_a_refusal"] = refused
    return out


def _auroc(y, p):
    import numpy as np
    order = np.argsort(p)
    ranks = np.empty(len(p), float)
    ranks[order] = np.arange(1, len(p) + 1)
    pos, neg = y == 1, y == 0
    n1, n0 = pos.sum(), neg.sum()
    return float((ranks[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


if __name__ == "__main__":
    main()
