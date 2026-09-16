#!/usr/bin/env python3
"""Generate docs/MODEL_CARD.md from the trained model -- never by hand.

A model card that is written by hand drifts from the model it describes. This
one is regenerated from ``MultiTaskModel.summary()`` every time, so the numbers
in it are the numbers the model actually reports.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond import claims
from quadcond.models.base import MultiTaskModel

KEY_METRIC = {"binary": "auroc", "multiclass": "balanced_accuracy", "regression": "r2"}


def fmt(v, nd=4):
    if v is None:
        return "—"
    if isinstance(v, float):
        return "—" if v != v else f"{v:.{nd}f}"
    return str(v)


def _version_of_package() -> str | None:
    """The package version, for the header's version note."""
    try:
        from quadcond import __version__
        return __version__
    except Exception:                                       # noqa: BLE001
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--ablation", default="artifacts/quadcond_model_seqonly.joblib")
    ap.add_argument("--out", default="docs/MODEL_CARD.md")
    a = ap.parse_args()

    m = MultiTaskModel.load(a.model)
    seq_only = MultiTaskModel.load(a.ablation) if Path(a.ablation).exists() else None
    when = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    L: list[str] = []
    L.append("# QuadCond model card\n")
    # Emitted by the generator, not pasted into the output afterwards: this
    # file says "do not edit by hand", so a note added by hand is a note that
    # disappears the next time anyone regenerates it.
    _app = _version_of_package()
    _model = getattr(m, "version", None)
    if _model and _app and _model != _app:
        L.append(f"> **VERSION NOTE** -- package {_app} ships the model "
                 f"artifact stamped **{_model}**. The application and the "
                 f"model are versioned separately; see `RELEASE_{_app}.md` "
                 f"for what that rests on.\n")
    L.append(f"*Generated {when} by `scripts/03_model_card.py` from "
             f"`{a.model}`. Do not edit by hand.*\n")

    L.append("## Summary\n")
    L.append(f"**{claims.headline(m.heads)}**\n")
    L.append("`claim basis` is the field to read. It is not the same question as "
             "\"was anything measured\": a CUT&Tag peak is an experimental "
             "observation whose target is antibody occupancy at a locus, so it "
             "carries `has_experimental_observation = true` and "
             "`biophysically_grounded = false`. Only the latter licenses a "
             "folding or stability claim.\n")
    L.append("| head | kind | task | claim basis | experimental? | biophysical? | rows | key metric |")
    L.append("|---|---|---|---|---|---|---|---|")
    for name, h in m.heads.items():
        tm_ = h.training_meta
        sem = claims.semantics(name, tm_, h.task)
        key = KEY_METRIC[h.task]
        L.append(
            f"| `{name}` | {h.kind} | {h.task} | {sem['target_semantics']} | "
            f"{'yes' if sem['has_experimental_observation'] else 'no'} | "
            f"{'**yes**' if sem['biophysically_grounded'] else 'no'} | "
            f"{tm_.get('n_rows', 0):,} | {key} = {fmt(h.metrics.get(key))} |"
        )
    L.append("")

    by = {}
    for name, h in m.heads.items():
        by.setdefault(claims.target_semantics(name, h.training_meta), []).append(name)
    grounded = by.get(claims.BIOPHYSICAL, [])
    ungrounded = [n for k, v in by.items() if k != claims.BIOPHYSICAL for n in v]
    if by.get(claims.GENOMIC_PROXY):
        L.append("> **Genomic proxy, not biophysics:** "
                 + ", ".join(f"`{n}`" for n in by[claims.GENOMIC_PROXY]) +
                 ". The positive label is antibody occupancy at a genomic locus "
                 "(iMab / BG4 CUT&Tag), which is an experimental observation of a "
                 "different quantity from folding in a defined buffer -- and one whose "
                 "interpretation is contested, since iMab specificity and the "
                 "accessible-chromatin contribution to targeted CUT&Tag signal are both "
                 "open questions. No buffer was measured for these rows; the attached "
                 "condition is nominal with every field flagged imputed. Read these "
                 "heads as resemblance to peak-enriched motifs, never as P(folds).\n")
    aux = by.get(claims.DERIVED, []) + by.get(claims.PREDICTED, [])
    if aux:
        L.append("> **Derived / prior-model auxiliary:** " + ", ".join(f"`{n}`" for n in aux) +
                 ". No measurement of any kind stands behind these labels. They exist to "
                 "exercise the pipeline and to supply a condition-response prior; their "
                 "numbers are not evidence about DNA.\n")

    L.append("## Evaluation protocol\n")
    L.append("- **Grouped cross-validation.** Sequences are clustered on 4-mer cosine "
             "similarity and folds are built over clusters, so near-duplicates cannot "
             "straddle a split. Group statistics are reported per head below.")
    L.append("- **Out-of-fold metrics only.** Nothing is reported on data a model saw.")
    L.append("- **Cross-fitted calibration.** The calibrator is itself cross-fitted over "
             "groups before ECE/Brier are computed, so calibration is not scored on the "
             "calibrator's own training data. The deployed calibrator is separately fitted "
             "on all out-of-fold predictions.")
    L.append("- **OOF-residual interval.** Regression heads report a split-conformal "
             "half-width from out-of-fold residuals, with the realised coverage stated.\n")

    if seq_only is not None:
        L.append("## Condition ablation\n")
        L.append("Same data, same protocol, condition and interaction features removed.\n")
        L.append("| head | metric | sequence + conditions | sequence only | delta |")
        L.append("|---|---|---|---|---|")
        for name, h in m.heads.items():
            if name not in seq_only.heads:
                continue
            key = KEY_METRIC[h.task]
            va, vb = h.metrics.get(key), seq_only.heads[name].metrics.get(key)
            if va is None or vb is None:
                continue
            L.append(f"| `{name}` | {key} | {fmt(va)} | {fmt(vb)} | {va - vb:+.4f} |")
        L.append("")
        L.append("A delta near zero is a finding about the *data*, not the method: it means "
                 "that head's training rows do not vary the buffer enough for conditioning "
                 "to be learnable. Ingest condition-resolved measurements and re-run.\n")

    L.append("## Heads in detail\n")
    for name, h in m.heads.items():
        tm_, mt = h.training_meta, h.metrics
        L.append(f"### `{name}`\n")
        L.append(f"**Claim.** {tm_.get('claim', '(none recorded)')}\n")
        L.append(f"- kind: {h.kind}; task: {h.task}; target column: `{h.target}`")
        L.append(f"- training rows: {tm_.get('n_rows', 0):,}; features: {tm_.get('n_features')}")
        L.append(f"- evidence tiers: {', '.join(tm_.get('tiers', []))}")
        L.append(f"- sources: {', '.join(f'`{s}`' for s in tm_.get('sources', []))}")
        g = tm_.get("grouping", {})
        L.append(f"- sequence groups: {g.get('n_groups', 0):,} from {g.get('n_items', 0):,} "
                 f"sequences (largest cluster {g.get('largest_group', 0)}, "
                 f"{g.get('singleton_groups', 0):,} singletons)")
        L.append(f"- CV: {tm_.get('n_folds')} grouped folds x {tm_.get('n_seeds')} seeds\n")

        L.append("| metric | value |")
        L.append("|---|---|")
        for k, v in mt.items():
            if k in {"reliability", "per_class", "calibration_note", "target_range"}:
                continue
            L.append(f"| {k} | {fmt(v)} |")
        L.append("")
        if mt.get("per_class"):
            L.append("| class | n | recall |")
            L.append("|---|---|---|")
            for c, d in mt["per_class"].items():
                L.append(f"| {c} | {d['n']} | {fmt(d.get('recall'), 3)} |")
            L.append("")

        L.append("**Applicability domain** — the model saw only this much variation, so "
                 "queries outside it are flagged at prediction time:\n")
        L.append("| variable | range | distinct values seen | enforced |")
        L.append("|---|---|---|---|")
        # The applicability dict holds range objects AND scalar annotations
        # (`sequence_allowlist`, `temperature_is_enforced`, `*_note`). Iterating
        # it as though every value were a range is how a new annotation broke
        # this generator; the ranges are the entries that look like ranges.
        ranges = {k: v for k, v in (h.applicability or {}).items()
                  if isinstance(v, dict) and "min" in v and "n_unique" in v}
        for k, v in ranges.items():
            # A field the domain check skips is documented as not enforced
            # rather than listed as though it gated a query. A Tm head predicts
            # a temperature; it does not consume one.
            enforced = (h.applicability or {}).get(f"{k}_is_enforced", True)
            note = "yes" if enforced else "**no** — see note below"
            L.append(f"| {k} | {v['min']:g} – {v['max']:g} | {v['n_unique']} | {note} |")
        L.append("")
        for k in ranges:
            n = (h.applicability or {}).get(f"{k}_note")
            if n:
                L.append(f"> `{k}`: {n}\n")
        allow = (h.applicability or {}).get("sequence_allowlist")
        if allow:
            L.append(f"> This head answers for {len(allow)} specific construct(s) and "
                     f"**returns no value for any other sequence**.\n")
        one_value = [k for k, v in ranges.items() if v["n_unique"] <= 1]
        if one_value:
            L.append(f"> This head saw a single value of: {', '.join(f'`{k}`' for k in one_value)}. "
                     f"It cannot resolve those variables, and predictions that vary them are "
                     f"extrapolation, not inference.\n")

    L.append("## Intended use\n")
    L.append("- Prioritising candidate G4/i-motif elements for experimental follow-up, "
             "*under a stated buffer*, with a calibrated probability rather than a rank.")
    L.append("- Asking how a prediction changes across a cation or pH gradient "
             "(`quadcond sweep`), including whether the model is entitled to answer.")
    L.append("- Retrieving the nearest real measurements to any query "
             "(`quadcond atlas neighbours`), so a score can be checked against evidence.\n")

    L.append("## Out of scope\n")
    L.append("- **Structure prediction.** No 3D coordinates, no ligand docking.")
    L.append("- **RNA.** The bundled data is DNA; RNA sequences are accepted (U→T) but no "
             "RNA-specific head has been trained.")
    L.append("- **In-cell occupancy.** Nothing here models chromatin, transcription, or "
             "protein binding. A high in-vitro folding probability is not a claim about a "
             "structure existing in a nucleus.")
    L.append("- **Clinical or diagnostic use.**\n")

    L.append("## Known biases\n")
    L.append("- The topology data is predominantly human and predominantly 20–25 nt; "
             "long-looped, bulged and multimeric species are under-represented.")
    L.append("- The topology classes are imbalanced (parallel is the plurality), so "
             "balanced accuracy and per-class recall are reported alongside accuracy.")
    L.append("- Distillation heads inherit every bias of the model they distil, including "
             "its training-buffer coverage.\n")

    L.append("## Reproducing\n")
    L.append("```bash\npython scripts/01_build_atlas.py\npython scripts/02_train.py\n"
             "python scripts/03_model_card.py\nquadcond report\n```\n")

    L.append("<details><summary>Raw model summary (JSON)</summary>\n")
    L.append("```json")
    L.append(json.dumps(m.summary(), indent=2, default=str))
    L.append("```\n</details>")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L), encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
