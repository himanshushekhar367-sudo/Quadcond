#!/usr/bin/env python3
"""v0.4.3 corrections applied to the trained artifacts. No refitting.

Three things a reviewer found that are wrong in the *metadata*, not in the
fitted functions, and can therefore be corrected without another training run.
Every change is a rename or a string; the estimators are untouched, which is why
the ablation and the reported metrics stay valid.

**1. `locus_state` is named for something it does not measure.**

Zanin et al. ran iMab and BG4 CUT&Tag as *parallel reactions on separate
aliquots* of the same cell population, and reported peaks that overlap between
the two experiments. That supports population-level peak co-localisation. It
does not support co-occupancy: nothing in the design observes both antibodies
bound to the same DNA molecule, and a name like "joint state" invites exactly
that reading. Renamed to ``locus_peak_overlap_state``, with the claim rewritten
to say what the label is.

The head also carries limits the old claim did not state:

- every training window is **201 nt**, so the applicability table's length range
  is 201-201 and any oligo-length query is out of domain by construction;
- balanced accuracy is **0.450** against a 0.250 floor;
- **iM-only recall is 0.156** -- the class the tool is most often asked about is
  the one it recovers worst;
- the "neither" class is a composition-matched shuffle, not an observed
  unoccupied locus, so three of the four classes are observations and the fourth
  is generated.

**2. The condition heads' allowlist needs to read as a refusal.**

Already enforced in ``_domain_check``; this records the matching prose in the
head's own claim, so the artifact says it too.

**3. Two condition fields were never in the applicability tables.**

``crowder_pct`` and ``strand_conc`` are on every ``Condition`` and were absent
from every head's applicability table, so the domain check had nothing to
compare a query against and skipped them silently. Their ranges are recomputed
here from the atlas rows each head was trained on -- a measurement of the
training data, not a refit.

    python scripts/14_v043_corrections.py
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import joblib

from quadcond import __version__
from quadcond.models.base import MultiTaskModel

OLD, NEW = "locus_state", "locus_peak_overlap_state"

PEAK_OVERLAP_CLAIM = (
    "GENOMIC PEAK-OVERLAP PROXY FOR A 201-NT WINDOW. NOT CO-OCCUPANCY AND NOT "
    "COMPETITION. Calibrated posterior over how a 201-nt genomic window was "
    "classified by two CUT&Tag experiments run as PARALLEL REACTIONS on separate "
    "aliquots of the same HEK293T population: reproducible peaks from both "
    "antibodies, from BG4 only, from iMab only, or from neither (a "
    "composition-matched shuffle of the same windows -- generated, not an "
    "observed unoccupied locus). Because the two antibodies never touched the "
    "same aliquot, a 'both' call is population-level co-localisation of peaks; "
    "nothing here observes two structures on one DNA molecule, and nothing here "
    "is a free energy. Balanced accuracy 0.450 against a 0.250 four-class floor, "
    "iM-only recall 0.156 -- the class most users are asking about is the one it "
    "recovers worst. Trained exclusively on 201-nt windows, so a query of "
    "oligonucleotide length is outside the applicability domain by construction. "
    "Read it as an exploratory genomic-resemblance score."
)

CONDITION_CLAIM_SUFFIX = (
    " REFUSES ANY SEQUENCE OUTSIDE ITS TWO-CONSTRUCT ALLOWLIST: no value is "
    "returned, because this head gives essentially the same answer whatever "
    "sequence it is given, so there is no prediction to report rather than an "
    "uncertain one."
)


def condition_ranges(db: Path, sources) -> dict | None:
    """min/max/percentiles for the two missing fields, over a head's own rows."""
    import sqlite3

    import numpy as np

    if not db.exists():
        return None
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        if sources:
            q = ("SELECT crowder_pct, strand_conc FROM records WHERE source IN "
                 f"({','.join('?' * len(sources))})")
            rows = con.execute(q, list(sources)).fetchall()
        else:
            rows = con.execute("SELECT crowder_pct, strand_conc FROM records").fetchall()
    finally:
        con.close()
    if not rows:
        return None

    def rng(vals):
        a = np.asarray([v for v in vals if v is not None], dtype=float)
        if not a.size:
            return None
        return {"min": float(a.min()), "max": float(a.max()),
                "p05": float(np.percentile(a, 5)), "p95": float(np.percentile(a, 95)),
                "n_unique": int(len(np.unique(np.round(a, 3))))}

    out = {"crowder_pct": rng([r[0] for r in rows]),
           "strand_conc": rng([r[1] for r in rows])}
    return {k: v for k, v in out.items() if v}


def correct(path: Path, *, backup: bool = True, db: Path = Path("data/atlas.db")) -> list[str]:
    model = MultiTaskModel.load(path, verify=False)
    notes: list[str] = []

    for name, head in model.heads.items():
        missing = [f for f in ("crowder_pct", "strand_conc") if f not in head.applicability]
        if not missing:
            continue
        got = condition_ranges(db, head.training_meta.get("sources"))
        if got:
            head.applicability.update({k: v for k, v in got.items() if k in missing})
            notes.append(f"{name}: applicability gained {', '.join(missing)}")

    # --- the canonical artifact, not the presentation layer ------------------
    #
    # Two corrections had been applied only where a human reads them. The report
    # and the README were regenerated with the right numbers while the model's
    # own sidecar -- the file the model card is built from and the API serves --
    # still carried the old ones. A number is corrected when the artifact
    # carries it, not when the page does.

    # 1. The atlas snapshot's per-source `n_conditions` was computed with a key
    #    that omitted Mg2+, giving 255 for g4stab where every other surface now
    #    says 261.
    if db.exists():
        from quadcond.atlas import Atlas

        atlas = Atlas(db)
        try:
            fresh = atlas.summary().to_dict(orient="records")
        finally:
            atlas.close()
        old_rows = (model.atlas_snapshot or {}).get("summary") or []
        old_by_source = {r.get("source"): r for r in old_rows}
        changed = [r["source"] for r in fresh
                   if old_by_source.get(r["source"], {}).get("n_conditions")
                   not in (None, r["n_conditions"])]
        if changed:
            model.atlas_snapshot = {**(model.atlas_snapshot or {}), "summary": fresh}
            notes.append(f"atlas snapshot: n_conditions refreshed for "
                         f"{len(changed)} source(s) (Mg2+ was missing from the key)")

    # 2. `_domain_check` skips the query temperature for a Tm-target head -- it
    #    predicts a temperature, it does not consume one -- but the field stayed
    #    in the head's applicability table, which /info serves verbatim. The
    #    served contract described a gate the code does not apply.
    for name, head in model.heads.items():
        ap = head.applicability or {}
        if head.target == "tm" and "temperature" in ap and "temperature_is_enforced" not in ap:
            ap["temperature_is_enforced"] = False
            ap["temperature_note"] = (
                "NOT a domain gate for this head: it PREDICTS a melting "
                "temperature and does not take one as an input. Every training "
                "row records temperature=25 as a metadata default on a Tm "
                "measurement, so this range describes the records rather than a "
                "limit on what you may ask.")
            notes.append(f"{name}: temperature marked non-enforced in applicability")

    if OLD in model.heads:
        head = model.heads.pop(OLD)
        head.name = NEW
        head.training_meta["claim"] = PEAK_OVERLAP_CLAIM
        head.training_meta["renamed_from"] = OLD
        head.training_meta["window_nt"] = 201
        model.heads[NEW] = head
        notes.append(f"{OLD} -> {NEW}, claim rewritten")

    for name in ("im_tm_condition", "im_pht_condition"):
        h = model.heads.get(name)
        if h is None:
            continue
        claim = h.training_meta.get("claim", "")
        if "REFUSES ANY SEQUENCE" not in claim:
            h.training_meta["claim"] = claim.rstrip() + CONDITION_CLAIM_SUFFIX
            notes.append(f"{name}: refusal recorded in the claim")

    if model.version != __version__:
        notes.append(f"version {model.version} -> {__version__}")
        model.version = __version__

    if notes:
        if backup:
            keep = path.with_suffix(path.suffix + ".pre-v0.4.3")
            if not keep.exists():
                shutil.copy2(path, keep)
        # `model.save()`, never a bare joblib.dump: save() rewrites the sidecar
        # .json with the checksum of what it actually wrote. A dump leaves the
        # old checksum in place and the loader then refuses the file -- which is
        # the integrity check behaving correctly against a stale record, and is
        # exactly how this script failed the first time it ran.
        model.save(path)
    return notes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("artifacts", nargs="*",
                    default=["artifacts/quadcond_model.joblib",
                             "artifacts/quadcond_model_seqonly.joblib"])
    a = ap.parse_args()

    for raw in a.artifacts:
        path = Path(raw)
        if not path.exists():
            print(f"  [skip] {path} not present")
            continue
        notes = correct(path)
        print(f"{path} -> {__version__}")
        for n in notes:
            print(f"    {n}")
        if not notes:
            print("    already corrected")

    print("\nNothing was refitted. Every change is a name or a string; the "
          "estimators, the folds and the reported metrics are the same objects, "
          "which is what makes it safe to do this without retraining.")


if __name__ == "__main__":
    main()
