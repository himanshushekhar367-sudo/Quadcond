#!/usr/bin/env python3
"""v0.4.2: four integrity fixes, none of which needs a retrain.

All four were found by review rather than by a failing test, and all four share
a shape: the model was right and the *description* of it was not, in a direction
that flattered the tool.

**1. Heads with two training sequences claimed to be in domain for any sequence.**
``im_tm_condition`` and ``im_pht_condition`` learned how Tel21C and C9 respond to
salt, pH and temperature. Asked about a different C-tract sequence at an in-range
buffer, the applicability check returned ``in_domain=True`` and a confident
44.82 C -- against 44.73 C for Tel21C itself. The value was not merely
unlicensed, it was uninformative: a two-sequence panel gives the model almost no
way to depend on sequence at all, so the head returns nearly the same number for
anything. Every head trained on fewer than ``--allowlist-threshold`` distinct
sequences now carries the list of those sequences, and refuses anything else.

**2. The cation-identity check was in-sample.** The K+/Na+ comparison scored
411 matched pairs with the final model, which was trained on those same records.
It is a mechanistic consistency check and was reported as validation. The
corrected version is in ``scripts/10_cation_identity_check.py``: grouped
hold-out predictions plus a sequence-level bootstrap interval.

**3. "Split-conformal" survived in the prose.** The frozen artifact says
``oof_residual_halfwidth`` and ``in_sample_coverage_of_residual_pool``, which is
accurate. The README and the model-card generator still said conformal, which
promises a finite-sample guarantee the computation never had.

**4. The test count was typed, not counted.** The report said 52; the archive
contains 56 test functions. A number a reader would reasonably take as measured.

    python scripts/12_v042_integrity.py --version 0.4.2
"""
from __future__ import annotations

import argparse
import re
import shutil
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond.models.base import MultiTaskModel
from quadcond.motifs import clean

STALE_CONFORMAL = [
    (re.compile(r"split-conformal (half-width|interval|coverage)", re.I),
     "OOF-residual \\1"),
    (re.compile(r"\bconformal (interval|half-width)s?\b", re.I),
     "OOF-residual \\1"),
    (re.compile(r"finite-sample coverage", re.I),
     "in-sample coverage of the residual pool"),
]


def training_sequences(db: Path, sources: tuple[str, ...]) -> list[str]:
    if not sources:
        return []
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        q = ",".join("?" for _ in sources)
        rows = con.execute(
            f"SELECT DISTINCT sequence FROM records WHERE source IN ({q})",
            tuple(sources)).fetchall()
    finally:
        con.close()
    return sorted({clean(r[0]) for r in rows})


def add_allowlists(model: MultiTaskModel, db: Path, threshold: int) -> list[str]:
    notes = []
    for name, head in model.heads.items():
        srcs = tuple(head.training_meta.get("sources") or ())
        seqs = training_sequences(db, srcs)
        if not seqs or len(seqs) > threshold:
            head.applicability.pop("sequence_allowlist", None)
            continue
        head.applicability["sequence_allowlist"] = seqs
        head.applicability["sequence_allowlist_note"] = (
            f"This head was trained on {len(seqs)} distinct sequences. It learned "
            "how those constructs respond to conditions, not how sequence affects "
            "the response. Any other sequence is out of domain and is refused."
        )
        notes.append(f"{name}: allowlist of {len(seqs)} sequence(s) added")
    return notes


def fix_prose(paths: list[Path]) -> list[str]:
    notes = []
    for p in paths:
        if not p.exists():
            continue
        text = original = p.read_text()
        for pattern, repl in STALE_CONFORMAL:
            text = pattern.sub(repl, text)
        if text != original:
            p.write_text(text)
            notes.append(f"{p}: stale conformal wording corrected")
    return notes


def count_tests(root: Path) -> int:
    """Count test functions rather than trusting a number in a document."""
    n = 0
    for p in (root / "tests").glob("test_*.py"):
        n += len(re.findall(r"^def test_", p.read_text(), re.M))
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="0.4.2")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--allowlist-threshold", type=int, default=20,
                    help="heads with at most this many distinct training "
                         "sequences get a sequence allowlist")
    ap.add_argument("--artifacts", nargs="*", default=[
        "artifacts/quadcond_model.joblib",
        "artifacts/quadcond_model_seqonly.joblib"])
    a = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    db = Path(a.db)

    print(f"tests present in this tree: {count_tests(root)} "
          "(counted, not typed)\n")

    for path in a.artifacts:
        p = Path(path)
        if not p.exists():
            print(f"  ! {p} not found; skipping")
            continue
        keep = p.with_suffix(p.suffix + ".v0.4.1")
        if not keep.exists():
            shutil.copy2(p, keep)
        model = MultiTaskModel.load(p, verify=True)
        notes = add_allowlists(model, db, a.allowlist_threshold)
        model.version = a.version
        model.save(p)
        print(f"{p} -> v{model.version}")
        for n in notes:
            print(f"    {n}")

    for n in fix_prose([root / "README.md",
                        root / "docs" / "QUICKSTART.md",
                        root / "scripts" / "03_model_card.py",
                        root / "scripts" / "report_templates" / "body.html"]):
        print(f"    {n}")

    print("\nNothing was refitted. The estimators are the v0.4.1 objects; what "
          "changed is what they are willing to answer.")


if __name__ == "__main__":
    main()
