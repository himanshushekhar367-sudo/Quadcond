#!/usr/bin/env python3
"""Prove the panel harness runs end to end, without leaving fake data behind.

The panels in panels/ are empty until someone extracts the sequences from the
sources they name, so the scoring path is untested until that happens -- and
discovering it is broken *after* doing the extraction is the wrong order. This
builds a throwaway panel of well-known promoter sequences in a temporary
directory, scores it, checks the output has the columns the write-up needs, and
deletes it.

The observations in here are invented. That is why this lives in a temp
directory and never in panels/: a file of made-up measurements sitting beside
real ones is an accident waiting to be cited.

    ..\..\.venv\Scripts\python.exe selftest.py        (Windows, the model environment)
    python selftest.py                                  (any environment that can import quadcond)

It scores with whichever interpreter started it, so running it with the wrong
one fails the same way a real panel would -- which is the point.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Sequences are real and public; the observed_* values are placeholders that
# exist only to exercise the metric code paths.
SEQS = [
    ("MYC_Pu27", "TGGGGAGGGTGGGGAGGGTGGGGAAGG", 1, "parallel"),
    ("VEGFA_Pu22", "GGGGCGGGCCGGGGGCGGGG", 1, "parallel"),
    ("hTERT_hT21", "AGGGGAGGGGCTGGGAGGGCC", 1, "parallel"),
    ("scrambled_1", "ATCATGCATGACTAGCATCGA", 0, ""),
    ("scrambled_2", "CATAGCATCAGTACTAGCATC", 0, ""),
    ("scrambled_3", "TACGATCAGTACGATCAGTAC", 0, ""),
]
HEADER = ("sequence_id\tsequence\tmolecule\tbuffer_k_mm\tbuffer_na_mm\tbuffer_mg_mm\tph\t"
          "temperature_c\tstrand_conc_um\tobserved_folded\tobserved_topology\tobserved_tm_c\t"
          "source_table\tnotes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT.parent.parent / "artifacts" / "quadcond_model.joblib"))
    a = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        panel = Path(td) / "selftest_synthetic.tsv"
        lines = [
            "# panel: selftest_synthetic",
            "# title: THROWAWAY -- invented observations, for exercising the harness only",
            "# measurements_published: 1970-01-01",
            "# why_this_panel: not a panel. See selftest.py.",
            HEADER,
        ]
        for sid, seq, folded, topo in SEQS:
            lines.append(f"{sid}\t{seq}\tDNA\t5\t132\t1.5\t7.4\t25\t20\t{folded}\t{topo}\t\t-\t-")
        panel.write_text("\n".join(lines) + "\n", encoding="utf-8")

        out = subprocess.run([sys.executable, str(ROOT / "run_panel.py"), str(panel),
                              "--model", a.model],
                             capture_output=True, text=True, cwd=str(ROOT))
        if out.returncode:
            print(out.stdout); print(out.stderr, file=sys.stderr)
            raise SystemExit("harness failed on the synthetic panel; fix that before extracting "
                             "real sequences")
        csv_path = ROOT / "results" / "selftest_synthetic.csv"
        json_path = ROOT / "results" / "selftest_synthetic_summary.json"
        # Cleanup goes in a finally. It used to sit after the assertions, so a
        # failing check left a file of invented measurements in results/ beside
        # the real ones -- which is the exact accident this script exists to
        # avoid, and it happened: two selftest_synthetic files survived from a
        # run on 22 September.
        try:
            rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
            assert len(rows) == len(SEQS), f"expected {len(SEQS)} rows, got {len(rows)}"
            cols = set(rows[0])
            need = {"sequence_id", "sequence", "condition", "condition_imputed"}
            missing = need - cols
            assert not missing, f"output is missing {missing}"
            assert any(c.startswith("g4_fold__") for c in cols), "no g4_fold output at all"
            # Conditions must reach the model: nothing the panel supplied may be imputed.
            for r in rows:
                assert not r["condition_imputed"], (
                    f"{r['sequence_id']}: the harness imputed {r['condition_imputed']} even "
                    "though the panel supplied it -- the buffer is not reaching the model")
        finally:
            for leftover in (csv_path, json_path):
                leftover.unlink(missing_ok=True)
    print(f"harness OK: {len(SEQS)} sequences scored, conditions passed through, "
          "temporary outputs removed")


if __name__ == "__main__":
    main()
