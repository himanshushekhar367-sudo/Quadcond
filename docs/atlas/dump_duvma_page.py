#!/usr/bin/env python3
"""Show what pypdf actually extracts from the 5DUVMA table pages.

Table S9 parses to zero rows on a rebuild (docs/REBUILD_THE_ATLAS.md), costing
115 melting temperatures and with them any faithful retraining of
`im_tm_condition`. S9 is the only table the adapter treats as `multirow`, and its
row pattern is strict:

    ^(\\d+\\.\\d+)\\s+iM\\s+(.*)$

a decimal pH, whitespace, the literal "iM", then ten values. If the installed
pypdf lays that page out differently -- the label on its own line, a non-breaking
space, columns joined -- every line fails the match and the table silently yields
nothing. This prints the extracted lines with repr() so the whitespace is visible,
and reports which ones the pattern matches.

    python docs/atlas/dump_duvma_page.py                 # S9
    python docs/atlas/dump_duvma_page.py --table S10     # a table that works, to compare
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

MULTIROW = re.compile(r"^(\d+\.\d+)\s+iM\s+(.*)$")
PLAIN = re.compile(r"^(\d+(?:\.\d+)?)\s+(.*)$")
HEADER_IS = re.compile(r"(\d+(?:\.\d+)?)\s*M\b")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default=str(REPO / "data" / "external" / "downloads"
                                        / "gkag110_supplemental_file.pdf"))
    ap.add_argument("--table", default="S9")
    ap.add_argument("--lines", type=int, default=30)
    a = ap.parse_args()

    from quadcond.atlas.ingest.duvma import TABLES
    spec = next((t for t in TABLES if t["table"] == a.table), None)
    if spec is None:
        raise SystemExit(f"no table {a.table}; known: {[t['table'] for t in TABLES]}")

    import pypdf
    print(f"pypdf {pypdf.__version__}")
    reader = pypdf.PdfReader(a.pdf)
    print(f"{len(reader.pages)} pages; table {a.table} is expected on index {spec['page']}\n")

    text = reader.pages[spec["page"]].extract_text() or ""
    lines = text.split("\n")
    pat = MULTIROW if spec.get("multirow") else PLAIN
    header = HEADER_IS.findall(text)
    print(f"ionic-strength header tokens found on the page: {header[:12]}")
    print(f"pattern: {pat.pattern}\n")
    hits = 0
    for i, raw in enumerate(lines[: a.lines]):
        line = raw.strip()
        m = pat.match(line)
        hits += bool(m)
        print(f"{i:3} {'MATCH' if m else '     '} {line[:110]!r}")
    total = sum(bool(pat.match(l.strip())) for l in lines)
    print(f"\n{total} of {len(lines)} lines match (first {a.lines} shown: {hits})")
    if not total:
        print("Nothing matches. Compare the repr above against the pattern -- the usual "
              "causes are the label separated from the pH, a non-breaking space, or the "
              "pH and first value joined without a space.")


if __name__ == "__main__":
    main()
