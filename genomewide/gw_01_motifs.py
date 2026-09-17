#!/usr/bin/env python3
"""Step 1: every canonical G4 and i-motif in hg38, plus a matched control window for each.

G4:  G{3+} N{1-7} x4 on either strand (reported in G-rich orientation).
iM:  C{3+} N{1-12} x4 on the forward strand, and the same on the reverse strand
     (i.e. the C-rich partner of a G-rich motif is a separate iM record).
Control: a window of identical length 2-20 kb away, same chromosome, no N, no
canonical G4/iM on either strand within +/-50 nt, GC within 0.05 of the motif.

Output: $GW_OUT/motifs/<chrom>.tsv.gz
  id chrom start end strand kind sequence gc ctrl_start ctrl_end ctrl_gc
  (0-based half-open coordinates on the + strand; `sequence` in motif orientation)
"""
from __future__ import annotations

import argparse
import gzip
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gw_common import CHROMS, out_dir, read_chrom, revcomp  # noqa: E402

G4 = re.compile(r"(?=(G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}))")
IM = re.compile(r"(?=(C{3,}[ACGT]{1,12}C{3,}[ACGT]{1,12}C{3,}[ACGT]{1,12}C{3,}))")
ANY = re.compile(r"G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}|"
                 r"C{3,}[ACGT]{1,12}C{3,}[ACGT]{1,12}C{3,}[ACGT]{1,12}C{3,}")
MAXLEN = 60   # longer hits are split by the regex into their first canonical core


def non_overlapping(pat, s):
    last = -1
    for m in pat.finditer(s):
        a = m.start()
        b = a + len(m.group(1))
        if a < last or b - a > MAXLEN:
            continue
        last = b
        yield a, b


def gc(s):
    return (s.count("G") + s.count("C")) / max(1, len(s))


def find_control(seq, a, b, g, rng, tries=40):
    L = b - a
    for _ in range(tries):
        off = rng.randint(2000, 20000) * rng.choice((-1, 1))
        c = a + off
        if c - 50 < 0 or c + L + 50 > len(seq):
            continue
        w = seq[c:c + L]
        if "N" in w or abs(gc(w) - g) > 0.05:
            continue
        if ANY.search(seq[c - 50:c + L + 50]):
            continue
        return c, c + L, gc(w)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="indexed GRCh38 FASTA (.fai beside it)")
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    od = out_dir() / "motifs"
    od.mkdir(exist_ok=True)
    for chrom in a.chroms:
        dest = od / f"{chrom}.tsv.gz"
        if dest.exists():
            print(f"{chrom}: exists, skipped")
            continue
        seq = read_chrom(a.fasta, chrom)
        rng = random.Random(f"{a.seed}-{chrom}")
        rows = []
        for kind, pat in (("G4", G4), ("iM", IM)):
            for strand in "+-":
                s = seq if strand == "+" else revcomp(seq)
                n = len(s)
                for x, y in non_overlapping(pat, s):
                    a0, b0 = (x, y) if strand == "+" else (n - y, n - x)
                    mot = s[x:y]
                    g = gc(mot)
                    ctrl = find_control(seq, a0, b0, g, rng)
                    rows.append((chrom, a0, b0, strand, kind, mot, round(g, 3),
                                 *(ctrl if ctrl else ("", "", ""))))
        rows.sort(key=lambda r: (r[1], r[4], r[3]))
        tmp = dest.with_suffix(".part")
        with gzip.open(tmp, "wt") as fh:
            fh.write("id\tchrom\tstart\tend\tstrand\tkind\tsequence\tgc\tctrl_start\tctrl_end\tctrl_gc\n")
            for i, r in enumerate(rows):
                ctrl_gc = r[9] if r[9] == "" else round(r[9], 3)
                fh.write(f"{chrom}_{i}\t" + "\t".join(map(str, (*r[:9], ctrl_gc))) + "\n")
        tmp.rename(dest)
        n_ctrl = sum(1 for r in rows if r[7] != "")
        print(f"{chrom}: {sum(r[4]=='G4' for r in rows)} G4, {sum(r[4]=='iM' for r in rows)} iM, "
              f"{n_ctrl} with a matched control", flush=True)


if __name__ == "__main__":
    main()
