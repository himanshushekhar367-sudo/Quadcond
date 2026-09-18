#!/usr/bin/env python3
"""Step 3: AlphaGenome AVI scores for every motif SNV and every control SNV.

Reads the downloaded AVI Tabix bundle by region (no API key, no network). The
column layout is taken from the file's own header by quadcond.alphagenome, so a
publisher's column order is never assumed. `--avi` may be one file, or a pattern
containing {chrom} when the bundle is split per chromosome.

Output: $GW_OUT/avi/<chrom>.tsv.gz
  set(motif|control) motif_id chrom pos ref alt <one column per AVI score column>
Control SNVs are all 3 alternates at every base of the matched control window.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gw_common import CHROMS, out_dir  # noqa: E402


class FastTabix:
    """quadcond.alphagenome.TabixAtlas column logic with one persistent handle."""

    def __init__(self, path):
        import pysam
        from quadcond.alphagenome import TabixAtlas
        self.meta = TabixAtlas.from_path(path)           # validates header + index
        self.tf = pysam.TabixFile(str(path))
        self.contigs = set(self.tf.contigs)
        self.cols = self.meta._columns
        self.scores = self.meta._score_columns

    def fetch(self, chrom, start1, end1):
        name = chrom if chrom in self.contigs else chrom.replace("chr", "")
        if name not in self.contigs:
            return
        c = self.cols
        for line in self.tf.fetch(name, start1 - 1, end1):
            p = line.split("\t")
            try:
                pos = int(p[c["position"]])
            except (ValueError, IndexError):
                continue
            ref, alt = p[c["reference"]].upper(), p[c["alternate"]].upper()
            if len(ref) != 1 or len(alt) != 1 or not (start1 <= pos <= end1):
                continue
            yield pos, ref, alt, [p[i] if i < len(p) else "" for i in self.scores.values()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--avi", required=True, help="AVI .tsv.gz (with .tbi/.csi), or a pattern with {chrom}")
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--flank", type=int, default=100,
                    help="nt each side of a motif, written as set=flank: the within-locus control")
    a = ap.parse_args()
    base = out_dir()
    od = base / "avi"
    od.mkdir(exist_ok=True)
    handles = {}
    for chrom in a.chroms:
        dest = od / f"{chrom}.tsv.gz"
        if dest.exists():
            print(f"{chrom}: exists, skipped")
            continue
        path = a.avi.format(chrom=chrom)
        if path not in handles:
            handles[path] = FastTabix(path)
        tb = handles[path]
        m = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t")
        tmp = dest.with_suffix(".part")
        n = 0
        with gzip.open(tmp, "wt") as fh:
            fh.write("set\tmotif_id\tchrom\tpos\tref\talt\t" + "\t".join(tb.scores) + "\n")
            for r in m.itertuples(index=False):
                m1, m2 = int(r.start) + 1, int(r.end)
                wins = [("motif", max(1, m1 - a.flank), m2 + a.flank, (m1, m2))]
                if pd.notna(r.ctrl_start) and str(r.ctrl_start) != "":
                    wins.append(("control", int(r.ctrl_start) + 1, int(r.ctrl_end), None))
                for tag, s1, e1, bounds in wins:
                    for pos, ref, alt, vals in tb.fetch(chrom, s1, e1):
                        row_tag = tag
                        if bounds is not None:
                            row_tag = "motif" if bounds[0] <= pos <= bounds[1] else "flank"
                        fh.write(f"{row_tag}\t{r.id}\t{chrom}\t{pos}\t{ref}\t{alt}\t"
                                 + "\t".join(vals) + "\n")
                        n += 1
        tmp.rename(dest)
        print(f"{chrom}: {n} AVI records written", flush=True)


if __name__ == "__main__":
    main()
