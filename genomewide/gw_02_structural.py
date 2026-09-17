#!/usr/bin/env python3
"""Step 2: QuadCond structural effect of every SNV inside every motif.

For each motif (in its own orientation) every single-base substitution is made,
the motif rule is re-applied, and the structural head is evaluated on wild type
and mutant under one stated buffer:
  G4 -> g4_tm  (degC; condition-aware)
  iM -> im_pht (pH units; reference-buffer head, no condition response)
Motif loss is recorded as state `motif_lost` with no delta, exactly as the
package does. Variants are written on the + strand (VCF convention, 1-based).

Output: $GW_OUT/structural/<chrom>.tsv.gz
  chrom pos ref alt motif_id kind strand state wt delta
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gw_common import CHROMS, out_dir, revcomp  # noqa: E402

COND = dict(k=140.0, na=10.0, mg=1.0, ph=7.4, temperature=37.0)
_P = None


def _pred():
    global _P
    if _P is None:
        from quadcond import assets
        from quadcond.models.predict import Predictor
        _P = Predictor.load(assets.resolve("model"), None)
    return _P


def _value(entry):
    if not entry or entry.get("refused") or (entry.get("applicability") or {}).get("refused"):
        return None
    return entry.get("value")


def score_chunk(args):
    chunk, cond_kw = args
    from quadcond.conditions import Condition
    from quadcond.motifs import find_g4, find_im
    pred = _pred()
    cond = Condition(**cond_kw)
    out = []
    for kind in ("G4", "iM"):
        sub = chunk[chunk.kind == kind]
        if sub.empty:
            continue
        head = "g4_tm" if kind == "G4" else "im_pht"
        finder = find_g4 if kind == "G4" else find_im
        wts, muts, meta = [], [], []
        for r in sub.itertuples(index=False):
            s = r.sequence
            wts.append(s)
            L = len(s)
            for i in range(L):
                for b in "ACGT":
                    if b == s[i]:
                        continue
                    m = s[:i] + b + s[i + 1:]
                    lost = not finder(m)
                    if r.strand == "+":
                        pos, ref, alt = r.start + i + 1, s[i], b
                    else:
                        pos, ref, alt = r.end - i, revcomp(s[i]), revcomp(b)
                    meta.append((r.chrom, pos, ref, alt, r.id, kind, r.strand, lost, len(wts) - 1))
                    muts.append(m)
        wv = [_value(p["predictions"].get(head)) for p in
              pred.predict(wts, cond, heads=[head], n_neighbours=0)]
        keep = [k for k, mt in enumerate(meta) if not mt[7]]
        mv = [None] * len(muts)
        if keep:
            res = pred.predict([muts[k] for k in keep], cond, heads=[head], n_neighbours=0)
            for k, p in zip(keep, res):
                mv[k] = _value(p["predictions"].get(head))
        for mt, v in zip(meta, mv):
            chrom, pos, ref, alt, mid, kd, strand, lost, wi = mt
            w = wv[wi]
            if lost:
                out.append((chrom, pos, ref, alt, mid, kd, strand, "motif_lost", w, ""))
            elif w is None or v is None:
                out.append((chrom, pos, ref, alt, mid, kd, strand, "refused", w, ""))
            else:
                out.append((chrom, pos, ref, alt, mid, kd, strand, "motif_retained", w, round(v - w, 3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--chunk", type=int, default=400, help="motifs per task")
    ap.add_argument("--sample", type=float, default=1.0, help="fraction of motifs (pilot runs)")
    for k, v in COND.items():
        ap.add_argument(f"--{k}", type=float, default=v)
    a = ap.parse_args()
    cond = {k: getattr(a, k) for k in COND}
    base = out_dir()
    od = base / "structural"
    od.mkdir(exist_ok=True)
    with ProcessPoolExecutor(a.workers) as ex:
        for chrom in a.chroms:
            dest = od / f"{chrom}.tsv.gz"
            if dest.exists():
                print(f"{chrom}: exists, skipped")
                continue
            m = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t")
            if a.sample < 1:
                m = m.sample(frac=a.sample, random_state=0).sort_values("start")
            chunks = [(m.iloc[i:i + a.chunk], cond) for i in range(0, len(m), a.chunk)]
            tmp = dest.with_suffix(".part")
            n = 0
            with gzip.open(tmp, "wt") as fh:
                fh.write("chrom\tpos\tref\talt\tmotif_id\tkind\tstrand\tstate\twt\tdelta\n")
                for k, rows in enumerate(ex.map(score_chunk, chunks)):
                    for r in rows:
                        fh.write("\t".join("" if x is None else str(x) for x in r) + "\n")
                    n += len(rows)
                    if k % 20 == 0:
                        print(f"  {chrom}: {k + 1}/{len(chunks)} chunks", flush=True)
            tmp.rename(dest)
            print(f"{chrom}: {len(m)} motifs, {n} SNVs scored", flush=True)


if __name__ == "__main__":
    main()
