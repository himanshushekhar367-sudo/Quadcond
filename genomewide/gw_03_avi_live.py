#!/usr/bin/env python3
"""Step 3, fallback: AVI scores from the live AlphaGenome API instead of the bundle.

Use this only for a pilot. One interval query per motif window and one per
control window is ~66,000 calls for chr22 alone and ~2.2 M genome-wide, which is
not a reasonable load on the service; `--max-motifs` caps it and defaults to a
small number. The published Tabix bundle (gw_03_avi.py) is the way to do this at
genome scale.

Motifs are sampled by structural interest, so the pilot is not wasted on windows
where nothing happens: every motif with at least one motif-destroying SNV is
eligible, ranked by its largest |delta|, and `--max-motifs` are taken at random
from that pool with a fixed seed (`--strategy random`) or from the top of it
(`--strategy top`, which biases the sample and is only for looking at examples).

Output is the same format as gw_03_avi.py, so gw_04_stats.py reads either.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gw_common import out_dir  # noqa: E402


def pick_motifs(base: Path, chrom: str, n: int, strategy: str, seed: int) -> pd.DataFrame:
    m = pd.read_csv(base / "motifs" / f"{chrom}.tsv.gz", sep="\t")
    sp = base / "structural" / f"{chrom}.tsv.gz"
    if not sp.exists():
        raise SystemExit(f"{sp} is missing: run gw_02_structural.py first")
    s = pd.read_csv(sp, sep="\t", usecols=["motif_id", "state", "delta"])
    s["delta"] = pd.to_numeric(s.delta, errors="coerce")
    agg = s.groupby("motif_id").agg(n_lost=("state", lambda v: (v == "motif_lost").sum()),
                                    max_abs_delta=("delta", lambda v: v.abs().max()))
    agg = agg[agg.n_lost > 0].sort_values("max_abs_delta", ascending=False)
    ids = agg.index[:n] if strategy == "top" else agg.sample(min(n, len(agg)), random_state=seed).index
    return m[m.id.isin(set(ids))].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=["chr22"])
    ap.add_argument("--max-motifs", type=int, default=500, help="per chromosome")
    ap.add_argument("--strategy", choices=["random", "top"], default="random")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sleep", type=float, default=0.0, help="seconds between queries")
    ap.add_argument("--api-key", help="defaults to $ALPHAGENOME_API_KEY")
    ap.add_argument("--scorers", default="AVI_SCORE",
                    help="comma-separated Atlas scorers. The default is AVI_SCORE alone, "
                         "which is the axis this analysis uses; asking for the full default "
                         "set returns megabytes of per-track scores per window and is what "
                         "made these queries fail")
    ap.add_argument("--max-message-mb", type=int, default=256,
                    help="gRPC receive limit for the Atlas channel")
    a = ap.parse_args()
    from quadcond import alphagenome as ag
    src = ag.LiveAtlas.from_api_key(
        a.api_key, scorers=[x.strip() for x in a.scorers.split(",") if x.strip()],
        max_message_mb=a.max_message_mb)
    base = out_dir()
    od = base / "avi"
    od.mkdir(exist_ok=True)
    for chrom in a.chroms:
        dest = od / f"{chrom}.tsv.gz"
        if dest.exists():
            print(f"{chrom}: exists, skipped")
            continue
        m = pick_motifs(base, chrom, a.max_motifs, a.strategy, a.seed)
        print(f"{chrom}: {len(m)} motifs selected ({a.strategy}); "
              f"{2 * len(m)} interval queries", flush=True)
        tmp = dest.with_suffix(".part")
        n, failed, scorers = 0, 0, []
        with gzip.open(tmp, "wt") as fh:
            header_written = False
            for j, r in enumerate(m.itertuples(index=False)):
                wins = [("motif", int(r.start) + 1, int(r.end))]
                if str(r.ctrl_start) not in ("", "nan"):
                    wins.append(("control", int(r.ctrl_start) + 1, int(r.ctrl_end)))
                for tag, s1, e1 in wins:
                    try:
                        recs = src.records_for_interval(chrom, s1, e1)
                    except ag.AtlasUnavailable as exc:
                        failed += 1
                        if failed <= 3:
                            print(f"  query failed ({tag} {chrom}:{s1}-{e1}): {exc}", flush=True)
                            if "RESOURCE_EXHAUSTED" in str(exc) or "larger than max" in str(exc):
                                print("  hint: that is the gRPC receive cap, not the service. "
                                      "Keep --scorers to AVI_SCORE (the default) and/or raise "
                                      "--max-message-mb.", flush=True)
                        if failed > 25:
                            raise SystemExit("too many failed queries; stopping rather than "
                                             "writing a table with silent holes")
                        continue
                    for key, rec in recs.items():
                        if not header_written:
                            scorers = sorted(rec.scores)
                            fh.write("set\tmotif_id\tchrom\tpos\tref\talt\t" + "\t".join(scorers) + "\n")
                            header_written = True
                        vals = [f"{rec.scores.get(s, '')}" for s in scorers]
                        fh.write(f"{tag}\t{r.id}\t{chrom}\t{key.position}\t{key.reference}\t"
                                 f"{key.alternate}\t" + "\t".join(vals) + "\n")
                        n += 1
                    if a.sleep:
                        time.sleep(a.sleep)
                if j % 25 == 0:
                    print(f"  {chrom}: {j + 1}/{len(m)} motifs, {n} records", flush=True)
            if not header_written:
                raise SystemExit("no records returned; nothing written")
        tmp.rename(dest)
        print(f"{chrom}: {n} AVI records from the live API ({failed} failed queries)", flush=True)
    print("NOTE: this is a sampled pilot, not the genome-wide set. The statistics "
          "step will say so only if you do -- record --max-motifs and --strategy "
          "in your methods.")


if __name__ == "__main__":
    main()
