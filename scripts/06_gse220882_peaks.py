#!/usr/bin/env python3
"""Call CUT&Tag peaks from the GSE220882 bigWigs and extract their DNA sequences.

Zanin, I. et al. (2023) *Genome-wide mapping of i-Motifs reveals their association
with transcription regulation in live human cells.* Nucleic Acids Research
51(16):8309-8321. doi:10.1093/nar/gkad626   GEO: GSE220882

This is the sequence-diversity source for the i-motif side of QuadCond. The
folding head trained on 160 designed C-tract oligos does not transfer to
peak-derived sequences (``scripts/07_validate_external.py`` measures the
separation at -0.020), and the reason is that the training set contains a
handful of parent constructs and their variants. These peaks are tens of
thousands of independent genomic loci.

Three things about the deposited data decide the design of this script.

**There are no BAMs, and no raw reads.** GSE220882_RAW contains 17 bigWigs and
three TPM tables. ``MACS2 -f BAMPE`` cannot run on a bigWig at all -- it needs
alignments, which means going back to SRA. So the peak caller has to work from
coverage, which is what SEACR does and what the authors used.

**SEACR itself cannot run here either.** It is a bash wrapper around an R script
and bedtools, and neither R nor bedtools is installed on the machine holding the
data, which has no network route to install them. What this script implements is
SEACR's *no-control* logic, in Python, and it says so on every record it writes:
contiguous non-zero coverage blocks, ranked by area under the curve. That is the
part of SEACR that does the work. It is a reimplementation, not SEACR, and the
records are flagged ``peak_caller=seacr_style_auc_blocks_reimplementation``.

**The cutoff does not have to be guessed.** SEACR run without an IgG control
takes the top *n* fraction of blocks, and choosing *n* is normally the arbitrary
step. Table S1 of the paper reports the exact number of SEACR peaks the authors
obtained for each of the twelve samples. This script ranks blocks by AUC exactly
as SEACR does and cuts at the authors' own published count for that sample, so
the free parameter is taken from the paper rather than invented here. The
implied top-fraction is reported per sample so the choice stays visible.

After calling, per (cell line, antibody):

* a region is kept only if it is called in **at least 2 of the 3 replicates**;
* regions overlapping the ENCODE hg38 blacklist v2 are dropped;
* the summit is the highest-signal position of the best-supported replicate block;
* the sequence is read out of a reference genome (see ``--genome-star``) as
  ``summit +/- flank``, and any window containing an N is dropped.

No lift-over is performed and none is needed: every bigWig in the series is
already on hg38 (chr1 = 248,956,422).

Output is JSONL, one record per peak, for ``quadcond.atlas.ingest.gse220882``.

Runs where the data is. It needs only numpy and pyBigWig.

    python scripts/06_gse220882_peaks.py \
        --bigwig-dir GSE220882_RAW \
        --genome-star /path/to/STAR_index \
        --blacklist hg38-blacklist.v2.bed \
        --out gse220882_peaks.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# The twelve CUT&Tag samples, with the authors' own peak counts (Table S1 of
# the supplementary information). The count is the cutoff, not a check.
# --------------------------------------------------------------------------
SAMPLES = [
    # gsm,          file stem,              cell,      antibody, target, rep, seacr_peaks
    ("GSM7507860", "GSM7507860_HEK_G4_rep1",   "HEK293T", "BG4",  "G4", 1, 45551),
    ("GSM7507861", "GSM7507861_HEK_G4_rep2",   "HEK293T", "BG4",  "G4", 2, 57905),
    ("GSM7507862", "GSM7507862_HEK_G4_rep3",   "HEK293T", "BG4",  "G4", 3, 63981),
    ("GSM7507863", "GSM7507863_HEK_iM_rep1",   "HEK293T", "iMab", "iM", 1, 41698),
    ("GSM7507864", "GSM7507864_HEK_iM_rep2",   "HEK293T", "iMab", "iM", 2, 60804),
    ("GSM7507865", "GSM7507865_HEK_iM_rep3",   "HEK293T", "iMab", "iM", 3, 53117),
    ("GSM7507866", "GSM7507866_WDLPS_G4_rep1", "WDLPS",   "BG4",  "G4", 1, 12458),
    ("GSM7507867", "GSM7507867_WDLPS_G4_rep2", "WDLPS",   "BG4",  "G4", 2, 14153),
    ("GSM7507868", "GSM7507868_WDLPS_G4_rep3", "WDLPS",   "BG4",  "G4", 3, 8138),
    ("GSM7507869", "GSM7507869_WDLPS_iM_rep1", "WDLPS",   "iMab", "iM", 1, 9031),
    ("GSM7507870", "GSM7507870_WDLPS_iM_rep2", "WDLPS",   "iMab", "iM", 2, 5469),
    ("GSM7507871", "GSM7507871_WDLPS_iM_rep3", "WDLPS",   "iMab", "iM", 3, 7263),
]

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DOI = "10.1093/nar/gkad626"
GEO = "GSE220882"


# --------------------------------------------------------------------------
# Reference sequence, read out of a STAR index
# --------------------------------------------------------------------------
class StarGenome:
    """Random-access reference sequence from a STAR genome directory.

    STAR stores the reference in ``Genome`` as one byte per base -- 0=A, 1=C,
    2=G, 3=T, 4=N, higher values for inter-chromosome padding -- with each
    chromosome starting at the offset given in ``chrStart.txt``. That makes the
    index a perfectly good random-access FASTA substitute, which matters here
    because the machine holding the data has a STAR index but no plain genome
    file and no network route to UCSC or Ensembl.

    :meth:`self_test` verifies the decoding against the index's own annotated
    splice junctions: real introns begin GT and end AG (or CT/AC read on the
    minus strand), so if a random sample of them decodes canonically the byte
    offsets and the alphabet are both right. Random bytes score ~0.4%; a
    correct decode scores ~96-99%, the shortfall being genuine GC-AG and AT-AC
    introns plus annotation noise.
    """

    CODE = "ACGTN"

    def __init__(self, star_dir: str | Path):
        d = Path(star_dir)
        names = [l.strip() for l in (d / "chrName.txt").read_text().split("\n") if l.strip()]
        starts = [int(l) for l in (d / "chrStart.txt").read_text().split() if l.strip()]
        self.offset = {n: starts[i] for i, n in enumerate(names)}
        self.length = {
            l.split()[0]: int(l.split()[1])
            for l in (d / "chrNameLength.txt").read_text().split("\n") if l.strip()
        }
        self._dir = d
        self._fh = open(d / "Genome", "rb")

    def fetch(self, chrom: str, start: int, end: int) -> str | None:
        """0-based half-open. Returns None if out of bounds."""
        n = self.length.get(chrom)
        if n is None or start < 0 or end > n or end <= start:
            return None
        self._fh.seek(self.offset[chrom] + start)
        raw = self._fh.read(end - start)
        return "".join(self.CODE[b] if b < 5 else "N" for b in raw)

    def self_test(self, n: int = 2000, seed: int = 0) -> dict:
        import random

        path = self._dir / "sjdbList.out.tab"
        if not path.exists():
            return {"skipped": "no sjdbList.out.tab"}
        rows = [l.split() for l in path.read_text().split("\n") if l.strip()]
        rng = random.Random(seed)
        rows = rng.sample(rows, min(n, len(rows)))
        canon = other = 0
        motifs: dict[str, int] = {}
        for c, s, e, strand in rows:
            if c not in self.offset:
                continue
            s0, e0 = int(s) - 1, int(e)
            d, a = self.fetch(c, s0, s0 + 2), self.fetch(c, e0 - 2, e0)
            if d is None or a is None:
                continue
            pair = f"{d}-{a}"
            motifs[pair] = motifs.get(pair, 0) + 1
            if (strand == "+" and pair == "GT-AG") or (strand == "-" and pair == "CT-AC") \
               or (strand not in "+-" and pair in ("GT-AG", "CT-AC")):
                canon += 1
            else:
                other += 1
        total = canon + other
        top = sorted(motifs.items(), key=lambda kv: -kv[1])[:5]
        return {
            "junctions_checked": total,
            "canonical_GT_AG": canon,
            "fraction_canonical": round(canon / total, 4) if total else None,
            "top_dinucleotide_pairs": top,
            "passes": bool(total and canon / total > 0.90),
        }


# --------------------------------------------------------------------------
# SEACR-style block calling
# --------------------------------------------------------------------------
def call_blocks(bw, chroms: list[str]) -> dict[str, np.ndarray]:
    """Contiguous non-zero coverage blocks with their AUC, max and summit.

    A bigWig stores coverage as runs of constant value. Runs that touch
    (``end == next start``) are one block; any gap starts a new one. This is
    SEACR's definition of a candidate peak.
    """
    ch, st, en, auc, mx, summit = [], [], [], [], [], []
    for c in chroms:
        try:
            iv = bw.intervals(c)
        except (RuntimeError, TypeError):
            iv = None
        if not iv:
            continue
        a = np.asarray(iv, dtype=np.float64)      # (n, 3): start, end, value
        s, e, v = a[:, 0], a[:, 1], a[:, 2]
        # a new block starts wherever this run does not begin where the last ended
        new = np.empty(len(s), dtype=bool)
        new[0] = True
        new[1:] = s[1:] != e[:-1]
        bid = np.cumsum(new) - 1
        nb = int(bid[-1]) + 1
        area = (e - s) * v
        block_auc = np.bincount(bid, weights=area, minlength=nb)
        block_max = np.zeros(nb)
        np.maximum.at(block_max, bid, v)
        first = np.flatnonzero(new)
        last = np.append(first[1:], len(s)) - 1
        # summit: midpoint of the highest-valued run in the block. Sorting by
        # (block, value) puts each block's peak run last within its group, so
        # the group end offsets index the summits directly -- a Python loop over
        # millions of blocks here costs more than the rest of the script.
        order = np.lexsort((v, bid))
        group_end = np.cumsum(np.bincount(bid, minlength=nb)) - 1
        top_run = order[group_end]
        ch.append(np.full(nb, c, dtype=object))
        st.append(s[first].astype(np.int64))
        en.append(e[last].astype(np.int64))
        auc.append(block_auc)
        mx.append(block_max)
        summit.append(((s[top_run] + e[top_run]) // 2).astype(np.int64))
    if not ch:
        return {k: np.array([]) for k in ("chrom", "start", "end", "auc", "max", "summit")}
    return {
        "chrom": np.concatenate(ch), "start": np.concatenate(st),
        "end": np.concatenate(en), "auc": np.concatenate(auc),
        "max": np.concatenate(mx), "summit": np.concatenate(summit),
    }


def top_k(blocks: dict[str, np.ndarray], k: int) -> dict[str, np.ndarray]:
    n = len(blocks["auc"])
    k = min(k, n)
    idx = np.argpartition(blocks["auc"], n - k)[n - k:]
    idx = idx[np.argsort(-blocks["auc"][idx])]
    return {key: val[idx] for key, val in blocks.items()}


# --------------------------------------------------------------------------
# Reproducibility across replicates
# --------------------------------------------------------------------------
def reproducible(peak_sets: list[dict[str, np.ndarray]], min_reps: int) -> list[dict]:
    """Merged regions covered by at least ``min_reps`` of the replicate peak sets.

    Per chromosome, a +1/-1 sweep over replicate intervals gives the depth of
    replicate support at every position; runs where depth >= min_reps are the
    reproducible regions. Each replicate contributes at most 1 to the depth at a
    position, which is what "called in 2 of 3 replicates" has to mean -- two
    overlapping blocks in the same replicate must not stand in for two
    replicates.
    """
    by_chrom: dict[str, list[list[tuple[int, int, int]]]] = {}
    for ri, ps in enumerate(peak_sets):
        for c, s, e in zip(ps["chrom"], ps["start"], ps["end"]):
            by_chrom.setdefault(c, [[] for _ in peak_sets])[ri].append((int(s), int(e), ri))

    out: list[dict] = []
    for c, per_rep in by_chrom.items():
        events: list[tuple[int, int]] = []
        for rep_intervals in per_rep:
            # flatten this replicate's own overlaps first, so it counts once
            merged: list[list[int]] = []
            for s, e, _ in sorted(rep_intervals):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            for s, e in merged:
                events.append((s, 1))
                events.append((e, -1))
        events.sort()
        depth = 0
        run_start = None
        for pos, delta in events:
            was = depth
            depth += delta
            if was < min_reps <= depth:
                run_start = pos
            elif was >= min_reps > depth and run_start is not None:
                out.append({"chrom": c, "start": run_start, "end": pos})
                run_start = None
    return out


def annotate(regions: list[dict], peak_sets: list[dict[str, np.ndarray]]) -> list[dict]:
    """Attach the best overlapping replicate block's summit, AUC and max.

    Both sides are sorted once and walked together; the regions came out of a
    coordinate sweep so they are already in order per chromosome.
    """
    index: dict[str, list[tuple[int, int, float, float, int, int]]] = {}
    for ri, ps in enumerate(peak_sets):
        for c, s, e, a, m, su in zip(ps["chrom"], ps["start"], ps["end"],
                                     ps["auc"], ps["max"], ps["summit"]):
            index.setdefault(c, []).append((int(s), int(e), float(a), float(m), int(su), ri))
    for v in index.values():
        v.sort()

    by_chrom: dict[str, list[dict]] = {}
    for r in regions:
        by_chrom.setdefault(r["chrom"], []).append(r)

    out: list[dict] = []
    for c, rs in by_chrom.items():
        cand = index.get(c, [])
        rs.sort(key=lambda r: r["start"])
        active: list[tuple[int, int, float, float, int, int]] = []
        i = 0
        for r in rs:
            while i < len(cand) and cand[i][0] < r["end"]:
                active.append(cand[i])
                i += 1
            active = [b for b in active if b[1] > r["start"]]
            best = None
            reps = set()
            for b in active:
                if b[0] < r["end"] and r["start"] < b[1]:
                    reps.add(b[5])
                    if best is None or b[2] > best[2]:
                        best = b
            if best is None:
                continue
            # The reproducible region is the >=2-replicate intersection, which is
            # narrower than the block it came from, so that block's summit can
            # fall outside it. Clamp, rather than centring a sequence window on a
            # coordinate the region does not contain.
            summit = min(max(best[4], r["start"]), r["end"] - 1)
            out.append({**r, "summit": summit, "summit_clamped": summit != best[4],
                        "auc": best[2], "max": best[3], "n_replicates": len(reps)})
    return out


class IntervalSet:
    """Merged per-chromosome intervals with an O(log n) overlap test.

    The start array is built once at construction. An earlier version rebuilt it
    inside the query, which turns a binary search into a linear scan and made the
    blacklist and cross-target steps quadratic -- on this data that is the
    difference between two seconds and never finishing.
    """

    def __init__(self, by_chrom: dict[str, list[tuple[int, int]]]):
        self.iv: dict[str, list[tuple[int, int]]] = {}
        self.starts: dict[str, list[int]] = {}
        for c, v in by_chrom.items():
            merged: list[list[int]] = []
            for s, e in sorted(v):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.iv[c] = [(s, e) for s, e in merged]
            self.starts[c] = [s for s, _ in merged]

    def hits(self, chrom: str, s: int, e: int) -> bool:
        import bisect

        v = self.iv.get(chrom)
        if not v:
            return False
        i = bisect.bisect_right(self.starts[chrom], e) - 1
        # merged intervals are disjoint, so only the one starting at or before e
        # can overlap [s, e)
        return i >= 0 and v[i][0] < e and s < v[i][1]

    @classmethod
    def from_bed(cls, path: str | Path) -> "IntervalSet":
        by: dict[str, list[tuple[int, int]]] = {}
        for line in Path(path).read_text().split("\n"):
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.split("\t")
            by.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        return cls(by)

    @classmethod
    def from_regions(cls, regions: list[dict]) -> "IntervalSet":
        by: dict[str, list[tuple[int, int]]] = {}
        for r in regions:
            by.setdefault(r["chrom"], []).append((r["start"], r["end"]))
        return cls(by)


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigwig-dir", required=True)
    ap.add_argument("--genome-star", required=True, help="STAR index directory (hg38)")
    ap.add_argument("--blacklist", required=True)
    ap.add_argument("--out", default="gse220882_peaks.jsonl")
    ap.add_argument("--flank", type=int, default=100, help="summit +/- flank -> 2*flank+1 nt")
    ap.add_argument("--min-reps", type=int, default=2)
    ap.add_argument("--relaxed-multiple", type=float, default=2.0,
                    help="the relaxed tier takes this multiple of the authors' peak count")
    ap.add_argument("--chroms", default=",".join(MAIN_CHROMS))
    a = ap.parse_args()

    import pyBigWig

    chroms = a.chroms.split(",")
    genome = StarGenome(a.genome_star)
    st = genome.self_test()
    print(f"genome self-test: {st}", flush=True)
    if not st.get("passes"):
        sys.exit("reference decoding failed its splice-junction check -- refusing to "
                 "write sequences that may be wrong")

    bl = IntervalSet.from_bed(a.blacklist)
    bwdir = Path(a.bigwig_dir)

    # ---- per-sample block calling -------------------------------------
    called: dict[tuple[str, str, str], list[dict]] = {}
    manifest = []
    for gsm, stem, cell, ab, target, rep, k in SAMPLES:
        path = next(iter(bwdir.glob(stem + ".*")), None)
        if path is None:
            print(f"  MISSING {stem}", flush=True)
            continue
        t0 = time.time()
        bw = pyBigWig.open(str(path))
        blocks = call_blocks(bw, chroms)
        bw.close()
        n_blocks = len(blocks["auc"])
        stringent = top_k(blocks, k)
        relaxed = top_k(blocks, int(k * a.relaxed_multiple))
        called.setdefault((cell, ab, "stringent"), []).append(stringent)
        called.setdefault((cell, ab, "relaxed"), []).append(relaxed)
        manifest.append({
            "gsm": gsm, "cell_line": cell, "antibody": ab, "target": target,
            "replicate": rep, "blocks": n_blocks,
            "authors_seacr_peaks": k,
            "implied_top_fraction": round(k / n_blocks, 5) if n_blocks else None,
            "auc_cutoff": float(stringent["auc"][-1]) if len(stringent["auc"]) else None,
        })
        print(f"  {stem:26s} blocks={n_blocks:>9,}  top-{k:,} "
              f"(={k/max(n_blocks,1):.3%})  {time.time()-t0:.0f}s", flush=True)

    # ---- reproducibility, blacklist, sequence -------------------------
    regions: dict[tuple[str, str, str], list[dict]] = {}
    for key, sets in called.items():
        cell, ab, tier = key
        reps = reproducible(sets, a.min_reps)
        reps = annotate(reps, sets)
        kept = [r for r in reps if not bl.hits(r["chrom"], r["start"], r["end"])]
        regions[key] = kept
        print(f"  {cell:8s} {ab:5s} {tier:9s} >= {a.min_reps}/3 reps: "
              f"{len(reps):>7,}  after blacklist: {len(kept):>7,}", flush=True)

    # cross-target overlap: an iM peak that is also a G4 peak is ambiguous
    other_idx = {}
    for (cell, ab, tier), rs in regions.items():
        other = "BG4" if ab == "iMab" else "iMab"
        other_idx[(cell, ab, tier)] = IntervalSet.from_regions(
            regions.get((cell, other, tier), []))

    n_written = n_dropped_n = 0
    with open(a.out, "w") as fh:
        for (cell, ab, tier), rs in sorted(regions.items()):
            target = "iM" if ab == "iMab" else "G4"
            for r in rs:
                s = r["summit"] - a.flank
                e = r["summit"] + a.flank + 1
                seq = genome.fetch(r["chrom"], s, e)
                if seq is None or "N" in seq:
                    n_dropped_n += 1
                    continue
                fh.write(json.dumps({
                    "sequence": seq,
                    "chrom": r["chrom"], "window_start": s, "window_end": e,
                    "peak_start": r["start"], "peak_end": r["end"], "summit": r["summit"],
                    "cell_line": cell, "antibody": ab, "target": target,
                    "tier": tier, "n_replicates": r["n_replicates"],
                    "auc": round(r["auc"], 3), "max_signal": round(r["max"], 3),
                    "summit_clamped": r["summit_clamped"],
                    "shares_locus_with_other_target":
                        other_idx[(cell, ab, tier)].hits(
                            r["chrom"], r["start"], r["end"]),
                    "geo": GEO, "doi": DOI, "genome_build": "hg38",
                }) + "\n")
                n_written += 1

    meta = {
        "geo": GEO, "doi": DOI, "genome_build": "hg38",
        "lift_over": "none required; the deposited bigWigs are already hg38",
        "peak_caller": "seacr_style_auc_blocks_reimplementation",
        "peak_caller_note":
            "SEACR itself (bash + R + bedtools) could not be run on the machine "
            "holding the data. Contiguous non-zero coverage blocks are ranked by "
            "AUC, as SEACR does without a control; the cutoff per sample is the "
            "peak count the authors report in Table S1 rather than a chosen "
            "top-fraction.",
        "min_replicates": a.min_reps,
        "blacklist": str(a.blacklist),
        "flank": a.flank,
        "genome_self_test": st,
        "samples": manifest,
        "regions": {f"{c}|{ab}|{t}": len(v) for (c, ab, t), v in sorted(regions.items())},
        "records_written": n_written,
        "windows_dropped_for_N": n_dropped_n,
    }
    Path(a.out).with_suffix(".meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote {n_written:,} records to {a.out}")
    print(f"wrote {Path(a.out).with_suffix('.meta.json')}")


if __name__ == "__main__":
    main()
