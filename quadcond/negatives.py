"""Negative-set construction.

A G4/i-motif classifier is only as honest as its negatives.  Random or
genome-background negatives make any model look excellent because the task
collapses to "is this sequence G-rich".  QuadCond therefore builds negatives
that are *composition-matched* to the positives:

``dinucleotide_shuffle``  Altschul-Erikson Eulerian-path shuffle, which
                          preserves exact mononucleotide **and** dinucleotide
                          counts.  A shuffled hTelo repeat still has the same
                          G content and the same GG frequency; what it loses is
                          the four-tract register.

``matched_negatives``     one or more shuffles per positive, optionally
                          rejecting shuffles that themselves still contain a
                          canonical motif (``reject_motif=True``) so the
                          negative class is genuinely non-forming.

The rejection option is a deliberate trade-off and is reported in the model
card: rejecting motif-positive shuffles makes the task cleaner but slightly
easier; keeping them makes the classifier learn architecture rather than
composition.  Default is to keep them (``reject_motif=False``) and let the
model earn its score.
"""
from __future__ import annotations

import random
from collections import defaultdict
from typing import Iterable, Sequence

from .motifs import clean, find_g4, find_im


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erikson dinucleotide-preserving shuffle."""
    s = clean(seq)
    if len(s) < 4:
        return s
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(s, s[1:]):
        edges[a].append(b)
    first, last = s[0], s[-1]

    # Build a random Eulerian path: choose a last-edge per vertex forming a
    # tree rooted at `last`, then shuffle the remainder.
    for _ in range(100):
        last_edge: dict[str, str] = {}
        ok = True
        for v, outs in edges.items():
            if v == last:
                continue
            last_edge[v] = rng.choice(outs)
        # connectivity check: following last_edge from every vertex must reach `last`
        for v in last_edge:
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in last_edge:
                    ok = False
                    break
                seen.add(cur)
                cur = last_edge[cur]
            if not ok:
                break
        if ok:
            break
    else:  # pragma: no cover - fall back to mononucleotide shuffle
        chars = list(s)
        rng.shuffle(chars)
        return "".join(chars)

    order: dict[str, list[str]] = {}
    for v, outs in edges.items():
        rest = list(outs)
        if v in last_edge:
            rest.remove(last_edge[v])
            rng.shuffle(rest)
            rest.append(last_edge[v])
        else:
            rng.shuffle(rest)
        order[v] = rest

    out = [first]
    cur = first
    idx: dict[str, int] = defaultdict(int)
    for _ in range(len(s) - 1):
        nxt = order[cur][idx[cur]]
        idx[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def has_motif(seq: str, kind: str) -> bool:
    if kind == "G4":
        return bool(find_g4(seq))
    return bool(find_im(seq))


def matched_negatives(
    positives: Sequence[str],
    *,
    kind: str = "G4",
    n_per_positive: int = 1,
    reject_motif: bool = False,
    max_tries: int = 25,
    seed: int = 0,
) -> list[str]:
    """Composition-matched negatives, one or more per positive sequence."""
    rng = random.Random(seed)
    out: list[str] = []
    for p in positives:
        made = 0
        for _ in range(max_tries * n_per_positive):
            if made >= n_per_positive:
                break
            cand = dinucleotide_shuffle(p, rng)
            if cand == clean(p):
                continue
            if reject_motif and has_motif(cand, kind):
                continue
            out.append(cand)
            made += 1
        while made < n_per_positive:  # give up on rejection, keep balance
            out.append(dinucleotide_shuffle(p, rng))
            made += 1
    return out


def composition_report(pos: Iterable[str], neg: Iterable[str]) -> dict[str, float]:
    """Sanity check that negatives really are composition-matched."""
    import numpy as np

    def stats(seqs):
        seqs = [clean(s) for s in seqs]
        gc = [(s.count("G") + s.count("C")) / max(len(s), 1) for s in seqs]
        g = [s.count("G") / max(len(s), 1) for s in seqs]
        ln = [len(s) for s in seqs]
        return float(np.mean(gc)), float(np.mean(g)), float(np.mean(ln))

    pgc, pg, pl = stats(pos)
    ngc, ng, nl = stats(neg)
    return {
        "pos_gc": pgc, "neg_gc": ngc, "delta_gc": pgc - ngc,
        "pos_g": pg, "neg_g": ng, "delta_g": pg - ng,
        "pos_len": pl, "neg_len": nl, "delta_len": pl - nl,
    }
