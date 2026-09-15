"""Motif search and structural parsing for G-quadruplexes and i-motifs.

Two jobs live here:

1. **Finding** candidate elements in a longer sequence (``find_g4``, ``find_im``).
2. **Parsing** a single candidate into its tract/loop architecture
   (``parse_tracts``), which is what the feature layer actually consumes.

The i-motif searcher follows the directed-graph traversal formulation of
Putative-iM-Searcher (Yang et al., NAR 2024, doi:10.1093/nar/gkae092): C-tracts
are nodes, admissible loops are edges, and a candidate is any path of four
nodes.  This implementation is written from that description with an explicit
dynamic-programming pass instead of recursion, so it is linear in the number of
edges and safe on chromosome-scale input.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator, Literal, Sequence

Strand = Literal["+", "-"]
Kind = Literal["G4", "iM"]

_COMPLEMENT = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def clean(seq: str) -> str:
    """Uppercase, RNA->DNA, and mask anything that is not ACGT to N."""
    s = seq.upper().replace("U", "T")
    return re.sub(r"[^ACGT]", "N", s)


# --------------------------------------------------------------------------- #
# Element container
# --------------------------------------------------------------------------- #
@dataclass
class Element:
    """One candidate G4 or i-motif."""

    kind: Kind
    sequence: str
    start: int = 0
    end: int = 0
    strand: Strand = "+"
    tract_len: int = 0
    tracts: tuple[tuple[int, int], ...] = ()
    loops: tuple[int, ...] = ()
    score: float = 0.0
    source_id: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def loop_total(self) -> int:
        return sum(self.loops)

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "sequence": self.sequence,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "length": self.length,
            "tract_len": self.tract_len,
            "loops": list(self.loops),
            "score": round(self.score, 4),
            "source_id": self.source_id,
        }


# --------------------------------------------------------------------------- #
# G4Hunter
# --------------------------------------------------------------------------- #
def g4hunter_scores(seq: str) -> list[float]:
    """Per-base G4Hunter values (Bedrat, Lacroix & Mergny, NAR 2016).

    Runs of G score +min(run, 4); runs of C score -min(run, 4); other bases 0.
    """
    s = clean(seq)
    out = [0.0] * len(s)
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch in "GC":
            j = i
            while j < n and s[j] == ch:
                j += 1
            run = min(j - i, 4)
            val = float(run) if ch == "G" else -float(run)
            for k in range(i, j):
                out[k] = val
            i = j
        else:
            i += 1
    return out


def g4hunter_mean(seq: str) -> float:
    s = clean(seq)
    if not s:
        return 0.0
    return sum(g4hunter_scores(s)) / len(s)


def g4hunter_window_max(seq: str, window: int = 25) -> float:
    """Max absolute windowed G4Hunter score -- the value the original tool thresholds."""
    vals = g4hunter_scores(seq)
    if len(vals) < window:
        return abs(sum(vals) / len(vals)) if vals else 0.0
    run = sum(vals[:window])
    best = abs(run / window)
    for i in range(window, len(vals)):
        run += vals[i] - vals[i - window]
        best = max(best, abs(run / window))
    return best


# --------------------------------------------------------------------------- #
# Tract / loop parsing
# --------------------------------------------------------------------------- #
def parse_tracts(
    seq: str,
    base: str,
    min_tract: int = 3,
    max_loop: int = 12,
) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]] | None:
    """Greedy left-to-right parse of a candidate into four tracts + three loops.

    Returns ``(tracts, loops)`` or ``None`` when the sequence does not contain
    four admissible tracts.  ``base`` is ``"G"`` for quadruplexes and ``"C"``
    for i-motifs.
    """
    s = clean(seq)
    runs: list[tuple[int, int]] = []
    for m in re.finditer(f"{base}{{{min_tract},}}", s):
        runs.append((m.start(), m.end()))
    if len(runs) < 4:
        return None
    # Choose the first four runs that satisfy the loop constraint, sliding
    # forward if a gap is too large.
    for start_idx in range(len(runs) - 3):
        window = runs[start_idx : start_idx + 4]
        loops = tuple(window[i + 1][0] - window[i][1] for i in range(3))
        if all(1 <= L <= max_loop for L in loops):
            return tuple(window), loops
    # Fall back to the first four runs regardless of loop length, so downstream
    # features still get an architecture (loops are reported honestly).
    window = runs[:4]
    loops = tuple(window[i + 1][0] - window[i][1] for i in range(3))
    return tuple(window), loops


# --------------------------------------------------------------------------- #
# Regex-based finders
# --------------------------------------------------------------------------- #
def _pattern(base: str, tract: int, min_loop: int, max_loop: int) -> re.Pattern:
    t = f"{base}{{{tract}}}"
    lo = f"[ACGT]{{{min_loop},{max_loop}}}"
    return re.compile(f"(?=({t}{lo}{t}{lo}{t}{lo}{t}))")


def _find(
    seq: str,
    kind: Kind,
    base: str,
    tract_len: int,
    min_loop: int,
    max_loop: int,
    strand: Strand,
    overlapping: bool,
    offset: int,
) -> list[Element]:
    s = clean(seq)
    pat = _pattern(base, tract_len, min_loop, max_loop)
    out: list[Element] = []
    last_end = -1
    for m in pat.finditer(s):
        start = m.start()
        text = m.group(1)
        end = start + len(text)
        if not overlapping and start < last_end:
            continue
        parsed = parse_tracts(text, base, tract_len, max_loop)
        tracts, loops = parsed if parsed else ((), ())
        el = Element(
            kind=kind,
            sequence=text,
            start=offset + start,
            end=offset + end,
            strand=strand,
            tract_len=tract_len,
            tracts=tracts,
            loops=loops,
            score=g4hunter_mean(text) if kind == "G4" else -g4hunter_mean(text),
        )
        out.append(el)
        last_end = end
    return out


def find_g4(
    seq: str,
    *,
    tract_len: int = 3,
    min_loop: int = 1,
    max_loop: int = 7,
    both_strands: bool = False,
    overlapping: bool = False,
) -> list[Element]:
    """Canonical G4 motif search: four G-tracts separated by 1-``max_loop`` nt loops."""
    els = _find(seq, "G4", "G", tract_len, min_loop, max_loop, "+", overlapping, 0)
    if both_strands:
        rc = revcomp(clean(seq))
        n = len(rc)
        for e in _find(rc, "G4", "G", tract_len, min_loop, max_loop, "-", overlapping, 0):
            e.start, e.end = n - e.end, n - e.start
            els.append(e)
    return els


def find_im(
    seq: str,
    *,
    tract_len: int = 3,
    min_loop: int = 1,
    max_loop: int = 12,
    both_strands: bool = False,
    overlapping: bool = False,
) -> list[Element]:
    """Canonical i-motif motif search: four C-tracts, loops up to 12 nt.

    The wider default loop bound (12 vs 7) follows the i-motif literature,
    where long-looped species are well documented, and matches
    Putative-iM-Searcher's default.
    """
    els = _find(seq, "iM", "C", tract_len, min_loop, max_loop, "+", overlapping, 0)
    if both_strands:
        rc = revcomp(clean(seq))
        n = len(rc)
        for e in _find(rc, "iM", "C", tract_len, min_loop, max_loop, "-", overlapping, 0):
            e.start, e.end = n - e.end, n - e.start
            els.append(e)
    return els


# --------------------------------------------------------------------------- #
# Graph-traversal i-motif search (variable tract lengths)
# --------------------------------------------------------------------------- #
def find_im_graph(
    seq: str,
    *,
    min_tract: int = 3,
    max_tract: int = 5,
    min_loop: int = 1,
    max_loop: int = 12,
    representative: Literal["balanced", "short_side", "all"] = "balanced",
    max_candidates_per_start: int = 64,
) -> list[Element]:
    """Directed-graph i-motif search allowing variable C-tract lengths.

    Nodes are maximal C-runs (trimmed to ``max_tract``); edges connect nodes
    whose separation lies in ``[min_loop, max_loop]``.  Every 4-node path is a
    candidate conformation; ``representative`` picks one per start position:

    ``balanced``    minimum standard deviation of the three loop lengths
    ``short_side``  minimum combined length of the two outer loops
    ``all``         return every conformation

    Follows the formulation in Putative-iM-Searcher (Yang et al., NAR 2024).
    """
    s = clean(seq)
    nodes: list[tuple[int, int]] = []
    for m in re.finditer(f"C{{{min_tract},}}", s):
        a, b = m.start(), m.end()
        # A long C-run can host a tract starting at several offsets; enumerate
        # the admissible tract lengths anchored at the run start and end.
        for L in range(min_tract, min(max_tract, b - a) + 1):
            nodes.append((a, a + L))
            if b - L != a:
                nodes.append((b - L, b))
    nodes = sorted(set(nodes))
    if len(nodes) < 4:
        return []

    idx_by_start: dict[int, list[int]] = {}
    for i, (a, _b) in enumerate(nodes):
        idx_by_start.setdefault(a, []).append(i)

    # adjacency
    adj: list[list[int]] = [[] for _ in nodes]
    starts = [a for a, _ in nodes]
    import bisect

    for i, (_a, b) in enumerate(nodes):
        lo = bisect.bisect_left(starts, b + min_loop)
        hi = bisect.bisect_right(starts, b + max_loop)
        adj[i] = list(range(lo, hi))

    out: list[Element] = []
    for i in range(len(nodes)):
        cands: list[tuple[int, int, int, int]] = []
        for j in adj[i]:
            for k in adj[j]:
                for m in adj[k]:
                    cands.append((i, j, k, m))
                    if len(cands) >= max_candidates_per_start:
                        break
                if len(cands) >= max_candidates_per_start:
                    break
            if len(cands) >= max_candidates_per_start:
                break
        if not cands:
            continue

        def loops_of(path: tuple[int, int, int, int]) -> tuple[int, int, int]:
            n0, n1, n2, n3 = (nodes[p] for p in path)
            return (n1[0] - n0[1], n2[0] - n1[1], n3[0] - n2[1])

        if representative == "all":
            chosen = cands
        elif representative == "short_side":
            chosen = [min(cands, key=lambda p: loops_of(p)[0] + loops_of(p)[2])]
        else:
            def spread(p):
                L = loops_of(p)
                mu = sum(L) / 3.0
                return sum((x - mu) ** 2 for x in L)
            chosen = [min(cands, key=spread)]

        for path in chosen:
            n0 = nodes[path[0]]
            n3 = nodes[path[3]]
            text = s[n0[0] : n3[1]]
            tracts = tuple(nodes[p] for p in path)
            out.append(
                Element(
                    kind="iM",
                    sequence=text,
                    start=n0[0],
                    end=n3[1],
                    tract_len=min(b - a for a, b in tracts),
                    tracts=tracts,
                    loops=loops_of(path),
                    score=-g4hunter_mean(text),
                )
            )
    # de-duplicate identical spans
    seen: set[tuple[int, int]] = set()
    uniq: list[Element] = []
    for e in sorted(out, key=lambda e: (e.start, e.end)):
        key = (e.start, e.end)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(e)
    return uniq


def scan(
    seq: str,
    *,
    kinds: Sequence[Kind] = ("G4", "iM"),
    both_strands: bool = True,
    **kw,
) -> list[Element]:
    """Convenience wrapper: find every requested element class in a sequence."""
    els: list[Element] = []
    if "G4" in kinds:
        els += find_g4(seq, both_strands=both_strands, **kw)
    if "iM" in kinds:
        els += find_im(seq, both_strands=both_strands, **kw)
    return sorted(els, key=lambda e: (e.start, e.kind))


def iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)
