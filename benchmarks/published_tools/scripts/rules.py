"""Faithful re-implementations of rule-based G4 / iM scores (documented algorithms)."""
import re, itertools, numpy as np

def g4hunter_vec(seq):
    s = seq.upper(); out = np.zeros(len(s)); i = 0
    while i < len(s):
        c = s[i]
        if c in "GC":
            j = i
            while j < len(s) and s[j] == c: j += 1
            v = min(j - i, 4) * (1 if c == "G" else -1)
            out[i:j] = v; i = j
        else:
            i += 1
    return out

def g4hunter(seq, window=25):
    """Bedrat et al. 2016. Returns (whole-sequence mean, max window, min window)."""
    v = g4hunter_vec(seq)
    if len(v) == 0: return 0.0, 0.0, 0.0
    w = min(window, len(v))
    cs = np.concatenate([[0], np.cumsum(v)])
    win = (cs[w:] - cs[:-w]) / w
    return float(v.mean()), float(win.max()), float(win.min())

G4_CANON = re.compile(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))")
G4_LONG = re.compile(r"(?=(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,}))")
G2_REL = re.compile(r"(?=(G{2,}\w{1,7}G{2,}\w{1,7}G{2,}\w{1,7}G{2,}))")
IM_CANON = re.compile(r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))")
IM_SEEKER_LIKE = re.compile(r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))")

def regex_hits(seq, pat):
    return len(pat.findall(seq.upper()))

def qgrs_gscore(seq, max_len=30, min_g=2):
    """Kikin et al. 2006 QGRS Mapper G-score (re-implementation of qgrs-cpp scoring):
    G = gmax - gavg + gmax*(tetrads-2); gmax = max_len - (min_g*4 + 1);
    gavg = mean pairwise |loop_i - loop_j|. Returns best score over all QGRS."""
    s = seq.upper(); n = len(s); best = 0.0
    gmax = max_len - (min_g * 4 + 1)
    gr = [(m.start(), m.end()) for m in re.finditer(r"G{%d,}" % min_g, s)]
    starts = sorted({p for a, b in gr for p in range(a, b)})
    sset = set(starts)
    def is_run(p, t):
        return p + t <= n and s[p:p + t] == "G" * t
    for t in range(min_g, 8):
        L4 = 4 * t
        if L4 > max_len: break
        cand = [p for p in starts if is_run(p, t)]
        cset = set(cand)
        for p1 in cand:
            for l1 in range(0 if t == 2 else 1, max_len - L4 + 1):
                p2 = p1 + t + l1
                if p2 not in cset: continue
                for l2 in range(0 if t == 2 else 1, max_len - L4 - l1 + 1):
                    p3 = p2 + t + l2
                    if p3 not in cset: continue
                    for l3 in range(0 if t == 2 else 1, max_len - L4 - l1 - l2 + 1):
                        p4 = p3 + t + l3
                        if p4 not in cset: continue
                        loops = (l1, l2, l3)
                        if sum(1 for x in loops if x == 0) > 1: continue
                        gavg = (abs(l1 - l2) + abs(l2 - l3) + abs(l1 - l3)) / 3.0
                        sc = gmax - gavg + gmax * (t - 2)
                        if sc > best: best = sc
    return best

def gc(seq):
    s = seq.upper(); return (s.count("G") + s.count("C")) / max(1, len(s))

def revcomp(s):
    return s.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
