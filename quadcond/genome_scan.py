"""Genomic G-quadruplex window scanner, trained on G4-seq rather than on oligos.

Why this exists
---------------
``g4_fold`` is trained on short, purified oligonucleotides (<= 99 nt) against
dinucleotide-shuffled negatives. It answers "does this oligo fold, compared with
a shuffle of itself". It does not answer "which windows of a genome form a G4":
benchmarked on human G4-seq it scored AUROC 0.66-0.69 against random genomic
windows and 0.29-0.48 against motif-matching windows G4-seq did not observe,
while G4Hunter and pqsfinder reached 0.92-0.97 on the first task
(``benchmarks/published_tools``).

This module is the genomic path. Its label is an experimental observation of a
different kind: polymerase stalling on purified human genomic DNA in K+ buffer
(G4-seq, Marsico et al. 2019, GEO GSE110582), as released in the G4detector
benchmark. The positives are observed G4-seq windows; the negatives are random
genomic windows and -- the hard class -- canonical-motif windows that G4-seq did
not observe. Chromosomes 2, 8 and 17 are held out of training, and a mouse
G4-seq set is kept as a cross-species test.

What a score means
------------------
P(window resembles an observed K+ G4-seq window) at the class balance of the
training set (50 % positives). It is a genomic in-vitro observation, not a
folding measurement under a user-chosen buffer, so the condition vector is not
an input and no condition claim is made. Use ``g4_tm`` / ``g4_topology`` on the
motifs this scanner finds for buffer-dependent stability and topology.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .conditions import Condition
from .features import featurize
from .motifs import clean, g4hunter_scores, revcomp

#: Window length of the G4-seq benchmark the scanner was trained on.
WINDOW = 124

_G3 = re.compile(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))")
_G3L = re.compile(r"(?=(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,}))")
_G2 = re.compile(r"(?=(G{2,}\w{1,7}G{2,}\w{1,7}G{2,}\w{1,7}G{2,}))")
_RUNS = re.compile(r"G{3,}")

CLAIM = (
    "P(window resembles an observed G4-seq window, human genomic DNA, K+ buffer). "
    "GENOMIC IN-VITRO OBSERVATION: polymerase stalling on purified genomic DNA "
    "(G4-seq, GSE110582), not a folding measurement under a chosen buffer. "
    "Calibrated at 50 % prevalence; recalibrate to your own background before "
    "reading it as an absolute probability."
)


def _window_stats(s: str) -> list[float]:
    v = np.asarray(g4hunter_scores(s), float)
    if len(v) == 0:
        return [0.0] * 10
    w = min(25, len(v))
    cs = np.concatenate([[0.0], np.cumsum(v)])
    win = (cs[w:] - cs[:-w]) / w
    runs = [len(m.group(0)) for m in _RUNS.finditer(s)]
    return [
        float(win.max()), float(win.min()), float(np.mean(v)),
        float(len(_G3.findall(s))), float(len(_G3L.findall(s))), float(len(_G2.findall(s))),
        float(len(runs)), float(max(runs) if runs else 0), float(sum(runs)),
        float((win >= 1.2).sum()) / len(win),
    ]


WINDOW_STAT_NAMES = ["g4h_win_max25", "g4h_win_min25", "g4h_mean_all", "n_g3_canonical",
                     "n_g3_longloop", "n_g2_relaxed", "n_g3_runs", "max_g_run",
                     "total_g3_run_nt", "frac_windows_ge_1p2"]


def window_features(seqs: Sequence[str]) -> np.ndarray:
    """Strand-symmetric features: the G-richer strand first, the other second."""
    seqs = [clean(s) for s in seqs]
    fw, rv = [], []
    for s in seqs:
        r = revcomp(s)
        a, b = (s, r) if sum(g4hunter_scores(s)) >= 0 else (r, s)
        fw.append(a)
        rv.append(b)
    conds = [Condition(k=100)] * len(seqs)
    Xa = featurize(fw, conds, kind="G4", use_conditions=False)
    Xb = featurize(rv, conds, kind="G4", use_conditions=False)
    Sa = np.array([_window_stats(s) for s in fw])
    Sb = np.array([_window_stats(s) for s in rv])
    return np.hstack([Xa, Sa, Xb, Sb]).astype(np.float32)


@dataclass
class G4SeqScanner:
    estimators: list
    calibrator: object | None
    metrics: dict = field(default_factory=dict)
    training_meta: dict = field(default_factory=dict)
    version: str = "g4seq-scanner-1"
    window: int = WINDOW

    # ------------------------------------------------------------ persistence
    def save(self, path: str | Path) -> Path:
        import joblib
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, p, compress=3)
        return p

    @classmethod
    def load(cls, path: str | Path) -> "G4SeqScanner":
        import joblib
        obj = joblib.load(path)
        if not isinstance(obj, cls):
            raise TypeError(f"{path} does not hold a G4SeqScanner")
        return obj

    # --------------------------------------------------------------- scoring
    def _raw(self, X: np.ndarray) -> np.ndarray:
        return np.mean([m.predict_proba(X)[:, 1] for m in self.estimators], axis=0)

    def score_windows(self, seqs: Sequence[str]) -> np.ndarray:
        """Calibrated score for each window (best at ~124 nt; 60-200 nt accepted)."""
        X = window_features(seqs)
        p = self._raw(X)
        if self.calibrator is not None:
            p = self.calibrator.transform(p)
        return np.clip(p, 0.0, 1.0)

    def scan(self, sequence: str, *, step: int = 25, threshold: float = 0.5,
             name: str = "query", rescue_g4hunter: float = 1.2) -> list[dict]:
        """Slide a 124-nt window along ``sequence``; merge windows above ``threshold``.

        Returns one record per merged region with its best window score, its
        ``call_basis`` and the canonical motifs inside it (both strands), ready to
        hand to ``g4_tm``. Canonical motifs with |G4Hunter| >= ``rescue_g4hunter``
        that no scanner region covers are added with call_basis
        ``canonical_motif_rescue`` (set it above 4 to disable).
        """
        from .motifs import find_g4

        s = clean(sequence)
        if len(s) < 60:
            raise ValueError("genome scan needs at least 60 nt; use `quadcond predict` for oligos")
        starts = list(range(0, max(1, len(s) - self.window + 1), step))
        if starts[-1] != max(0, len(s) - self.window):
            starts.append(max(0, len(s) - self.window))
        wins = [s[i:i + self.window] for i in starts]
        scores = self.score_windows(wins)
        regions: list[dict] = []
        for st, sc in zip(starts, scores):
            en = min(len(s), st + self.window)
            if sc < threshold:
                continue
            if regions and st <= regions[-1]["end"]:
                r = regions[-1]
                r["end"] = max(r["end"], en)
                r["best_score"] = max(r["best_score"], float(sc))
                r["n_windows"] += 1
            else:
                regions.append(dict(name=name, start=st, end=en, best_score=float(sc), n_windows=1))
        for r in regions:
            r["call_basis"] = "g4seq_scanner"
        # Canonical-motif rescue. The scanner learned what G4-seq observes in the
        # human genome, which is mostly G-rich neighbourhoods; an isolated
        # canonical motif in a G-poor context (an implanted human telomeric
        # repeat, for example) scores low even though it folds readily in vitro.
        # Such motifs are added as their own regions, labelled by call basis,
        # with the scanner score of the window centred on them reported as is.
        covered = [(r["start"], r["end"]) for r in regions]
        for e in find_g4(s, both_strands=True):
            if any(a <= e.start and e.end <= b for a, b in covered):
                continue
            if abs(e.score) < rescue_g4hunter:
                continue
            c = (e.start + e.end) // 2
            ws = max(0, min(len(s) - self.window, c - self.window // 2))
            sc_c = float(self.score_windows([s[ws:ws + self.window]])[0])
            regions.append(dict(name=name, start=e.start, end=e.end, best_score=sc_c,
                                n_windows=1, call_basis="canonical_motif_rescue"))
            covered.append((e.start, e.end))
        regions.sort(key=lambda r: r["start"])
        for r in regions:
            sub = s[r["start"]:r["end"]]
            r["start_1based"] = r["start"] + 1
            r["motifs"] = [
                {**e.as_dict(), "start": e.start + r["start"], "end": e.end + r["start"]}
                for e in find_g4(sub, both_strands=True)
            ]
        return regions

    def describe(self) -> dict:
        return {"name": "g4seq_scanner", "version": self.version, "window": self.window,
                "claim": CLAIM, "target_semantics": "genomic_in_vitro",
                "biophysically_grounded": False, "metrics": self.metrics,
                "training": self.training_meta}


def iter_windows(fasta_records: Iterable[tuple[str, str]], step: int = 25):
    for name, seq in fasta_records:
        s = clean(seq)
        for i in range(0, max(1, len(s) - WINDOW + 1), step):
            yield name, i, s[i:i + WINDOW]
