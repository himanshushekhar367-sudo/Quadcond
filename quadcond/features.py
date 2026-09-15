"""Feature construction.

Three blocks, kept separable so the ablation study (sequence-only vs
sequence+condition) is a one-line change:

``sequence_features``   composition, k-mers, tract/loop architecture
``condition_features``  the buffer, transformed into learnable coordinates
``interaction_features``the products that carry the physics we expect

Everything returns plain ``list[float]`` with a matching name vector, so the
same code path feeds gradient-boosted trees and any neural backbone.
"""
from __future__ import annotations

import itertools
import math
from typing import Iterable, Sequence

import numpy as np

from .conditions import Condition
from .motifs import clean, g4hunter_mean, g4hunter_window_max, parse_tracts

# Cytosine N3 pKa: the hemiprotonated C.C+ pair that builds an i-motif titrates
# around here, so (pKa - pH) is the natural pH coordinate for i-motif models.
CYTOSINE_N3_PKA = 4.6

_K1 = ["A", "C", "G", "T"]
_K2 = ["".join(p) for p in itertools.product(_K1, repeat=2)]
_K3 = ["".join(p) for p in itertools.product(_K1, repeat=3)]


# --------------------------------------------------------------------------- #
# Sequence block
# --------------------------------------------------------------------------- #
def _kmer_freq(s: str, kmers: Sequence[str], k: int) -> list[float]:
    n = max(len(s) - k + 1, 1)
    counts = dict.fromkeys(kmers, 0)
    for i in range(len(s) - k + 1):
        sub = s[i : i + k]
        if sub in counts:
            counts[sub] += 1
    return [counts[km] / n for km in kmers]


def _run_stats(s: str, base: str) -> tuple[int, int, float]:
    """(number of runs >= 2, longest run, mean run length) for a given base."""
    runs, cur = [], 0
    for ch in s:
        if ch == base:
            cur += 1
        else:
            if cur:
                runs.append(cur)
            cur = 0
    if cur:
        runs.append(cur)
    long_runs = [r for r in runs if r >= 2]
    return (
        len(long_runs),
        max(runs) if runs else 0,
        float(np.mean(runs)) if runs else 0.0,
    )


def _base_feature_names(kind: str) -> list[str]:
    base = "G" if kind == "G4" else "C"
    names = [
        "length", "log_length", "gc", "frac_A", "frac_C", "frac_G", "frac_T",
        "purine_frac", "g4h_mean", "g4h_absmean", "g4h_win_max",
        f"n_{base}runs", f"max_{base}run", f"mean_{base}run",
        "n_Gruns", "max_Grun", "n_Cruns", "max_Crun",
        "has_4tracts", "tract_len_min", "tract_len_max", "tract_len_mean",
        "loop1", "loop2", "loop3", "loop_total", "loop_mean", "loop_sd",
        "loop_min", "loop_max", "loop_asymmetry", "core_span", "flank_len",
        "loop_gc", "loop_frac_T", "loop_frac_A",
        "n_bulges", "longest_homopolymer", "shannon_entropy",
    ]
    names += [f"k1_{x}" for x in _K1]
    names += [f"k2_{x}" for x in _K2]
    names += [f"k3_{x}" for x in _K3]
    return names


def _base_features(seq: str, kind: str = "G4") -> list[float]:
    s = clean(seq)
    n = len(s)
    if n == 0:
        return [0.0] * len(_base_feature_names(kind))
    base = "G" if kind == "G4" else "C"

    counts = {b: s.count(b) for b in "ACGT"}
    frac = {b: counts[b] / n for b in "ACGT"}
    gc = frac["G"] + frac["C"]
    purine = frac["A"] + frac["G"]

    g4h = g4hunter_mean(s)
    g4h_abs = abs(g4h)
    g4h_win = g4hunter_window_max(s, window=min(25, n))

    nb, mb, avb = _run_stats(s, base)
    ng, mg, _ = _run_stats(s, "G")
    nc, mc, _ = _run_stats(s, "C")

    parsed = parse_tracts(s, base, min_tract=3, max_loop=12)
    if parsed:
        tracts, loops = parsed
        tract_lens = [b - a for a, b in tracts]
        l1, l2, l3 = (list(loops) + [0, 0, 0])[:3]
        core_start, core_end = tracts[0][0], tracts[-1][1]
        loop_seq = "".join(s[tracts[i][1] : tracts[i + 1][0]] for i in range(3))
        has4 = 1.0
    else:
        tract_lens = [0]
        l1 = l2 = l3 = 0
        core_start, core_end = 0, n
        loop_seq = ""
        has4 = 0.0

    loops3 = [l1, l2, l3]
    loop_total = sum(loops3)
    loop_mean = loop_total / 3.0
    loop_sd = float(np.std(loops3))
    loop_asym = abs(l1 - l3) / max(l1 + l3, 1)
    ln = max(len(loop_seq), 1)
    loop_gc = (loop_seq.count("G") + loop_seq.count("C")) / ln if loop_seq else 0.0
    loop_T = loop_seq.count("T") / ln if loop_seq else 0.0
    loop_A = loop_seq.count("A") / ln if loop_seq else 0.0

    # bulges: tracts longer than the minimum, i.e. non-uniform tract lengths
    n_bulges = float(sum(1 for t in tract_lens if t != min(tract_lens)))

    longest_hp = max(
        (_run_stats(s, b)[1] for b in "ACGT"), default=0
    )
    probs = [frac[b] for b in "ACGT" if frac[b] > 0]
    entropy = -sum(p * math.log2(p) for p in probs)

    feats = [
        float(n), math.log(n), gc, frac["A"], frac["C"], frac["G"], frac["T"],
        purine, g4h, g4h_abs, g4h_win,
        float(nb), float(mb), avb,
        float(ng), float(mg), float(nc), float(mc),
        has4, float(min(tract_lens)), float(max(tract_lens)), float(np.mean(tract_lens)),
        float(l1), float(l2), float(l3), float(loop_total), loop_mean, loop_sd,
        float(min(loops3)), float(max(loops3)), loop_asym,
        float(core_end - core_start), float(n - (core_end - core_start)),
        loop_gc, loop_T, loop_A,
        n_bulges, float(longest_hp), entropy,
    ]
    feats += _kmer_freq(s, _K1, 1)
    feats += _kmer_freq(s, _K2, 2)
    feats += _kmer_freq(s, _K3, 3)
    return feats


# --------------------------------------------------------------------------- #
# Condition block
# --------------------------------------------------------------------------- #
CONDITION_FEATURE_NAMES = [
    "k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct", "strand_conc",
    "log1p_k", "log1p_na", "log1p_li_nh4", "log1p_mg",
    "monovalent", "log1p_monovalent", "ionic_strength", "log1p_ionic_strength",
    "k_fraction", "na_fraction", "li_fraction",
    "g4_competent_cation", "log1p_g4_competent",
    "ph_minus_7", "pka_minus_ph", "protonation_frac",
    "log_strand_conc", "has_crowder", "no_supporting_cation",
]


# --------------------------------------------------------------------------- #
# The 'locus' feature set
# --------------------------------------------------------------------------- #
# A locus is a duplex position, and the question asked of it -- does this site
# support a G-quadruplex, an i-motif, both, or neither -- needs BOTH tract
# systems in the feature vector. The G4 and iM feature sets share their
# composition and k-mer block and differ only in which base the tract/loop
# architecture is parsed on, so the locus set is the G4 vector plus the i-motif
# architecture columns, rather than two nearly-identical vectors concatenated.
_ARCH_SUFFIXES = (
    "runs", "run", "has_4tracts", "tract_len_min", "tract_len_max",
    "tract_len_mean", "loop1", "loop2", "loop3", "loop_total", "loop_mean",
    "loop_sd", "loop_min", "loop_max", "loop_asymmetry", "core_span",
    "flank_len", "loop_gc", "loop_frac_T", "loop_frac_A", "n_bulges",
)


_IM_ARCH_NAMES = (
    "n_Cruns", "max_Crun", "mean_Crun", "has_4tracts",
    "tract_len_min", "tract_len_max", "tract_len_mean",
    "loop1", "loop2", "loop3", "loop_total", "loop_mean", "loop_sd",
    "loop_min", "loop_max", "loop_asymmetry", "core_span", "flank_len",
    "loop_gc", "loop_frac_T", "loop_frac_A", "n_bulges",
)


def _im_architecture_index() -> list[int]:
    """Positions in the iM feature vector that describe C-tract architecture.

    Deduplicated by name: the C-run columns appear twice in the base list --
    once as the kind-parameterised ``n_{base}runs`` and once as the explicit
    ``n_Cruns`` -- and carrying the same number twice would just give the
    ensemble two identical columns to split on.
    """
    names = _base_feature_names("iM")
    keep, seen = [], set()
    for i, nm in enumerate(names):
        if nm in _IM_ARCH_NAMES and nm not in seen:
            seen.add(nm)
            keep.append(i)
    return keep




def sequence_feature_names(kind: str) -> list[str]:
    if kind != "locus":
        return _base_feature_names(kind)
    im = _base_feature_names("iM")
    return _base_feature_names("G4") + [f"im_{im[i]}" for i in _im_architecture_index()]


def sequence_features(seq: str, kind: str = "G4") -> list[float]:
    if kind != "locus":
        return _base_features(seq, kind)
    g = _base_features(seq, "G4")
    im = _base_features(seq, "iM")
    return g + [im[i] for i in _im_architecture_index()]


def condition_features(c: Condition) -> list[float]:
    mono = c.monovalent
    prot = 1.0 / (1.0 + 10.0 ** (c.ph - CYTOSINE_N3_PKA))  # fraction of C protonated
    return [
        c.k, c.na, c.li_nh4, c.mg, c.ph, c.temperature, c.crowder_pct, c.strand_conc,
        math.log1p(c.k), math.log1p(c.na), math.log1p(c.li_nh4), math.log1p(c.mg),
        mono, math.log1p(mono), c.ionic_strength, math.log1p(c.ionic_strength),
        c.k_fraction,
        (c.na / mono) if mono else 0.0,
        (c.li_nh4 / mono) if mono else 0.0,
        c.g4_competent_cation, math.log1p(c.g4_competent_cation),
        c.ph - 7.0, CYTOSINE_N3_PKA - c.ph, prot,
        math.log1p(c.strand_conc), 1.0 if c.crowder_pct > 0 else 0.0,
        1.0 if c.g4_competent_cation < 1.0 else 0.0,
    ]


# --------------------------------------------------------------------------- #
# Interaction block
# --------------------------------------------------------------------------- #
INTERACTION_FEATURE_NAMES = [
    "gtracts_x_logK", "gc_x_logIonic", "loopsum_x_logK", "len_x_logK",
    "ctracts_x_protonation", "cfrac_x_protonation", "loopsum_x_protonation",
    "gc_x_temp", "g4h_x_logKcompetent", "tracts_x_pH",
]


def interaction_features(seqfeat: Sequence[float], c: Condition, kind: str) -> list[float]:
    # the shared block is at the front for every kind, so G4 names index it
    names = sequence_feature_names("G4" if kind == "locus" else kind)
    ix = {nm: i for i, nm in enumerate(names)}
    logK = math.log1p(c.k)
    logI = math.log1p(c.ionic_strength)
    logKc = math.log1p(c.g4_competent_cation)
    prot = 1.0 / (1.0 + 10.0 ** (c.ph - CYTOSINE_N3_PKA))
    ng = seqfeat[ix["n_Gruns"]]
    nc = seqfeat[ix["n_Cruns"]]
    gc = seqfeat[ix["gc"]]
    loops = seqfeat[ix["loop_total"]]
    length = seqfeat[ix["length"]]
    cfrac = seqfeat[ix["frac_C"]]
    g4h = seqfeat[ix["g4h_mean"]]
    tracts = seqfeat[ix["has_4tracts"]]
    return [
        ng * logK, gc * logI, loops * logK, length * logK,
        nc * prot, cfrac * prot, loops * prot,
        gc * c.temperature / 100.0, g4h * logKc, tracts * c.ph,
    ]


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def feature_names(kind: str = "G4", *, use_conditions: bool = True) -> list[str]:
    names = list(sequence_feature_names(kind))
    if use_conditions:
        names += list(CONDITION_FEATURE_NAMES) + list(INTERACTION_FEATURE_NAMES)
    return names


def featurize(
    sequences: Iterable[str],
    conditions: Iterable[Condition] | None = None,
    *,
    kind: str = "G4",
    use_conditions: bool = True,
) -> np.ndarray:
    """Build the design matrix for a batch."""
    sequences = list(sequences)
    if use_conditions:
        conds = list(conditions) if conditions is not None else [Condition()] * len(sequences)
        if len(conds) != len(sequences):
            raise ValueError("sequences and conditions must be the same length")
    rows: list[list[float]] = []
    for i, s in enumerate(sequences):
        sf = sequence_features(s, kind)
        if use_conditions:
            c = conds[i]
            rows.append(sf + condition_features(c) + interaction_features(sf, c, kind))
        else:
            rows.append(sf)
    return np.asarray(rows, dtype=np.float32)
