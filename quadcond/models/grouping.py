"""Sequence grouping, so cross-validation does not lie to us.

G4/i-motif datasets are full of near-duplicates: the same telomeric repeat with
one extra flanking base, mutational series around one promoter, the same PDB
sequence deposited twice.  A random CV split puts siblings on both sides of the
fold and reports an accuracy the model will not reproduce on anything new.

``cluster_sequences`` performs greedy single-pass clustering on k-mer cosine
similarity and returns a group label per sequence, for use with
``GroupKFold``/``StratifiedGroupKFold``.  Every reported metric in QuadCond is
grouped; the ungrouped number is also reported so the size of the illusion is
visible.
"""
from __future__ import annotations

import itertools

import numpy as np

from ..motifs import clean


def kmer_matrix(sequences, k: int = 4) -> np.ndarray:
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    idx = {km: i for i, km in enumerate(kmers)}
    M = np.zeros((len(sequences), len(kmers)), dtype=np.float32)
    for r, s in enumerate(sequences):
        s = clean(s)
        for i in range(len(s) - k + 1):
            j = idx.get(s[i : i + k])
            if j is not None:
                M[r, j] += 1
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return M / norms


def cluster_sequences(
    sequences,
    threshold: float = 0.90,
    k: int = 4,
    *,
    exact_max: int = 6000,
) -> np.ndarray:
    """Group near-identical sequences; returns an integer group per sequence.

    Exact greedy clustering is O(n^2) and is used up to ``exact_max`` sequences.
    Above that it becomes the dominant cost, so we fall back to MiniBatchKMeans
    on the same k-mer profiles with a cluster count chosen to give a comparable
    mean group size.  The approximation only affects how conservative the folds
    are, never the labels, and the realised group statistics are reported in
    ``training_meta.grouping`` either way.
    """
    seqs = list(sequences)
    if len(seqs) > exact_max:
        from sklearn.cluster import MiniBatchKMeans

        M = kmer_matrix(seqs, k=k)
        n_clusters = max(2, min(len(seqs) // 4, 4000))
        km = MiniBatchKMeans(n_clusters=n_clusters, random_state=0, n_init=3,
                             batch_size=2048)
        return km.fit_predict(M).astype(int)
    return _greedy_clusters(seqs, threshold, k)


def _greedy_clusters(sequences, threshold: float, k: int) -> np.ndarray:
    seqs = list(sequences)
    if not seqs:
        return np.zeros(0, dtype=int)
    M = kmer_matrix(seqs, k=k)
    order = np.argsort([-len(clean(s)) for s in seqs])
    groups = -np.ones(len(seqs), dtype=int)
    centroids: list[int] = []
    for i in order:
        if groups[i] >= 0:
            continue
        gid = len(centroids)
        centroids.append(i)
        groups[i] = gid
        sims = M[order] @ M[i]
        for pos, j in enumerate(order):
            if groups[j] < 0 and sims[pos] >= threshold:
                groups[j] = gid
    return groups


def group_report(groups: np.ndarray) -> dict:
    _, counts = np.unique(groups, return_counts=True)
    return {
        "n_items": int(len(groups)),
        "n_groups": int(len(counts)),
        "largest_group": int(counts.max()) if len(counts) else 0,
        "singleton_groups": int((counts == 1).sum()),
        "mean_group_size": float(counts.mean()) if len(counts) else 0.0,
    }
