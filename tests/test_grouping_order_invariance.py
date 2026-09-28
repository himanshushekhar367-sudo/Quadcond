"""The sequence grouping must not depend on the order rows arrive in.

A grouped cross-validation partition is only reproducible if the grouping is a
function of the data alone. It was not: `_greedy_clusters` ordered sequences with
an unstable sort over their lengths, and in the G4 melting collection almost every
sequence shares a length with another, so tied rows were traversed in an arbitrary
order and a greedy pass seeded different groups. Shuffling one real input moved
the group count between 1,201 and 1,203, and the folds behind the shipped v0.5.1
heads cannot be replayed exactly as a result.

These tests fail if that regresses.
"""
from __future__ import annotations

import random

import numpy as np

from quadcond.models.grouping import cluster_sequences


def _partition(labels, seqs):
    """The grouping as a set of frozensets of sequences, independent of labelling.

    Comparing label arrays directly would fail on a harmless renumbering. What
    matters is which sequences share a group.
    """
    out = {}
    for lab, s in zip(labels, seqs):
        out.setdefault(int(lab), set()).add(s)
    return frozenset(frozenset(v) for v in out.values())


def _corpus(n=300, seed=0):
    """Sequences that mostly share lengths, which is the condition that broke it."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        # Same length for most of them, so the sort keys tie.
        length = rng.choice([21, 21, 21, 22, 22, 27])
        out.append("".join(rng.choice("ACGT") for _ in range(length)))
    return out


def test_grouping_is_repeatable_on_identical_input():
    seqs = _corpus()
    first = cluster_sequences(seqs)
    second = cluster_sequences(seqs)
    assert _partition(first, seqs) == _partition(second, seqs)


def test_grouping_is_invariant_to_input_order():
    seqs = _corpus()
    reference = _partition(cluster_sequences(seqs), seqs)
    rng = np.random.default_rng(0)
    for _ in range(5):
        idx = rng.permutation(len(seqs))
        shuffled = [seqs[i] for i in idx]
        assert _partition(cluster_sequences(shuffled), shuffled) == reference, (
            "shuffling the input changed which sequences share a group; the "
            "partition is order-dependent again and no published fold assignment "
            "can be trusted")


def test_group_count_is_invariant_to_input_order():
    seqs = _corpus(n=400, seed=1)
    n_ref = len(set(cluster_sequences(seqs).tolist()))
    rng = np.random.default_rng(1)
    for _ in range(5):
        idx = rng.permutation(len(seqs))
        n = len(set(cluster_sequences([seqs[i] for i in idx]).tolist()))
        assert n == n_ref, f"group count moved from {n_ref} to {n} under a shuffle"


def test_identical_sequences_share_a_group():
    seqs = ["GGGTTAGGGTTAGGGTTAGGG"] * 3 + ["ATCATGCATGACTAGCATCGA"]
    labels = cluster_sequences(seqs)
    assert labels[0] == labels[1] == labels[2]


if __name__ == "__main__":
    # pytest is not in every environment here, and this check is worth being able
    # to run anywhere quadcond imports.
    import sys
    failures = []
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures.append((name, exc))
                print(f"FAIL {name}: {exc}")
    print()
    print(f"{'FAILED' if failures else 'OK'}: {len(failures)} failing")
    sys.exit(1 if failures else 0)
