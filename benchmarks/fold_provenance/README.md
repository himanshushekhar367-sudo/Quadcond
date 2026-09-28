# Fold provenance

Two scripts that answer one question: **was the cross-validation partition used to
train the shipped v0.5.1 models ever recoverable, and if not, what follows?**

    python compare_sidecar_vs_manifest.py     # no heavy dependencies
    python check_determinism.py               # needs numpy; run in qcgw or .venv

## Why this matters

Folds were built at training time and discarded. The model sidecars record how
many groups a head's rows fell into, but not which row went where.
`quadcond/models/foldmanifest.py` reconstructs an assignment by replaying
`build_folds` on rows read today, and labels itself `reconstructed` rather than
`historical` — correctly, because a matching row fingerprint shows only that the
rows are the same, not that the split is.

## What the comparison found

`results/sidecar_vs_manifest.json`, regenerable at any time:

| Head | group_by | rows (record / replay) | groups (record / replay) | reproduced |
| --- | --- | --- | --- | --- |
| g4_fold | sequence | 6,744 / 6,744 | 1,619 / 1,624 | no |
| g4_topology | sequence | 1,005 / 1,005 | 617 / 618 | no |
| g4_tm | sequence | 2,274 / 2,274 | 400 / 407 | no |
| im_fold | sequence | 320 / 320 | 175 / 171 | no |
| im_pht | sequence | 160 / 160 | 85 / 85 | **yes** |
| im_pht_condition | condition | 558 / 558 | 10 / 10 | yes |
| im_tm_condition | condition | 379 / 379 | 10 / 10 | yes |

Row counts match everywhere, so the rows being read are the right rows and the
disagreement is in the grouping alone. Four of the five sequence-grouped heads
fail. `im_pht` — the smallest at 160 rows — reproduces, and both condition-grouped
heads reproduce exactly, which is what you would expect if the mechanism is
tie-breaking: fewer rows and fewer tied lengths mean fewer chances to diverge.

## The suspected mechanism

`_greedy_clusters` orders sequences with `np.argsort([-len(clean(s)) for s in seqs])`.
NumPy's default sort is not stable, so sequences of equal length come out in an
arbitrary order. The clustering is greedy and seeds a new group from whichever
tied sequence it reaches first, so a different traversal produces a different
number of groups from identical input.

### Confirmed by measurement

`check_determinism.py`, run on 2,274 distinct G4 sequences from `atlas_core.db`
(`results/determinism_check.json`):

| Check | Result |
| --- | --- |
| Same input, twice | 1,202 and 1,202 groups — deterministic within a process |
| Same sequences, shuffled three times | 1,201 / 1,203 / 1,203 groups — **order-dependent** |
| Sequences sharing a length with another | 2,210 of 2,274 |

So the clustering is a pure function of *(sequence list, order)*, and the order
changes the answer. With almost every sequence tied on length, an unstable sort
had ample room to permute them.

Two cautions on reading those numbers. The 1,202 is not comparable with the 400 or
407 in the table above: this check reads 2,274 *distinct* sequences, whereas the
`g4_tm` head has 2,274 *rows* over a smaller and less diverse sequence set. And
the shuffle moved the count by 1–2 groups in ~1,200, while the record-versus-replay
gap for `g4_tm` is 7 in 400 — proportionally larger. Order-dependence is therefore
demonstrated as *a* cause; whether it is the *whole* cause of the 400 → 407 shift
is not established, and the clustering code or the atlas row set may also have
changed since training.

### The fix, applied

`_greedy_clusters` now orders sequences by `(-length, sequence)`, making the
traversal a pure function of the sequence content. A stable sort alone would only
have fixed run-to-run variation for one fixed input order, which is weaker than
what a published partition needs. `tests/test_grouping_order_invariance.py` fails
if this regresses, checking the partition as a set of groups rather than as a label
array so a harmless renumbering does not trip it.

**This changes fold assignments**, so it belongs to a deliberate v0.5.2 step: the
existing fold manifest must be regenerated under the new grouping, and any metric
recomputed under it will differ slightly from both the shipped record and the
earlier replay. That is the point — from here the partition is reproducible.

To verify after the change:

    conda activate qcgw
    python benchmarks/fold_provenance/check_determinism.py   # order_dependent should be false
    python -m pytest tests/test_grouping_order_invariance.py -q

## What follows for the paper

- Table 1's per-head metrics were computed under a partition that no longer exists
  and cannot be regenerated. They remain a truthful record of the training run;
  they are not a claim a third party can verify.
- Any metric computed under the reconstructed manifest is internally sound but is
  **not** the shipped artifact's historical out-of-fold figure. The two must not be
  presented as the same number.
- The head-to-head benchmark is unaffected: all tools shared whichever partition
  was used, and the claim is relative ordering.

## The fix

1. Break the ties deterministically — `kind="stable"` plus a canonical row sort —
   and confirm two consecutive replays agree.
2. Publish the full atlas with a fixed row order and its fingerprint.
3. Ship the fold manifest as a versioned artifact with its own SHA-256 from
   v0.5.2, and recompute reported metrics under it.
4. State in the manuscript that the released artifact was trained under a
   partition that could not be reconstructed exactly, and that reported metrics
   are computed under the published one.
