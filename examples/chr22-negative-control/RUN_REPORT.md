# End-to-end run: chr22:36,201,668–36,201,728

*Revised 15 September 2026 after review. Corrections are listed at the end.*

## What this run establishes

The export, the variant matching and the missing-value handling work correctly.
That is the claim; it is narrower than "the integration is validated", and this
run does not exercise a useful structural-effect prediction or say anything
about accuracy.

| | |
|---|---|
| Variants joined | **183 / 183**, exact coordinates and alleles |
| AVI values | every one matches the original Atlas export |
| Variants with a structural delta | **0** |
| Quadrant counts | 183 unclassified; 0 in every other cell |
| Structural candidates (gain/loss/count change) | 0 |
| Model | `0.4.6`, sha256 `493640193ad7…`, **does not match the release manifest** |
| Libraries | scikit-learn 1.8.0, xgboost 3.2.0 — the manifest's `built_with` |
| Condition | `K140 Na12 Mg1 pH7.4 37C`, stated by the caller |

## Why every structural cell is empty

**QuadCond did not identify an applicable canonical motif under its implemented
rules, on either strand, and withheld structural deltas.** G4Hunter mean is
0.262, below the 0.5 gate.

That is a model-applicability result. It is **not** an experimental
demonstration that this DNA cannot form a non-canonical structure, and the
earlier phrasing — "there is no structure here to perturb" — was a claim about
the molecule that the output does not support.

Keep this window as a **negative-control software case**: it tests that
regulatory scores stay available while unsupported structural predictions stay
missing, which is precisely what it demonstrated.

## The two heads are not equally condition-aware

The run record now reports this per head, per axis, from each head's own
training ranges:

| head | responds to | inert for | predicts |
|---|---|---|---|
| `g4_tm` | K⁺, Na⁺, Li⁺/NH₄⁺, Mg²⁺, pH, strand conc. | crowder | temperature |
| `im_pht` | *nothing* | every axis | — |

`im_pht` was trained on 160 rows in a single buffer. It has **no learned
response to any condition**, so the buffer set above is inert for it and its
value is a reference-condition prediction. `g4_tm` *predicts* temperature rather
than consuming it — the requested 37 °C is the evaluation point of the
downstream two-state model, not an input to the Tm prediction.

A "physiological preset" is a stated set of numbers for reproducibility. It is
not an intracellular folding prediction.

## Provenance, now recorded

The previous outputs carried `"model_artifact_sha256": ""` — a field satisfied
and nothing identified. A path names a location in a filesystem, not a set of
bytes. Every run record now carries:

- the SHA-256 of the model file actually opened, the version it reports, the
  manifest's expectation, and a three-valued `matches_manifest`;
- Python, platform and per-package versions, compared package-by-package against
  the manifest's `built_with`;
- the regulatory table's own checksum, size and modification time, plus the
  exporter's `provenance.json` where one sits beside it;
- the coordinate convention (1-based inclusive), the strand of the window
  sequence, and the assembly as stated by the caller.

The CSV is now a clean interchange file: column names on line 1, nothing before
them. The narrative and provenance live in
`variant_join_chr22_36201668.provenance.json`.

## Screening now covers motifs a variant would create

The screen previously looked at the reference sequence only. That silently
discards the most interesting case a variant workflow can find — a window with
no motif where one substitution creates one — and would have declared such a
window uninteresting and never queried it.

It now screens the reference **and every single-base alternate, on both
strands**, and reports four states rather than two: `motif_present`,
`motif_gain_possible`, `motif_loss_possible`, `no_motif`. On a deliberately
near-miss G4 (`GGGTTAGGGTTAGGGTTAGTG`) it returns `motif_gain_possible` with the
creating substitution named — the case the old screen would have thrown away.

## Motif-gain and motif-loss variants no longer fall off the shortlist

The combined ranking requires a delta on both axes. Motif loss and gain have no
delta by construction, so they landed in `unclassified` and sorted off the
bottom of any ranked table — a retrieval failure, not a formatting one.

They now get a separate categorical list, `structural_candidates`, ranked within
category by the regulatory axis where one exists, and rendered in AENNA as their
own view. They are **not** given an invented numerical effect.

## On the Pu27 numbers

Same binary and libraries, run on MYC Pu27 as a sequence with no coordinate
asserted: predicted wild-type `g4_tm` 77.0 °C, and the largest predicted
destabilisations at G11C, G20C, G9A, G9C.

These are **model predictions, not measured temperatures**. The
guanine-in-G-tract pattern is a sanity check on plausibility. Neither that
pattern nor the paired-estimator spread validates mutation-effect accuracy, and
neither has been independently reproduced.

## Corrections to the previous version of this report

1. "Validates the plumbing, completely" → the export, matching and missing-value
   handling work; the integration is not validated.
2. "There is no structure here to perturb" → QuadCond found no applicable
   canonical motif under its rules and withheld deltas.
3. "That first run cost 183 API calls" → **wrong**. The exporter uses
   `query_interval()`, which batches the region and paginates; the 183 is a
   count of variants, not requests. Correct statement: *the query returned
   scores for 183 variants.* Request count has not been measured.
4. Model identity: see `MODEL_IDENTITY.md`, which withdraws "the only difference
   is the version string".
5. Both heads described as running under the stated buffer → `im_pht` is inert
   for every condition axis.

## Next, in order

1. Resolve artifact identity and freeze a reproducible exploratory environment,
   even while release sign-off stays pending.
2. Keep chr22 as the negative control.
3. Pick a literature-supported motif-bearing locus. MYC Pu27 is a reasonable
   demonstration candidate, but its exact sequence, GRCh38 mapping and strand
   must be verified before use — not taken from this document.
4. Run the full comparison there, including retained, gained and lost motifs.
   Non-zero deltas are not a success condition.
5. Verify AENNA end to end, including that browser exports reproduce backend
   values, and that missing-value labels survive the round trip.
6. Add one tissue-specific regulatory output for a defined question, preserving
   signed scores, gene and track identifiers and scorer configuration rather
   than collapsing to a maximum.
7. Evaluate scientific benefit on independent data: Atlas alone vs QuadCond
   alone vs combined, locus-separated, with a training-overlap audit and a
   simple motif-disruption baseline.
