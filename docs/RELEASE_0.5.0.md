# Release 0.5.0 — the asset decision

This release records a model artifact that is **not** the one the 0.4.8 asset
manifest expected. That is a deliberate, documented substitution, and this file
is the record of what was established before it was made.

## The three numbers

| | version |
| --- | --- |
| `quadcond` package | 0.5.0 |
| AENNA-3D workbench | 0.5.0 |
| asset manifest | 0.5.0 |
| **model artifact** | **0.4.6** |

The first three are one number because the repository is now one product. The
fourth is deliberately different: the code moved and the estimators did not.
A single number would have to misdescribe one of them, and it would misdescribe
the one that produces the numbers. Every prediction response now carries
`application_version`, `model_version_reported_by_artifact` and a three-valued
`model_version_matches_manifest`, so a consumer cannot see one and assume the
other.

## What was expected, and what is here

| asset | expected by 0.4.8 | recorded by 0.5.0 |
| --- | --- | --- |
| `model` | `204b4608…`, 27,280,002 B | `49364019…`, 27,280,001 B |
| `model_seqonly` | `b4fca5d3…`, 25,729,065 B | `402e9db1…`, **same byte count** |
| `atlas_core` | `145711cd…` | `145711cd…` — unchanged, verified |

The expected hashes are not deleted. They are kept under each asset's
`supersedes` in `assets_manifest.json`, because a checksum that is simply
overwritten leaves no trace that it was ever expected, and that is precisely
the history an auditor needs.

`atlas_core` verifying unchanged is worth noting: the 0.4.8 manifest is not
wholesale stale. It correctly describes the atlas on disk. The model entry was
its single discrepancy.

## What was established

Both artifacts were loaded in the environment named in `built_with` — an exact
match on scikit-learn 1.8.0, xgboost 3.2.0, numpy 2.4.4, scipy 1.17.1 and
joblib 1.5.3 — and compared against `MODEL_CARD.md`:

- full model: **12 of 12 heads** agree on headline metric and training-row count
- sequence-only model: **12 of 12 heads** agree with the ablation column,
  including its most distinctive entries (`im_tm_condition` r² 0.8532 → 0.0634,
  `g4_tm` 0.6701 → 0.2542)

**24 of 24 comparisons, no mismatches.** Live predictions compute, and the
refusal paths return their reasons in full.

## What was NOT established

1. **These metrics were read back, not re-derived.** They are stored inside the
   artifact. Matching them establishes agreement with the model card's stored
   metrics, not unique artifact identity or estimator equivalence. Re-deriving them needs the 354 MB `atlas` asset and
   `scripts/02_train.py`; that atlas is not present, so this reproduction is
   still open.
2. **Whether the superseded hashes describe a different model or a different
   save of this one.** The payload is compressed, so the version stamp cannot
   be read from the raw bytes. The seq-only artifact having the *exact*
   expected byte count with a different hash is suggestive of a same-shape
   re-save, but suggestive is not proof and nothing is built on it.
3. **That the expected binary does not exist.** A search of `C:\Users`, `D:\`
   and `E:\` — excluding system directories, and not reading inside archives —
   did not find it. The four archives in the project folder were checked
   separately and contain no `.joblib`. Absence from the locations searched is
   not absence.

If that binary is recovered, the honest action is to restore it as the `model`
asset and demote this one, not to keep both.

## Reproducing this release

```bash
quadcond assets status          # all required assets ok
python -m quadcond.service --port 8765
```

`artifacts/quadcond_model.joblib`, `artifacts/quadcond_model_seqonly.joblib`
and `data/atlas_core.db` are resolved from the repository directory and are
gitignored. `quadcond assets fetch` cannot retrieve them until the `url` field
of each asset is filled in against a published release.
