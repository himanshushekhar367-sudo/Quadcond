# Where everything lives

A map of the project's artefacts, written because several of them are gitignored,
unpublished, or only meaningful alongside a provenance file. Paths are relative to
the repository root.

## Tracked in git

| Path | What it is |
| --- | --- |
| `quadcond/` | the package: heads, conditions, claims, motifs, genome scanner, variant workflow, Atlas client, service |
| `web/` | AENNA-3D workbench (React, TypeScript, Three.js) |
| `manuscript/` | v0.5.1 draft, builder, frozen inputs, `verify_numbers.py` |
| `release-v0.4.9/manuscript/` | the previous draft, kept unchanged for the record |
| `benchmarks/published_tools/` | head-to-head benchmark: scripts, result tables, figures |
| `benchmarks/fold_provenance/` | whether the training partition is recoverable, and what follows |
| `benchmarks/ingest_audit/` | overlap audits run before a candidate dataset is allowed near training |
| `benchmarks/external_validation/` | panels published after the artifact was frozen, and the harness that scores them |
| `genomewide/` | the genome-wide pipeline and its README, which carries the pilot result |
| `docs/` | release notes, verification records, `PROJECT_RECORD.md`, this file |
| `tests/` | including `test_fold_manifest.py` and the variant-join tests |

## On disk, not in git

| Path | Size | Why it is not tracked, and what to do |
| --- | --- | --- |
| `data/atlas.db` | ~354 MB | the full 398,375-record atlas the shipped models were trained on. **Unpublished, and this is a blocker** — the manifest records its fingerprint but has no release URL. Publish it to Zenodo with a fixed row order. |
| `data/atlas_core.db` | ~5.4 MB | the trimmed core that ships with the release. Do not use it for overlap audits: it is a subset, so overlap computed against it is a lower bound. |
| `artifacts/*.joblib` | ~27 MB, ~26 MB, ~1.8 MB | main model, sequence-only ablation, G4-seq scanner. Hashes are in `quadcond/assets_manifest.json`; attach them to the v0.5.1 GitHub release. |
| `gw_out/` | ~568 MB | genome-wide run output: per-chromosome motif, structural and AVI tables, plus `stats/` and `stats_noflank/`. Frozen copies of the two summary files the manuscript uses live in `manuscript/inputs/`. |
| `.venv/` | — | Windows model environment, Python 3.11, versions pinned in `requirements-model.txt`. The `qcgw` conda environment in WSL works too. |

## Provenance files worth knowing about

- `quadcond/assets_manifest.json` — SHA-256 of every shipped asset, verified before
  the file is opened; a mismatch is an error, not a warning.
- `benchmarks/published_tools/results/fold_assignments_meta.json` — the replayed
  fold manifest. Its `assignment_provenance` field says `reconstructed`, which
  means it must not be used to clear a head of leakage.
- `benchmarks/fold_provenance/results/sidecar_vs_manifest.json` — the measured
  comparison between the training record and that replay.
- `manuscript/inputs/` — frozen copies of every result file the manuscript quotes,
  so the document builds from this directory alone and cannot drift from the
  analysis.
- `manuscript/verify_numbers.py` — recomputes the manuscript's quantitative claims
  from those inputs and fails on any mismatch. Run it after any edit touching a
  number.

## Environments

The model artifact loads only where numpy 2.4.4, scipy 1.17.1, scikit-learn 1.8.0
and xgboost 3.2.0 are installed with `quadcond` alongside them: the `qcgw` conda
environment in WSL, or `.venv\Scripts\python.exe` on Windows. Conda `base` is not
one of them, and the failure looks like a missing numpy rather than a wrong
interpreter.

## Reproducing the paper

    cd manuscript
    python build_manuscript.py      # renders the .md and .docx; fits nothing
    python verify_numbers.py        # checks every quantitative claim

Neither step needs the atlas or the model, because the manuscript builds from
frozen result files. Reproducing the *results* does need both, and the atlas is
not yet published — see `benchmarks/fold_provenance/README.md` for what that
means and what to do about it.
