# Step 3 verification, 15 September 2026

The current backend and workbench were reviewed and tested locally. The existing
step-3 edits were preserved. This report records software verification, not
independent scientific validation of the trained model.

## Changes verified

- Scan run records include per-head condition metadata for the ordinary backend
  workflows. The interface distinguishes varied, fixed/unused, predicted target,
  and unknown fields. Variation in training inputs is not presented as evidence
  that the fitted model learned a response.
- Categorical candidates distinguish motif loss, gain, count change, and gain on
  the other strand. They carry the selected structural head, strand, units,
  explanation, delta, and explicit delta availability.
- The scan detects a motif gain on the strand not selected for numerical
  prediction and retains it as an unscored candidate. It does not numerically
  score both strands. Reversing the input alone is not a reliable workaround
  for the wild-type-based strand selection.
- The production HTTP variant route was tested on an ephemeral local port with
  a software fixture. The T32C reverse-strand candidate and condition metadata
  survive JSON serialization. These are software test inputs, not an asserted
  genomic locus or experimentally verified structure.
- The workbench in `web/dist` was rebuilt. An already-running Python service
  needs restarting to import the changed backend modules.

## Checks

- Python: 112 passed, 38 skipped. Skips include unavailable model/data assets,
  supporting source files, and optional pysam. This is not an all-assets pass.
- Frontend: 51 passed.
- TypeScript typecheck and Vite production build: passed. Vite reports a large
  visualization chunk; this is a bundle-size warning, not a failed build.
- Ruff on the changed Python files and `git diff --check`: passed.
- Asset-search regression tests: 3 passed, including distinct files with equal
  names and sizes, overlapping roots, quick-mode semantics, and invalid roots.

## Asset search and current server

The search script beside the repository was repaired to avoid argparse variable
shadowing and incorrect deduplication by filename and size. It deduplicates paths
and verifies candidate contents against the unchanged manifest.

The default search examined 25 candidate files under `C:\Users`, `D:\`, and
`E:\`. It excludes the directories listed in its output, files smaller than
1 MB, unrecognized extensions, and archive contents. It is not exhaustive.

Only `sandbox/atlas_core.db` matched a recorded identity. No expected model,
seq-only model, or full atlas matched in that scan. The seq-only reassembled
artifact matched the expected length but not its recorded identity.

Full output: `C:\Users\pc\Downloads\Quadcond\asset-search-20260915.txt`.

The running localhost:8765 health response separately reported `atlas_core`
as `ok`, `ready: false`, and `missing_required_assets: ["model"]`.
No model was repackaged, retrained, or relabeled; the manifest and its verification
gate were not changed. No claim of fitted-estimator equivalence follows from
matching stored metadata. Real-model release validation still requires recovery
of the expected artifact or a documented new release of the available artifact
with appropriate reproducibility checks.
