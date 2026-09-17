# Local real-model workbench verification

Verified 16 September 2026 on Windows, starting from commit `e6e0b71`.
The additional changes described below are local working-tree changes.

## Outcome

The application now runs at http://localhost:8765/ using the real model stamped
0.4.6, verified against the documented 0.5.0 manifest, and the verified core
atlas. The saved chr22 AlphaGenome score table is attached. No new AlphaGenome
query or key was needed. No model was retrained or modified during this check.

The background server was started from this checkout's `.venv` with
`bridge/serve_table.py`. Logs are in the parent folder:
`workbench-20260916.stdout.log` and `workbench-20260916.stderr.log`.

## Problems found and corrected

1. The machine's default Python had scikit-learn 1.6.1. The fitted calibrator
   raised `AttributeError: LogisticRegression has no attribute multi_class`.
   This reproduced as ten failing tests. A separate `.venv` now matches the
   five recorded model dependencies in `requirements-model.txt`. Global Python
   packages were not changed. The existing test suite also required openpyxl
   and jsonschema, now declared in the development extra.
2. Both bridge copies rejected application 0.5.0 and the runner defaulted to
   the old release folder. The inspected versions 0.4.9 and 0.5.0 are accepted,
   and the default root is the current repository. Unknown versions still fail.
3. The Variant join panel filtered away every unranked row. For chr22 this hid
   all 183 scored variants. It now shows coverage, the unclassified count,
   ranked and unranked rows, and an option to display all variants. Missing
   deltas and ranks remain missing. Its unsupported assembly-verification claim
   was removed: the panel does not fetch a reference genome.
4. The structural panel treated an empty response as an unreachable service.
   It now distinguishes connected-but-no-cards, waiting, and unavailable.
5. Asset search no longer silently swallows traversal failures. It reports
   directories visited per root, filenames encountered, and filesystem errors.
6. The release note and manifest narrative now describe matching stored
   metadata as agreement, not proof of unique artifact identity. Manifest
   hashes, expected versions, and verification rules were not changed.

## Browser observations

The real browser was exercised against the real service, not the fixture:

- Structure: `AGGGTTAGGGTTAGGGTTAGGG` displayed g4_tm 72.5 degrees C under the
  UI buffer K140 Na10 Mg1 pH7.4 37C. Applicability warnings remained visible.
- Mutation scan: 66 substitutions, numerical deltas and stated exclusions.
- Condition response: a potassium sweep displayed the g4_tm curve; g4_fold
  had too few supported points, and refused heads were identified.
- Batch: two input records were returned with applicability and refusal states.
- Chr22 Variant join: 183 variants, 183 AVI scores, zero structural deltas,
  and 183 unclassified. The full table contained 183 data rows plus its header.
  Model 0.4.6 and the verified hash appeared in the interface. The connected
  empty-structure message correctly described applicability under model rules.

These are software/inference checks, not validation of biological mechanisms or
mutation-effect accuracy. The chr22 window remains a negative-control software
case. Its lack of an applicable motif is not proof that DNA cannot fold.

## API evidence

`GET /ready` returned 200 and `ready: true`. Five concurrent `/predict` requests
with atlas neighbours enabled succeeded and returned identical g4_tm values
72.578 for the 22-nt control at the explicit physiological preset. That preset
has Na12, unlike the UI's Na10 above; the inputs are not identical.

The chr22 `/scan/variant` response independently confirmed the 183/183 join,
zero deltas and a matching artifact hash. The full response and summary are
saved in `../verification-20260916/chr22-real-service.json` and
`../verification-20260916/api-verification.json`.

## Verification

- Python suite including the real artifact and launcher regression checks:
  147 passed, 7 skipped.
- Frontend tests: 51 passed. Typecheck and production build passed; Vite still
  reports a large visualization chunk.
- Asset-search regression checks: 4 passed.
- Ruff for `quadcond` and `tests`, and `git diff --check`: passed.
- Skips remain for unavailable supporting data, full atlas, model JSON sidecar,
  and optional pysam. There is no claim that every data-dependent test ran.

## Search scope

The new search used the current 0.5.0 manifest, not the superseded hashes. It
examined 34 candidate files after encountering 314,891 filenames in 48,189
directories. Visited directories: C:\Users 44,411; D:\ 3,709; E:\ 69.
There were 67 filesystem errors. Directory exclusions, extension/size filters,
and exclusion of archive interiors still apply. Drive usage or elapsed time
alone does not establish traversal coverage. No full atlas matched.
See `../asset-search-20260916.txt` for the bounded scan report.

## Restart

Run from the repository directory, after stopping any existing server on 8765:

```powershell
.\.venv\Scripts\python.exe bridge\serve_table.py ..\alphagenome_bridge\results\chr22_example_20260915_105922_ce722e\quadcond_scores.csv
```

Alternatively, `scripts/start_workbench.ps1 -AtlasTable PATH_TO_CSV` uses the
same isolated Python. Without a table it starts the ordinary service.

For a new installation, create `.venv` and install `-e . -r requirements-model.txt`;
use `-e '.[dev]'` when running tests. The application/model dependency environment
is separate from the AlphaGenome exporter environment.

Remaining release/scientific work: publish downloadable assets and fill their
verified URLs, obtain the full training atlas for metric re-derivation, and
validate the integrated prioritization against independent measured variants.
Local readiness does not establish fresh-clone readiness or publication validity.
