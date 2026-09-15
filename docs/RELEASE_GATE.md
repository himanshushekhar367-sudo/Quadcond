# What has been run, and what has not

A release note that quotes a test count is a claim. This file is where the
claim is made checkable, because the previous one was not: v0.4.7 shipped a
JUnit file dated four days before the tree it described (74 cases against 117
collected, with the new mutation, condition-scan and provenance checks absent
from it), a changelog line reading "96 passing and 45 browser checks", and a
packaging script that excluded the entire `artifacts/` directory — so the
archive contained screenshots and a preview JSON and no machine-readable gate
result at all. One of the screenshots showed model 0.4.6.

None of that was a false test result. It was evidence that had come unstuck
from the thing it was evidence *for*, which is the failure mode a gate exists
to prevent, so it is treated here as a defect of the same kind as the ones in
the code.

## The rule

Three things must be true of any number quoted as a result in a changelog, a
paper draft or a README:

1. It came from a run against **this** tree, at or after the commit being
   described.
2. The artifact it came from is **in the release archive**, not on the machine
   that produced it.
3. The artifact **names what it was a run of** — source version, model version,
   and for the browser gate the service identity it connected to.

## What a full gate is

```
# QuadCond
python -m pytest tests -q --junitxml=artifacts/pytest-results.xml

# AENNA-3D
npm run verify        # typecheck → app tests → scaffold tests →
                      # render-check (dev) → build → render-check (production)
                      # → browser-smoke
```

`npm run verify` writes `artifacts/render-check/result.json` for the dev leg and
`artifacts/render-check/result-prod.json` for the production leg. Both are now
inside the source archive: `scripts/package-source.mjs` keeps them by name
despite the blanket `artifacts/` exclusion, and the archive is where a reviewer
looks.

The browser gate connects to a running QuadCond service. **It is not a test of
the application unless that service is serving the intended release model.**
Each verdict file records the service's reported `quadcond_version`,
`model_version` and synthetic-provenance banner for exactly this reason: a green
gate against a synthetic or stale model is a green gate about something else.

## Status for v0.4.8

| | run against this tree | evidence in the archive |
|---|---|---|
| QuadCond pytest | yes — 117 collected, 82 passed, 35 skipped, 0 failed | `artifacts/pytest-results.xml` |
| AENNA typecheck | yes — clean | — |
| AENNA app tests | yes — 77 passed | — |
| AENNA scaffold tests | yes — 200 passed | — |
| Production build | yes — succeeds | — |
| Browser gate (dev + production) | **no** | **absent** |
| Real-model inference and metric reproduction | **no** | **absent** |

The 35 skips are not silent. They are model-, atlas- or supplementary-file
gated, and each one names its reason in the pytest output; the trained
artifacts are not bundled and every `url` in `quadcond/assets_manifest.json` is
still null, so a source-only checkout cannot run them. That is a packaging
limitation, not a finding that the assets are unrecoverable.

The comparison workflows are no longer among the skips. `tests/
test_comparison_contracts.py` exercises strand routing, mutant identity,
refusal handling, batch reconciliation, the withdrawn competition combination
and the shared run record against sentinel estimators driving the real `Head`,
`Predictor` and `scans` code — no trained artifact required. Those are
control-flow properties, and control flow is what was wrong; pinning them
behind an asset nobody has is how they stayed wrong for four releases.

**The browser gate has not been run for v0.4.8 and no v0.4.8 render-check
verdict exists.** Nothing in this release quotes a browser-check count. Before
tagging, run `npm run verify` against a service serving the intended release
model, confirm both verdict files record that model's version, and only then
write the counts down.
