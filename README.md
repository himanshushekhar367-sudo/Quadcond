# QuadCond / AENNA-3D

Condition-aware prediction of G-quadruplex and i-motif propensity, with a
browser workbench and an optional join to AlphaGenome regulatory scores —
served by a single process.

The repository holds three parts that ship together:

| Part | What it is |
| --- | --- |
| `quadcond/` | Python package: 12 calibrated prediction heads, a sequence/conditions feature stack, a motif atlas, a CLI and an HTTP service |
| `web/` | AENNA-3D: a React + Three.js workbench that talks to that service |
| `bridge/` | Helpers for pulling AlphaGenome variant scores and joining them to a QuadCond mutation scan |

`quadcond serve` serves both the API and the built workbench on one port, so a
deployment is one process and one URL.

---

## Install

```bash
git clone https://github.com/himanshushekhar367-sudo/Quadcond.git
cd Quadcond
python -m pip install -e .
```

Python 3.11 or newer. The pinned scientific stack (scikit-learn 1.8.0,
xgboost 3.2.0) is what the released model files were built with; loading them
under other versions is allowed but reported in every run's provenance block.

### Assets

Model and atlas files are **not** in git — they are binary artefacts of a
training run, verified by SHA-256 against `quadcond/assets_manifest.json`.

```bash
quadcond assets status     # where each asset is expected and whether it verifies
quadcond assets fetch      # download and verify the recorded assets
quadcond assets path       # resolved model and atlas paths
```

`docs/GET_THE_DATA.md` lists every external file the project can use, its
licence, and what the software does when it is absent. Nothing is downloaded
implicitly and no asset is installed unverified.

---

## Run the web server

```bash
npm --prefix web install
npm --prefix web run build      # writes web/dist
quadcond serve --port 8765      # or, without installing: python -m quadcond.service --port 8765
```

Then open `http://localhost:8765/`. The startup banner says whether it found a
built workbench, so you know before opening the browser whether `/` will be the
application or the health document.

`--host 0.0.0.0` accepts connections from other machines. There are no
accounts and no sessions, so do that only on a network you trust.

The server resolves the built workbench in this order: `QUADCOND_WEB_ROOT`,
then `quadcond/_web`, then `web/dist` beside the repository root. When no build
is present the service is API-only and `/` returns the health document, so an
existing API deployment is unchanged by this feature.

For frontend development, run the API and the dev server separately:

```bash
quadcond serve --port 8765          # terminal 1
npm --prefix web run dev            # terminal 2 — http://127.0.0.1:8080, proxies the API
```

### Endpoints

| Route | Method | Purpose |
| --- | --- | --- |
| `/health` | GET | liveness — 200 whenever the process can reply |
| `/ready` | GET | readiness — 503 until every required asset is present and verifies |
| `/info` | GET | version, assets, heads, and what each head may claim |
| `/schema/prediction` | GET | response schema |
| `/predict` | POST | score one or more sequences |
| `/scan/mutations` | POST | every single-base substitution across a window |
| `/scan/conditions` | POST | response across a condition grid |
| `/scan/variant` | POST | mutation scan joined to AlphaGenome scores |
| `/batch` | POST | many sequences in one call |

`/health` and `/ready` are deliberately different questions: a process whose
required model file is missing is alive but not ready, and a load balancer
should be able to tell.

---

## Command line

```bash
quadcond info                                   # version, assets, heads, claim basis
quadcond predict GGGTTAGGGTTAGGGTTAGGG --k 100
quadcond scan --fasta window.fa --heads g4_tm,im_pht
quadcond competition <sequence>                 # both strands, scored separately
quadcond sweep <sequence> --vary k --start 0 --stop 150
quadcond atlas summary
quadcond report --out report.html
```

`competition` reports each strand on its own. It does not emit a combined
"which structure wins" score: the two motifs sit on opposite strands of the
same duplex, and no combination rule in this package is calibrated against
data that would support one.

---

## AlphaGenome

`quadcond variant` joins a QuadCond mutation scan to AlphaGenome variant
scores for the same edits, matching exactly on `(chromosome, position, ref,
alt)` — never by proximity.

```bash
# No key needed: the published AVI Tabix bundle
quadcond variant --fasta hg38.fa --chromosome chr22 --start 10510000 --end 10510200 \
    --avi-tabix /data/alphagenome/avi.tsv.gz

# Cheap pre-check: motif screening only — no model, no query, no key
quadcond variant --fasta hg38.fa --chromosome chr22 --start 10510000 --end 10510200 --screen

# A CSV/TSV export you already have
quadcond variant ... --atlas-table scores.tsv

# Live service (needs ALPHAGENOME_API_KEY)
quadcond variant ... --live
```

Four sources are supported behind one contract — Tabix bundle, table export,
live API, and a null source — and the resolver prefers the Tabix bundle
because it is the permissively licensed one.

`--screen` exists because most windows have no canonical motif at all. It
checks the reference **and every single-base alternate on both strands** and
returns one of `motif_present`, `motif_gain_possible`, `borderline` or
`no_motif`, so a query is only spent on a window where a structural delta is
possible.

`docs/ALPHAGENOME.md` covers the data files, their licences, and the exact
semantics of the join.

### Worked example

`examples/chr22-negative-control/` is a complete run: 183 variants, all 183
paired to a regulatory score, and **zero** structural deltas — the window
carries no canonical motif on either strand (G4Hunter 0.262). It is kept
precisely because it is negative. It shows the join, the provenance block and
the refusal path working on a window where the honest answer is "nothing to
report", which is the common case.

---

## What this software refuses to do

Refusal is a first-class result here, not an error path. A head that is asked
about a sequence outside the region its training data covers returns a refusal
with a reason, and the frontend renders that refusal rather than a number.
Downstream code reads it through a discriminated union, so treating a refused
head as a score is a compile error, not a silent zero.

Every run carries a provenance block: the model file's SHA-256, a three-valued
`matches_manifest` (`true` / `false` / `unknown` — never a bare `false` when
the manifest simply has no entry), the build environment, and whether the
current environment matches it.

`docs/CLAIMS.md` states, head by head, the evidence tier, what the target
actually measures, and whether the head is biophysically grounded. Only some
heads have a learned response to conditions; the ones that do not are labelled
and are not described as condition-responsive anywhere in the interface.

---

## Tests

```bash
python -m pytest tests -q          # backend
npm --prefix web run typecheck     # frontend types
npm --prefix web test              # frontend unit tests
```

The suite runs without the model binaries: sentinel estimators stand in for
them, so contract and refusal behaviour is testable on a clean checkout.

---

## Layout

```
quadcond/       prediction package, CLI, HTTP service, atlas, claim registry
web/            AENNA-3D workbench (React, TypeScript, Three.js, Vite)
bridge/         AlphaGenome export, table server, join runner
examples/       complete worked runs with provenance
docs/           claims, model card, data sources, release gate, quickstart
tests/          backend test suite
```

## Licence

MIT — see `LICENSE`. Third-party data files carry their own terms; the
AlphaGenome AVI bundle is permissively licensed, while the splicing and SHAP
bundles are not available for commercial use. `docs/DATA_SOURCES.md` records
the terms for each.

## Citation

See `CITATION.cff`.
