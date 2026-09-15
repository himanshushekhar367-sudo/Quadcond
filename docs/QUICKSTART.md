# Quickstart

## 1. Install

```bash
python -m venv .venv && source .venv/bin/activate      # or: conda create -n quadcond python=3.11
pip install -e ".[app,dev]"
```

## 2. Build the atlas

```bash
python scripts/01_build_atlas.py
quadcond atlas summary
quadcond atlas sources          # what is in there, at which evidence tier, with DOIs
```

The bundled sources live in `data/external/downloads/`. If you cloned this without them,
re-fetch:

```bash
mkdir -p data/external/downloads/g4stab_db
curl -sSLo "data/external/downloads/G4 Dataset.xlsx" \
  "https://raw.githubusercontent.com/donn-liew/G4ShapePredictor/main/g4sp%20supplementary/G4%20Dataset.xlsx"
for i in $(seq -w 0 18); do
  curl -sSLo "data/external/downloads/g4stab_db/pred_seqs_0${i}.csv" \
    "https://raw.githubusercontent.com/donn-liew/g4stab-web-database/main/data/pred_seqs_0${i}.csv"
done
```

## 3. Train and report

```bash
python scripts/02_train.py         # ~20 min on 2 cores; trains the ablation too
python scripts/03_model_card.py    # regenerates docs/MODEL_CARD.md from the model
quadcond report                    # artifacts/quadcond_report.html
python scripts/04_case_studies.py  # hTelo, c-MYC, c-KIT, BCL2, VEGF across six buffers
pytest -q
```

Faster iteration while developing:

```bash
python scripts/01_build_atlas.py --skip-predicted      # experimental + derived only
python scripts/02_train.py --tasks g4_fold g4_topology --seeds 2 --folds 3 --no-ablation
```

## 4. Use it

```bash
quadcond predict AGGGTTAGGGTTAGGGTTAGGG --preset physiological
quadcond predict AGGGTTAGGGTTAGGGTTAGGG --k 0 --na 100          # same sequence, Na+ buffer
quadcond sweep   AGGGTTAGGGTTAGGGTTAGGG --vary k --start 0 --stop 150 --out ktitration.csv
quadcond competition GGGTTAGGGTTAGGGTTAGGGTTAGGG --preset acidic_tumour
quadcond scan --fasta promoter.fa --preset nuclear --out elements.csv
quadcond atlas neighbours TGAGGGTGGGTAGGGTGGGTAA -n 8
quadcond info                                                    # metrics and provenance
streamlit run app/streamlit_app.py
```

## 5. Add your data

```bash
quadcond atlas ingest g4stab-supp Table1.xlsx  --dry-run     # check column matching
quadcond atlas ingest g4stab-supp Table1.xlsx
quadcond atlas ingest imseeker    TableS1.xlsx
quadcond atlas ingest lab cd_melts_2026.xlsx --source deep_lab_cd_2026
python scripts/02_train.py && python scripts/03_model_card.py && quadcond report
```

See `data/templates/lab_measurements_template.csv` for the expected shape of lab data.

## Reading the output

Every prediction carries three things beyond the number:

- **claim** — what that head is entitled to assert.
- **applicability** — whether your buffer falls inside the region the head actually saw.
  `[OUT OF DOMAIN]` means the answer is extrapolation; a variable that took one value in
  training is named explicitly, because the model cannot resolve it at all.
- **evidence** — the nearest real measurements in joint sequence × condition space.

A head marked `[not measurement-grounded]` was trained on derived or model-predicted
labels. Treat its output as a prior, never as evidence.

---

## What ships in the tarball

- Full source, tests and docs.

**What is *not* in the code archive.** The two trained models
(`artifacts/*.joblib`, ~28 MB together) and the atlas databases are distributed
separately, because a 250 KB source archive and a 300 MB one have different
homes and different lifetimes. `quadcond assets fetch` downloads them into a
cache directory and verifies the SHA-256 recorded in
`quadcond/assets_manifest.json` before anything is unpickled; `quadcond assets
status` says what is present. Point `--model` / `--db` at your own copies if you
already have them.

An earlier version of this page said the models and `data/atlas_core.db` were
included in the archive. They were not, and following the instructions on a
fresh clone failed at the first predict.
- `data/atlas_core.db` — the **experimental + derived** atlas (2,010 records: the 1,005
  G4SP topology measurements and their matched shuffles). Small, and enough for
  `quadcond predict`, `quadcond atlas neighbours` and the evidence panel to work the
  moment you unpack.
- `data/external/downloads/G4 Dataset.xlsx` — the experimental source table itself.
- The full 220,976-record atlas is **not** shipped (150 MB). Rebuild it with the two
  commands in step 2 above; it takes about 30 seconds once the shards are downloaded.

Point any command at the core atlas with `--db`:

```bash
quadcond --db data/atlas_core.db predict AGGGTTAGGGTTAGGGTTAGGG --preset k100
```
