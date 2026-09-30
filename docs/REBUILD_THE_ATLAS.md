# The full atlas: where it is, and how to get one

**It does not exist.** Not on this machine, not in any archive in the Quadcond
folder, and it was never published. Three searches agree:

- a filesystem search of `/mnt/c` and `/home` for `atlas*.db` and for files of the
  recorded size returns three copies of `atlas_core.db` (5,394,432 bytes each) and
  nothing else;
- `quadcond-repo.zip`, `release-v0.4.9.zip`, `sandbox.zip` and
  `cowork-v0.4.9-parallel.zip` contain `atlas_core.db` only;
- `quadcond/assets_manifest.json` records a fingerprint for `atlas.db` but its
  release URL is `null`, so nothing was ever uploaded.

`atlas_core.db` holds 8,938 records from 8 sources. The atlas the shipped v0.5.1
models were trained on held **398,375 records from 18 sources**. The core is 2.2%
of it, and is not a substitute.

## What survives

The frozen model sidecar recorded the original's full per-source composition, and
it is now broken out into:

- `docs/atlas/original_atlas_composition.csv` — every source, tier, label class,
  record count, unique sequences and condition count
- `docs/atlas/original_atlas_composition.json` — totals and the largest sources

By tier: 254,761 predicted, 88,643 derived, 54,971 experimental. Recorded dataset
fingerprint `4c5cca5a98e1259452824cc3f5aaad151a11671a4391b99f5abaa67e716d061b`.

That file is the specification. A rebuild is checked against it source by source;
without it a rebuild would be unverifiable.

## A rebuild is possible, and the derived rows are deterministic

`scripts/01_build_atlas.py` builds the atlas from raw sources. Every random step
in it is seeded — the shuffled negatives from `sha1(parent_source)`, the predicted
block with `seed=11`, the i-motif shuffles with `seed=13` — but **seeded is not the
same as reproducing what core holds, and it does not.** Comparing core against a
rebuild source by source (`docs/atlas/core_vs_rebuild.csv`):

| | result |
| --- | --- |
| The three experimental sources | **byte-identical content**, all three |
| The three shuffled sources | same row counts, **different content** |

So the adapters reproduce the measurements exactly, which is the part that
matters, and the derived negatives do not match the ones in core. Two
explanations, and they have different consequences:

- core predates the seeding fix (a comment in the builder records that the seed
  "changed on every run" before it was corrected), in which case the rebuild's
  shuffles are the deterministic ones and core's are legacy; or
- the shuffle is still nondeterministic for some other reason, which would mean no
  atlas containing derived rows can be reproduced.

The check that separates them is cheap: build twice to different paths and compare
the shuffled sources against each other. Stable across two runs means the first
explanation. That has not been done, and until it is, do not describe the derived
rows as reproducible — an earlier draft of this file did, on the strength of
reading the seeds rather than measuring the output.

    conda activate qcgw
    pip install openpyxl pypdf          # readers the adapters need, once
    python scripts/01_build_atlas.py --db data/atlas.db

Run it in `qcgw` or the checkout's `.venv`, not conda `base`: base happens to have
pandas, so the build starts and then fails partway, leaving a partial database
behind. Three readers are needed beyond the model environment's own pins —
`pandas` for the tables, `openpyxl` for the two `.xlsx` sources, and `pypdf` for
the 5DUVMA supplementary PDF. A missing one is not caught up front; it surfaces
as an ImportError in the middle of a source.

The builder deletes the target database before it starts, so a failed run leaves a
partial file and a rerun replaces it. A partial `data/atlas.db` is not usable and
not a smaller atlas — check the record count against
`docs/atlas/original_atlas_composition.csv` before using one.

Output lands at **`quadcond-repo/data/atlas.db`**.

## What is in place, and what is missing

Raw inputs go in `data/external/downloads/`. That directory did not exist; it now
holds the three that were findable in the Quadcond folder:

| Input | Status | What it contributes |
| --- | --- | --- |
| `Dataset (G4STAB) Supplementary Table 1.csv` | **in place** | 2,367 records, 2,274 usable Tm across 261 buffers — the source behind every condition-aware stability claim |
| `Supplementary Data Set.xlsx` (gkae092) | **in place** | 160 measured i-motif pH_T over 147 sequences |
| `gkag110_supplemental_file.pdf` (5DUVMA) | **in place** | 937 records: 379 Tm + 558 pH_T across 10 ionic strengths |
| `G4 Dataset.xlsx` | **missing** | 1,005 G4 topology records — the `g4_topology` head |
| `gse220882_peaks.jsonl` | **missing** | 8,482 CUT&Tag + 35,314 locus records. `GSE220882_RAW` is in the Quadcond folder but needs peak calling against hg38 first — see `scripts/06_gse220882_peaks.py` |
| `g4stab_db/pred_seqs_*.csv` | **missing** | 254,761 predicted records over 40,000 sequences: 64% of the atlas by row count, and the largest single block |

So a rebuild run today recovers the measurement-grounded sources — which is what
every quantitative claim in the manuscript rests on — and misses topology, the
genomic proxy heads and the predicted block. The builder skips a missing source
with a warning rather than failing, so check its output against
`original_atlas_composition.csv` before trusting the result.

## What a rebuild will and will not be

It will not carry the recorded fingerprint. Even with every input restored, the
partial and full builds differ, and SQLite file layout is not content. What can be
established is content equivalence per source, against the composition file.

Which means the honest statement for the paper does not change: the released
models were trained on an atlas that is no longer available, the reported per-head
metrics are a record of that training run, and everything published from here
should be built on an atlas that is itself published. See
`benchmarks/fold_provenance/README.md` for the matching problem on the fold side.


## Known regression: the 5DUVMA PDF loses table S9

A rebuild on 30 September recovered 6,698 of 398,375 rows, 6 of 18 sources. Most
of the gap is the three missing inputs, which is expected. One part is not.

`duvma_im_stability_landscape` came back **822 rows against 937, and 264 melting
temperatures against 379** — short by exactly 115, all of them Tm. The builder's
own per-table report says why:

    tables parsed: {'S1': 146, 'S4': 160, 'S9': 0, 'S10': 118, 'S12': 119, 'S13': 279}

S9 parsed zero rows, and 822 + 115 = 937. So S9 is the entire shortfall, and it
is where the missing melting temperatures live.

This matters more than its size. The `im_tm_condition` head is recorded in the
sidecar as trained on 379 condition-dependent Tm rows; a rebuild that yields 264
cannot retrain it faithfully. Fixing the S9 parse is therefore a prerequisite for
any recomputation involving that head, not a tidy-up.

The likely cause is the PDF reader. `duvma.load` calls `pypdf`, table extraction
from PDFs is notoriously version-sensitive, and the version used for the original
build was never pinned — `requirements-model.txt` covers the model stack only.
Next steps: check `duvma.py`'s S9 pattern against what `pypdf` 6.19.0 actually
returns for that page, and pin whichever version parses all six tables.

## Verifying any rebuild

    python docs/atlas/compare_rebuild.py --db data/atlas.db

Compares the rebuild against `docs/atlas/original_atlas_composition.csv` source by
source and writes `docs/atlas/original_vs_rebuild.csv`. Run it after every build:
the builder skips a missing source with a warning and exits successfully, so a
finished run is not evidence of a correct atlas.


## Assembling an atlas.db today

`docs/atlas/atlas_from_core.py` writes `data/atlas.db` from `atlas_core.db`, plus
any source a rebuild adds that core lacks. As of 30 September a rebuild adds
nothing: core is a superset of it, holding the same three experimental sources
byte-identically and 2,010 rows more (`g4sp_topology` and its shuffle, which
cannot be rebuilt because `G4 Dataset.xlsx` is not on disk).

    python docs/atlas/atlas_from_core.py --rebuild data/atlas_rebuild.db

The result is 8,938 records over 8 of the 18 sources — 2.2% — and it carries an
`assembly_provenance` row in `meta` saying so, so nothing downstream mistakes it
for the training atlas or quotes its fingerprint as the training fingerprint.

Note that the script writes to `data/atlas.db` by default, which will overwrite a
rebuild sitting at that path. Build rebuilds to their own filename.

What this means for the three heads whose data is now complete: `g4_tm`,
`im_pht`, `im_pht_condition` and `im_tm_condition` all have their full measured
training sets back, and `g4_topology` has core's copy. The genomic proxy heads and
the predicted block remain unavailable.


## 30 September: the GSE220882 peaks were on disk after all

`gse220882_peaks.jsonl.gz` (14 MB) exists in three places in the Quadcond folder —
`pipeline/`, `work/quadcond/data/external/downloads/`, and an uncompressed 89 MB
copy in `_to_delete/` — together with its `meta.json`. The builder's glob accepts
the `.gz`, and the adapter opens gzip directly. Copied into
`data/external/downloads/`, a rebuild recovers six sources **at exactly the
original counts**:

| Source | Original | Rebuild |
| --- | --- | --- |
| gse220882_hek293t_locus_state | 35,314 | 35,314 |
| shuffled::gse220882_hek293t_locus_state | 35,314 | 35,314 |
| gse220882_hek293t_im_cutandtag | 8,482 | 8,482 |
| shuffled::gse220882_hek293t_im_cutandtag | 7,725 | 7,725 |
| gse220882_hek293t_g4_cutandtag | 6,706 | 6,706 |
| shuffled::gse220882_hek293t_g4_cutandtag | 5,135 | 5,135 |

98,676 rows, and the peak calling reproduces exactly — same motif-carrying
fractions (G4 22.9% against a 2.1% genomic background, iM 37.3% against 3.9%) and
the same locus-state split. The rebuild reaches 105,604 rows; folding in core's
`g4sp_topology` pair brings it to **107,614 rows over 14 of 18 sources, 27.0%**.

## Everything still missing hangs off one file

`G4 Dataset.xlsx` is not on disk anywhere (the only G4-named workbooks are the
stG4 tables in `gkag821_supplemental_files`, a different dataset), but core holds
its 1,005 topology records, so that source is covered.

The other four are one dependency chain, and `build_complementary_im` is the link
that makes it so — it queries `atlas.query(kind="G4", tiers=("predicted",))`:

    g4stab_db/pred_seqs_*.csv   (absent everywhere; no CSV over 1 MB in the folder)
      -> g4stab_webdb_predicted              254,761
      -> complementary_strand_im_candidates   12,000
           -> im_architecture_positives       12,000
           -> im_architecture_negatives       12,000
                                             --------
                                             290,761   = 73% of the atlas

So the shards are not just the largest single block; they gate the two i-motif
architecture sets and the complementary-strand candidates as well. Nothing else
recovers them, and `--skip-predicted` exists precisely to build without them.

## Building the best atlas available today

    python scripts/01_build_atlas.py --db data/atlas_rebuild.db
    python docs/atlas/atlas_from_core.py --rebuild data/atlas_rebuild.db
    python docs/atlas/compare_rebuild.py --db data/atlas.db

Build to `atlas_rebuild.db`, not `atlas.db`: the assembly script writes the latter
and would overwrite a rebuild sitting there. Expect 107,614 rows, 14 sources, with
the four predicted-chain sources reported MISSING.


## Can the heads be retrained from this atlas?

Checked per head against the sources each one records in the sidecar
(`data/atlas.db`, 107,614 rows, 14 sources):

| Head | Trained on | Sources available |
| --- | --- | --- |
| g4_fold | 6,744 | yes |
| g4_topology | 1,005 | yes |
| g4_tm | 2,274 | yes |
| g4_tm_distilled | 90,000 | **no** — needs g4stab_webdb_predicted |
| im_fold | 320 | yes |
| im_pht | 160 | yes |
| im_pht_condition | 558 | yes |
| im_tm_condition | 379 | yes |
| im_fold_genomic | 16,207 | yes |
| g4_fold_genomic | 11,841 | yes |
| im_architecture | 16,000 | **no** — needs the two architecture sets |
| locus_peak_overlap_state | 70,628 | yes |

**Ten of twelve.** And the two that cannot be retrained are exactly the two the
manuscript already describes as "auxiliary heads built on predicted or derived
labels" and makes no performance claim about. Every head behind a claim in the
paper — the four regression heads in Table 1, the five compared in Table 4, the
genomic proxies — has its full training data.

### Sequence matters more than the decision

Retraining is worth doing: it replaces metrics that cannot be verified with
metrics that regenerate from published files. But two things have to be settled
first, or a retrain bakes in the same problem it is meant to fix.

1. **The shuffle determinism question.** `g4_fold` and `im_fold` train on shuffled
   negatives, and those still do not reproduce core's. Until two consecutive
   builds are shown to agree, their training data is not reproducible even with
   the atlas published.
2. **One head before twelve.** Retrain `g4_tm` alone on this atlas and compare
   against the recorded R² of 0.672. Its 2,274 experimental rows are
   byte-identical to the original, so the only thing that changed is the fold
   partition (the stable-sort fix). A result within a few hundredths says the
   full retrain is low-risk; a large move says the grouping change matters more
   than expected and should be understood first.

What a full retrain then costs: every QuadCond number in Tables 1 and 4 shifts,
the ablation and figures regenerate, and the competitor architectures have to be
retrained on the *new* folds or the comparison stops being matched. That last
part is the expensive one. The two auxiliary heads would ship untrained or be
dropped; dropping them is cleaner but changes "12 heads" throughout the package.


## 30 September, later: every input is now on disk

The G4STAB web-database shards (`pred_seqs_000..019.csv`, 20 files, 210 MB,
1,957,530 rows) and `G4 Dataset.xlsx` were both located in the Quadcond folder
and staged. **All 18 sources are now buildable, and core is no longer needed** —
the merge step exists only for the case where `G4 Dataset.xlsx` is missing.

A build with the shards but before the xlsx gave 398,199 rows over 16 sources,
with one discrepancy worth recording:

| Source | Original | Rebuild | |
| --- | --- | --- | --- |
| g4stab_webdb_predicted | 254,761 | **256,595** | +1,834 (+0.72%) |
| complementary_strand_im_candidates | 12,000 | 12,000 | exact |
| im_architecture_positives | 12,000 | 12,000 | exact |
| im_architecture_negatives | 12,000 | 12,000 | exact |

The sampling is seeded and takes 40,000 unique sequences over 5 conditions, so
the overage means the shard set itself differs slightly from the one used
originally — the web database was presumably updated in between. The arithmetic
closes: 398,199 + 2,010 (the g4sp pair) - 1,834 = 398,375.

This is `predicted` tier. It feeds `g4_tm_distilled` and nothing else, and no
measurement is affected. The honest description is equivalent data, not identical
data, and any rebuild-trained `g4_tm_distilled` should say so.

### Two bugs found in the process

**The assembly script corrupted its own output.** It copied the base database
with `shutil.copy2`, which copies the `.db` and its `-wal` as two separate,
non-atomic operations. The pair landed inconsistent, producing first a
UNIQUE-constraint error from a damaged index and then `database disk image is
malformed`. It now uses SQLite's online backup API, wraps the merge in one
transaction, and runs `PRAGMA integrity_check` before returning. The corrupt file
was renamed `data/atlas.db.CORRUPT` and should be deleted.

**Core and the rebuild disagree on schema again.** Their indexes are identical and
both report `schema_version = 4`, but core's `kind` CHECK is narrower. Two
databases claiming one schema version while differing in constraints is worth
fixing on its own: bump the version when a constraint changes.

### The build, now complete

    python scripts/01_build_atlas.py --db data/atlas.db
    python docs/atlas/compare_rebuild.py --db data/atlas.db

Straight to `data/atlas.db`. No merge, no copy, nothing to corrupt. Expect 18 of
18 sources and roughly 400,209 rows — the original 398,375 plus the 1,834-row
predicted overage.


## Determinism: settled, both halves

**Grouping** — fixed (`benchmarks/fold_provenance/`), and the regression test
`tests/test_grouping_order_invariance.py` passes.

**Shuffled negatives** — reproducible. Two levels of evidence:

- `matched_negatives` is deterministic on the same positives with the same seed,
  but consumes one RNG stream across them in order, so it *is* order-dependent:
  shuffling the positives gives a different negative set of the same size.
- It does not bite, because positives arrive in insertion order from a
  deterministic build. Two independently built atlases produced **byte-identical**
  shuffled sources for all five: g4stab (2,367), im_pht (160), duvma (937),
  gse220882 g4 CUT&Tag (5,135) and locus_state (35,314).

So core's shuffles differ because core predates the seeding fix, not because the
current pipeline is unreliable. `check_shuffle_determinism.py` re-runs both
levels. The consequence to remember: any change to ingest order silently changes
the negatives, so ingest order is now part of the reproducibility contract.

## The g4_tm probe

    python scripts/02_train.py --db data/atlas.db --tasks g4_tm \
      --out artifacts/probe_g4_tm.joblib --no-ablation

    g4_tm: n=2274 groups=407 r2=0.692 (13.7s)

Against the recorded 0.670 (sidecar) and 0.672 (recomputed OOF). The 2,274 rows
are byte-identical to the original; **the only thing that changed is the fold
partition** — 407 groups against the recorded 400.

**0.692 is not an improvement and must not be reported as one.** Same
architecture, same measurements, different partition. It is the size of the
partition-to-partition variation, and it is the honest reason Table 1's number
will move if the package is retrained: not a better model, a different split. The
0.15 R² gap between grouped and random splitting already in the manuscript is the
same phenomenon an order of magnitude larger.

`ece=nan` is expected: expected calibration error is not defined for a regression
head.

A ±0.02 move is the low-risk outcome the probe was run to check, so a full
retrain is safe to proceed with.

## Full retrain on the rebuilt atlas (2026-09-30)

`python scripts/02_train.py --db data/atlas.db --seeds 5 --folds 5` on the
400,209-row rebuild: 4,469 biophysical measurements, 50,502 genomic proxy,
256,595 predicted, 88,643 derived/shuffle. All twelve heads trained; no
failures.

### Every head against the artifact the manuscript was frozen from

Old values read from `manuscript/inputs/quadcond_model.json`
(`artifact_sha256: 204b4608…`, internal `version: 0.4.7`), which is the
sidecar the v0.5.1 model card and Tables 1 and 4 were generated from.

| head | metric | frozen | retrained | delta |
|---|---|---|---|---|
| `im_pht_condition` | r2 | 0.834554683480 | 0.834554683480 | **bit-identical** |
| `im_tm_condition` | r2 | 0.853153876542 | 0.853153876542 | **bit-identical** |
| `im_fold_genomic` | auroc | 0.778981 | 0.779557 | +0.000576 |
| `g4_fold` | auroc | 0.970190 | 0.970762 | +0.000572 |
| `locus_peak_overlap_state` | bal. acc. | 0.450356 | 0.450975 | +0.000619 |
| `g4_tm_distilled` | r2 | 0.963151 | 0.963666 | +0.000516 |
| `im_architecture` | auroc | 0.996608 | 0.996905 | +0.000297 |
| `im_fold` | auroc | 0.870430 | 0.871172 | +0.000742 |
| `g4_fold_genomic` | auroc | 0.764531 | 0.762804 | −0.001727 |
| `g4_topology` | bal. acc. | 0.723835 | 0.741277 | +0.017443 |
| `im_pht` | r2 | 0.591567 | 0.610239 | +0.018673 |
| `g4_tm` | r2 | 0.670059 | 0.692084 | +0.022025 |

**The two bit-identical heads are exactly the two condition-grouped heads.**
Their folds are partitioned on condition, not on sequence clustering, so the
grouping fix cannot touch them — and to twelve decimal places it does not.
This is independent confirmation of the fold-provenance diagnosis, arrived at
by a different route from the sidecar replay: the replay compared partitions,
this compares trained metrics, and the two agree on precisely which heads the
grouping controls.

Every head that moved is sequence-grouped, and the movement scales inversely
with n. The three movers are the three smallest sequence-grouped sets —
`im_pht` (160 rows), `g4_topology` (1,005), `g4_tm` (2,274) — at +0.019,
+0.017, +0.022. Every sequence-grouped head above 6,000 rows moved by less
than 0.002. That is partition variance behaving exactly as it should.

**None of this is an improvement, and none of it may be reported as one.**
Same architecture, same hyperparameters, same measurements — byte-identical
rows for every measured source. The only thing that changed is which
sequences share a fold. `g4_tm` 0.670 → 0.692 is the size of the
partition-to-partition variation on 2,274 rows in 407 groups, and it is the
honest reason Table 1's number moves on a retrain. The manuscript already
documents the same effect an order of magnitude larger, in the 0.15 R² gap
between grouped and random splitting.

Group counts shifted slightly everywhere, as the grouping fix intends:
`g4_fold` 1,619 → 1,616, `g4_topology` 617 → 616, `g4_tm` 400 → 407,
`g4_tm_distilled` 3,552 → 3,431, `im_fold` 175 → 174, `im_pht` 85 → 85,
`im_fold_genomic` 3,907 → 3,913, `g4_fold_genomic` 2,929 → 2,930,
`im_architecture` 3,849 → 3,822, `locus_peak_overlap_state` 3,940 → 3,934.

### The ablation reproduces

Conditions carry the thermodynamics and nothing else:

| head | delta (frozen) | delta (retrained) |
|---|---|---|
| `im_tm_condition` | +0.7898 | +0.7898 |
| `g4_tm` | +0.4158 | +0.3888 |
| `g4_tm_distilled` | +0.3344 | +0.3337 |
| `im_pht_condition` | +0.2092 | +0.2092 |
| everything else | \|delta\| < 0.007 | \|delta\| < 0.007 |

The four large deltas are the four heads with real buffer variation; the eight
near-zero deltas are heads whose rows do not vary the buffer enough for
conditioning to be learnable. That is a statement about the data, not the
method, and the model card says so.

### `locus_peak_overlap_state` balanced accuracy 0.451 is above chance

Four classes — `both`, `g4_only`, `im_only`, `neither` — so chance is 0.250.
0.451 is about 1.8x chance: a weak but real signal, not a broken head. Its
condition delta is +0.0001 because the attached condition for genomic-proxy
rows is nominal with every field flagged imputed. It also costs 391 s of the
~700 s total train time, which is worth weighing against what it contributes.

### `im_architecture` AUROC 0.997 is circular, and already labelled as such

Positives are sequences that pass a canonical four-C-tract motif call;
negatives are dinucleotide shuffles built with `reject_motif=True`, i.e.
filtered by that same call. The head therefore distinguishes sequences that
match a regex from sequences constructed not to match it, and 0.997 means the
features can recompute the regex. `scripts/01_build_atlas.py` flags every row
`architecture_label_only` / `not_folding_evidence`, and the model card's
per-head claim already opens "ARCHITECTURE ONLY, AND NEARLY DETERMINISTIC BY
CONSTRUCTION … it is NOT evidence of biological predictive power." No change
needed; recorded here so the audit is not repeated.

## Two artifact-provenance defects the retrain exposed

### 1. The version stamp went backwards

The regenerated model card header reads "package 0.5.1 ships the model
artifact stamped **0.4.1**", where the previous card read **0.4.6** and the
frozen manuscript sidecar reads **0.4.7**. `scripts/03_model_card.py` takes
this from `getattr(model, "version")`, which is the dataclass default
`version: str = "0.4.1"` in `quadcond/models/base.py:127` — unchanged since
the initial commit. So every artifact built from this repo stamps 0.4.1, and
the 0.4.6 / 0.4.7 artifacts were built before the import. The field is
hand-maintained, nothing validates it, and it currently ships a version
number lower than the artifact it replaces. Decide the intended model version
and set it before generating a card that ships.

### 2. The manuscript and the external-validation panels were frozen against
different binaries

- `manuscript/inputs/quadcond_model.json` → `artifact_sha256: 204b4608…`
- `benchmarks/external_validation/panels/plasma_cfdna_2026.tsv` and
  `panels/fam230_kcl_titration.tsv` → `model_artifact_sha256_at_freeze:
  493640193…`

Two different serializations. Both report `g4_tm` R² 0.670, and
`sandbox/quadcond_v0.4.6_results.html` independently reports 0.670, so the
two binaries are the same trained model written out twice under different
version stamps — the science agrees, the provenance record contradicts
itself. This predates the retrain and needs one of the two documents
corrected.

`493640193…` survives as
`sandbox/quadcond_model.reassembled-v0.4.6.joblib`, verified by hash, and is
now also in the repo as
`artifacts/quadcond_model.panel-freeze-49364019.joblib`, so the panels remain
scoreable.

`204b4608…` is **gone**. It lived only at `artifacts/quadcond_model.joblib`,
`*.joblib` is gitignored, and the retrain overwrote it. What that costs is
small and worth stating exactly: its cross-validated metrics survive in
`manuscript/inputs/quadcond_model.json` and in the git history of
`docs/MODEL_CARD.md`; its numbers were never reproducible anyway, because the
historical fold assignments were never saved — that is the hole this whole
exercise was opened to close; and v0.5.1 was never published, so no third
party holds it. What is lost is the ability to re-score anything with that
exact binary. The panels do not reference it.

### The structural fix

A single unversioned `artifacts/quadcond_model.joblib` that every script
loads by default, is gitignored, and is silently overwritten by a training
run is how both defects happened. Name shipped artifacts by version or
content hash, load them explicitly, and keep the bare path as a symlink at
most. Until then, this retrain is preserved as
`artifacts/quadcond_model.rebuilt-atlas-20260930.joblib` (+ `.json`) and
`artifacts/quadcond_model_seqonly.rebuilt-atlas-20260930.joblib`.

### What still has to happen before any of these numbers ship

The head-to-head benchmark in Tables 1 and 4 retrained every competitor
architecture on the *old* folds. Those folds no longer exist. Until
`benchmarks/published_tools` is re-run with the competitors on the new folds,
the retrained QuadCond numbers above cannot be placed in the same table as
the published competitor numbers — a matched comparison is the entire basis
of that table's fairness. Either re-run the benchmark, or keep reporting the
frozen v0.5.1 numbers and cite this section for the partition sensitivity.
Do not mix the two.
