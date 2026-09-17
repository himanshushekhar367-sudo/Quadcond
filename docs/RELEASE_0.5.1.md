# Release 0.5.1 — what the benchmark changed

0.5.1 is the release that acts on `benchmarks/published_tools/`: a head-to-head
against twelve published tools, run after 0.5.0 shipped. Three findings needed
code, not prose.

## 1. The folding heads do not transfer to genomic windows

`g4_fold` scores AUROC 0.66 against random genomic windows and **0.29** against
motif-matching windows G4-seq did not observe, where G4Hunter reaches 0.97 and
0.69. Its 0.97 on the atlas task is a statement about dinucleotide shuffles, not
about genomes.

- The binary folding heads (`g4_fold`, `g4_fold_genomic`, `im_fold`,
  `im_fold_genomic`) now **refuse** any window longer than their training range,
  with a reason that names the genomic path. A flagged number that is worse than
  chance is not a number to hand out.
- New model, new asset, new command: `quadcond genome-scan`, `POST /scan/genome`,
  `quadcond_g4seq_scanner.joblib`, trained on human G4-seq (K+, GSE110582) with
  chr2/8/17 held out and mouse kept as a cross-species test. On the hard
  negatives it scores 0.928 held-out and 0.775 cross-species, against 0.672 and
  0.554 for G4Hunter. `docs/GENOME_SCAN.md` has the full table, the claim, and
  the failure mode (isolated motifs in G-poor context, handled by an explicit
  `canonical_motif_rescue` call basis rather than silently).

## 2. Motif-destroying variants were being hidden

In the AlphaGenome 2×2, a substitution that destroys the motif has no delta, and
the quadrant needed two numeric coordinates — so the most disruptive class of
variant landed in `unclassified` and sank to the bottom of a table ranked on
|delta|. Such a variant is now placed structurally high by category:
`structural.basis = "motif_lost"`, structural rank 1.0, quadrant `both` or
`structural_only`, and still **no invented delta**. The categorical shortlist is
unchanged, and the workbench shows "motif lost" where the delta would be.

## 3. The assets were unreachable

Every `url` in `assets_manifest.json` was `null`, so a colleague who cloned the
public repository could not run the tool at all. The model, the sequence-only
ablation, the core atlas and the new scanner are now published as GitHub release
assets and recorded with their URLs; `quadcond assets fetch` works from a clean
checkout. The 354 MB training atlas is still not published (it is needed only to
retrain, and it is not present on the release machine) and its URL stays `null`
with that stated.

## Also in this release

- `genomewide/`: a resumable, multi-core pipeline that enumerates every
  canonical G4 and i-motif on chr1–22 and chrX, scores every SNV inside them,
  draws a matched out-of-motif control window per motif, reads AlphaGenome AVI
  scores from the published Tabix bundle, and runs the comparison — including a
  within-motif paired test that holds the locus fixed.
- `alphagenome_bridge/positive_controls.py`: the join on six published promoter
  G4s, with the coordinates derived from the UCSC hg38 API rather than typed in.
  The chr22 window remains the negative control it always was.
- The interpretation note on a variant scan no longer says mutation-effect
  prediction is unvalidated: it is now benchmarked on 1,020 measured
  single-substitution Tm pairs (Spearman 0.54, direction correct 81 % at
  |ΔTm| ≥ 2 °C). The *combined* ranking remains unvalidated, and says so.

## Publishing the assets

The manifest URLs point at the `v0.5.1` release, which has to exist for
`quadcond assets fetch` to work. With the GitHub CLI:

```powershell
winget install --id GitHub.cli        # if `gh` is not installed
gh auth login
gh release create v0.5.1 `
  artifacts\quadcond_model.joblib artifacts\quadcond_model_seqonly.joblib `
  artifacts\quadcond_g4seq_scanner.joblib data\atlas_core.db `
  --title "QuadCond 0.5.1" --notes-file docs\RELEASE_0.5.1.md
```

Without it, do the same thing in the browser: Releases -> Draft a new release ->
tag `v0.5.1` on `main` -> attach those four files -> publish. Either way the
file names must match `assets_manifest.json`, because the manifest checks each
download against the SHA-256 recorded there.

Then verify from a clean directory:

```powershell
quadcond assets fetch
quadcond assets status
```

## Versions

| | |
| --- | --- |
| `quadcond` package | 0.5.1 |
| AENNA-3D workbench | 0.5.1 |
| asset manifest | 0.5.1 |
| model artifact | **0.4.6** (unchanged: the estimators did not move) |
| genome scanner artifact | g4seq-scanner-1 (new) |

The scientific frozen core of 0.4.x is untouched. Nothing in this release
retrains or modifies the twelve heads.
