# QuadCond model card

*Generated 2026-09-09 by `scripts/03_model_card.py` from `artifacts/quadcond_model.joblib`. Do not edit by hand.*

## Summary

**7 biophysically anchored · 3 genomic-proxy · 2 derived/predicted auxiliary**

`claim basis` is the field to read. It is not the same question as "was anything measured": a CUT&Tag peak is an experimental observation whose target is antibody occupancy at a locus, so it carries `has_experimental_observation = true` and `biophysically_grounded = false`. Only the latter licenses a folding or stability claim.

| head | kind | task | claim basis | experimental? | biophysical? | rows | key metric |
|---|---|---|---|---|---|---|---|
| `g4_fold` | G4 | binary | biophysical | yes | **yes** | 6,744 | auroc = 0.9702 |
| `g4_topology` | G4 | multiclass | biophysical | yes | **yes** | 1,005 | balanced_accuracy = 0.7238 |
| `g4_tm` | G4 | regression | biophysical | yes | **yes** | 2,274 | r2 = 0.6701 |
| `g4_tm_distilled` | G4 | regression | predicted | no | no | 90,000 | r2 = 0.9632 |
| `im_fold` | iM | binary | biophysical | yes | **yes** | 320 | auroc = 0.8704 |
| `im_pht` | iM | regression | biophysical | yes | **yes** | 160 | r2 = 0.5916 |
| `im_pht_condition` | iM | regression | biophysical | yes | **yes** | 558 | r2 = 0.8346 |
| `im_tm_condition` | iM | regression | biophysical | yes | **yes** | 379 | r2 = 0.8532 |
| `im_fold_genomic` | iM | binary | genomic_proxy | yes | no | 16,207 | auroc = 0.7790 |
| `g4_fold_genomic` | G4 | binary | genomic_proxy | yes | no | 11,841 | auroc = 0.7645 |
| `im_architecture` | iM | binary | derived | no | no | 16,000 | auroc = 0.9966 |
| `locus_peak_overlap_state` | locus | multiclass | genomic_proxy | yes | no | 70,628 | balanced_accuracy = 0.4504 |

> **Genomic proxy, not biophysics:** `im_fold_genomic`, `g4_fold_genomic`, `locus_peak_overlap_state`. The positive label is antibody occupancy at a genomic locus (iMab / BG4 CUT&Tag), which is an experimental observation of a different quantity from folding in a defined buffer -- and one whose interpretation is contested, since iMab specificity and the accessible-chromatin contribution to targeted CUT&Tag signal are both open questions. No buffer was measured for these rows; the attached condition is nominal with every field flagged imputed. Read these heads as resemblance to peak-enriched motifs, never as P(folds).

> **Derived / prior-model auxiliary:** `im_architecture`, `g4_tm_distilled`. No measurement of any kind stands behind these labels. They exist to exercise the pipeline and to supply a condition-response prior; their numbers are not evidence about DNA.

## Evaluation protocol

- **Grouped cross-validation.** Sequences are clustered on 4-mer cosine similarity and folds are built over clusters, so near-duplicates cannot straddle a split. Group statistics are reported per head below.
- **Out-of-fold metrics only.** Nothing is reported on data a model saw.
- **Cross-fitted calibration.** The calibrator is itself cross-fitted over groups before ECE/Brier are computed, so calibration is not scored on the calibrator's own training data. The deployed calibrator is separately fitted on all out-of-fold predictions.
- **OOF-residual interval.** Regression heads report a split-conformal half-width from out-of-fold residuals, with the realised coverage stated.

## Condition ablation

Same data, same protocol, condition and interaction features removed.

| head | metric | sequence + conditions | sequence only | delta |
|---|---|---|---|---|
| `g4_fold` | auroc | 0.9702 | 0.9686 | +0.0016 |
| `g4_topology` | balanced_accuracy | 0.7238 | 0.7205 | +0.0033 |
| `g4_tm` | r2 | 0.6701 | 0.2542 | +0.4158 |
| `g4_tm_distilled` | r2 | 0.9632 | 0.6288 | +0.3344 |
| `im_fold` | auroc | 0.8704 | 0.8714 | -0.0010 |
| `im_pht` | r2 | 0.5916 | 0.5902 | +0.0013 |
| `im_pht_condition` | r2 | 0.8346 | 0.6254 | +0.2092 |
| `im_tm_condition` | r2 | 0.8532 | 0.0634 | +0.7898 |
| `im_fold_genomic` | auroc | 0.7790 | 0.7788 | +0.0002 |
| `g4_fold_genomic` | auroc | 0.7645 | 0.7642 | +0.0004 |
| `im_architecture` | auroc | 0.9966 | 0.9969 | -0.0003 |
| `locus_peak_overlap_state` | balanced_accuracy | 0.4504 | 0.4515 | -0.0011 |

A delta near zero is a finding about the *data*, not the method: it means that head's training rows do not vary the buffer enough for conditioning to be learnable. Ingest condition-resolved measurements and re-run.

## Heads in detail

### `g4_fold`

**Claim.** P(sequence adopts a G-quadruplex) against a composition-matched dinucleotide-shuffled background. Pinned to its four biophysical sources by name. Filtering on label_class alone was not enough: the GSE220882 negatives are label_class=catalog like every other shuffle, so they passed the filter while their positives did not, and the head silently gained 5,135 unpaired negatives and an AUROC to match. g4_fold_genomic is the head that carries the peak-derived rows.

- kind: G4; task: binary; target column: `folded`
- training rows: 6,744; features: 160
- evidence tiers: derived, experimental
- sources: `g4sp_topology`, `g4stab_experimental_tm`, `shuffled::g4sp_topology`, `shuffled::g4stab_experimental_tm`
- sequence groups: 1,619 from 6,744 sequences (largest cluster 71, 511 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| auroc | 0.9702 |
| auprc | 0.9706 |
| brier_uncalibrated | 0.0691 |
| ece_uncalibrated | 0.0206 |
| calibration_method | platt |
| ece_isotonic | 0.0041 |
| ece_platt | 0.0081 |
| n_distinct_calibrated_values | 2959 |
| brier | 0.0678 |
| ece | 0.0081 |
| mce | 0.0268 |
| accuracy_at_0.5 | 0.9045 |
| positive_rate | 0.5000 |
| calibration_scope | ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 0 – 1015 | 57 | yes |
| na | 0 – 1020 | 39 | yes |
| li_nh4 | 0 – 300 | 19 | yes |
| mg | 0 – 10 | 4 | yes |
| ph | 4 – 8 | 21 | yes |
| temperature | 25 – 25 | 1 | yes |
| length | 5 – 99 | 59 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 0.2 – 300 | 19 | yes |

> This head saw a single value of: `temperature`, `crowder_pct`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `g4_topology`

**Claim.** Calibrated posterior over folding topology, given that the sequence forms a G4. Trained on K+ buffer data only.

- kind: G4; task: multiclass; target column: `topology`
- training rows: 1,005; features: 160
- evidence tiers: experimental
- sources: `g4sp_topology`
- sequence groups: 617 from 1,005 sequences (largest cluster 126, 513 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| balanced_accuracy | 0.7238 |
| accuracy | 0.7522 |
| auroc_ovr_macro | 0.9049 |
| brier_uncalibrated | 0.3380 |
| ece_uncalibrated | 0.0722 |
| temperature | 1.3457 |
| brier | 0.3278 |
| ece | 0.0359 |
| calibration_scope | ECE is measured over the three topology classes on sequences already known to form a G4, in K+ buffer. It says nothing about whether a sequence folds, and nothing about topology in Na+. |

| class | n | recall |
|---|---|---|
| parallel | 487 | 0.864 |
| antiparallel | 216 | 0.694 |
| hybrid | 302 | 0.613 |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 100 – 100 | 1 | yes |
| na | 0 – 0 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 7 – 7 | 1 | yes |
| temperature | 25 – 25 | 1 | yes |
| length | 13 – 72 | 41 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `g4_tm`

**Claim.** Melting temperature under the supplied buffer. Requires experimental Tm records (G4STAB Supplementary Table 1) in the atlas.

- kind: G4; task: regression; target column: `tm`
- training rows: 2,274; features: 160
- evidence tiers: experimental
- sources: `g4stab_experimental_tm`
- sequence groups: 400 from 2,274 sequences (largest cluster 305, 185 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| r2 | 0.6701 |
| rmse | 7.7564 |
| mae | 5.5793 |
| spearman | 0.8163 |
| oof_residual_halfwidth | 12.3615 |
| oof_residual_alpha | 0.1000 |
| in_sample_coverage_of_residual_pool | 0.9002 |
| split_conformal_coverage | 0.8888 |
| split_conformal_coverage_sd | 0.0267 |
| split_conformal_halfwidth_mean | 12.2071 |
| split_conformal_n_splits | 20 |
| split_conformal_note | half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90% |
| interval_note | interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool |
| calibration_scope | Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 0 – 1015 | 56 | yes |
| na | 0 – 1020 | 39 | yes |
| li_nh4 | 0 – 300 | 19 | yes |
| mg | 0 – 10 | 4 | yes |
| ph | 4 – 8 | 21 | yes |
| temperature | 25 – 25 | 1 | **no** — see note below |
| length | 5 – 99 | 56 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 0.2 – 300 | 19 | yes |

> `temperature`: NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask.

> This head saw a single value of: `temperature`, `crowder_pct`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `g4_tm_distilled`

**Claim.** DISTILLED: reproduces the G4STAB ensemble's Tm response to cation composition. Useful as a condition-response prior and as a benchmark; it is NOT trained on measurements and inherits every G4STAB bias.

- kind: G4; task: regression; target column: `tm`
- training rows: 90,000; features: 160
- evidence tiers: predicted
- sources: `g4stab_webdb_predicted`
- sequence groups: 3,552 from 90,000 sequences (largest cluster 5355, 74 singletons)
- CV: 5 grouped folds x 3 seeds

| metric | value |
|---|---|
| r2 | 0.9632 |
| rmse | 1.9067 |
| mae | 1.4765 |
| spearman | 0.9785 |
| oof_residual_halfwidth | 3.0373 |
| oof_residual_alpha | 0.1000 |
| in_sample_coverage_of_residual_pool | 0.9000 |
| split_conformal_coverage | 0.8961 |
| split_conformal_coverage_sd | 0.0151 |
| split_conformal_halfwidth_mean | 3.0290 |
| split_conformal_n_splits | 20 |
| split_conformal_note | half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90% |
| interval_note | interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool |
| calibration_scope | Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 0 – 100 | 4 | yes |
| na | 0 – 100 | 4 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 7 – 7 | 1 | yes |
| temperature | 25 – 25 | 1 | **no** — see note below |
| length | 14 – 50 | 37 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> `temperature`: NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask.

> This head saw a single value of: `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_fold`

**Claim.** P(sequence forms an i-motif) against a composition-matched dinucleotide-shuffled background. Positives are sequences with a measured transitional pH, so this is grounded in folding measurements -- but the panel is small and heavily focused on designed C-tract variants.

- kind: iM; task: binary; target column: `folded`
- training rows: 320; features: 160
- evidence tiers: derived, experimental
- sources: `im_pht_biophysical`, `shuffled::im_pht_biophysical`
- sequence groups: 175 from 320 sequences (largest cluster 43, 137 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| auroc | 0.8704 |
| auprc | 0.8838 |
| brier_uncalibrated | 0.1532 |
| ece_uncalibrated | 0.0896 |
| calibration_method | platt |
| ece_isotonic | 0.0783 |
| ece_platt | 0.0552 |
| n_distinct_calibrated_values | 293 |
| brier | 0.1483 |
| ece | 0.0552 |
| mce | 0.1267 |
| accuracy_at_0.5 | 0.7875 |
| positive_rate | 0.5000 |
| calibration_scope | ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 100 – 100 | 1 | yes |
| na | 10 – 10 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 7 – 7 | 1 | yes |
| temperature | 25 – 25 | 1 | yes |
| length | 15 – 116 | 43 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_pht`

**Claim.** Transitional pH from sequence, at the iM-Seeker reference buffer (100 mM KCl + 10 mM Na cacodylate). Answers 'which of these sequences folds at higher pH?'. Grouped by sequence similarity, so the number is what to expect on a sequence you have not measured. It carries no salt response, because every row was measured in one buffer.

- kind: iM; task: regression; target column: `ph_t`
- training rows: 160; features: 160
- evidence tiers: experimental
- sources: `im_pht_biophysical`
- sequence groups: 85 from 160 sequences (largest cluster 19, 58 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| r2 | 0.5916 |
| rmse | 0.3494 |
| mae | 0.2549 |
| spearman | 0.7878 |
| oof_residual_halfwidth | 0.6317 |
| oof_residual_alpha | 0.1000 |
| in_sample_coverage_of_residual_pool | 0.9125 |
| split_conformal_coverage | 0.9002 |
| split_conformal_coverage_sd | 0.0599 |
| split_conformal_halfwidth_mean | 0.6139 |
| split_conformal_n_splits | 20 |
| split_conformal_note | half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90% |
| interval_note | interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool |
| calibration_scope | Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 100 – 100 | 1 | yes |
| na | 10 – 10 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 7 – 7 | 1 | yes |
| temperature | 25 – 25 | 1 | yes |
| length | 15 – 116 | 43 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_pht_condition`

**Claim.** Transitional pH as a function of buffer, for the two constructs 5DUVMA measured across ten ionic strengths. GROUPED BY BUFFER, not by sequence: folds hold out whole ionic strengths, so the number answers 'does this extrapolate to a salt concentration we did not measure?'. It says nothing about a new sequence -- there are only two here. REFUSES ANY SEQUENCE OUTSIDE ITS TWO-CONSTRUCT ALLOWLIST: no value is returned, because this head gives essentially the same answer whatever sequence it is given, so there is no prediction to report rather than an uncertain one.

- kind: iM; task: regression; target column: `ph_t`
- training rows: 558; features: 160
- evidence tiers: experimental
- sources: `duvma_im_stability_landscape`
- sequence groups: 10 from 558 sequences (largest cluster 68, 0 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| r2 | 0.8346 |
| rmse | 0.1891 |
| mae | 0.1398 |
| spearman | 0.8753 |
| oof_residual_halfwidth | 0.2942 |
| oof_residual_alpha | 0.1000 |
| in_sample_coverage_of_residual_pool | 0.9014 |
| split_conformal_coverage | 0.8805 |
| split_conformal_coverage_sd | 0.0716 |
| split_conformal_halfwidth_mean | 0.2958 |
| split_conformal_n_splits | 20 |
| split_conformal_note | half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90% |
| interval_note | interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool |
| calibration_scope | Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 10 – 1000 | 10 | yes |
| na | 0 – 0 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 7 – 7 | 1 | yes |
| temperature | 30 – 60 | 31 | yes |
| length | 21 – 45 | 2 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head answers for 2 specific construct(s) and **returns no value for any other sequence**.

> This head saw a single value of: `na`, `li_nh4`, `mg`, `ph`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_tm_condition`

**Claim.** i-motif melting temperature as a function of pH and ionic strength, for the two 5DUVMA constructs. Grouped by buffer, same caveat as im_pht_condition: a condition-response model, not a sequence model. REFUSES ANY SEQUENCE OUTSIDE ITS TWO-CONSTRUCT ALLOWLIST: no value is returned, because this head gives essentially the same answer whatever sequence it is given, so there is no prediction to report rather than an uncertain one.

- kind: iM; task: regression; target column: `tm`
- training rows: 379; features: 160
- evidence tiers: experimental
- sources: `duvma_im_stability_landscape`
- sequence groups: 10 from 379 sequences (largest cluster 40, 0 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| r2 | 0.8532 |
| rmse | 6.4586 |
| mae | 4.7768 |
| spearman | 0.9265 |
| oof_residual_halfwidth | 11.1802 |
| oof_residual_alpha | 0.1000 |
| in_sample_coverage_of_residual_pool | 0.9024 |
| split_conformal_coverage | 0.9043 |
| split_conformal_coverage_sd | 0.0429 |
| split_conformal_halfwidth_mean | 11.3924 |
| split_conformal_n_splits | 20 |
| split_conformal_note | half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90% |
| interval_note | interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool |
| calibration_scope | Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 10 – 1000 | 10 | yes |
| na | 0 – 0 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 4.35 – 7.05 | 19 | yes |
| temperature | 25 – 25 | 1 | **no** — see note below |
| length | 21 – 45 | 2 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> `temperature`: NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask.

> This head answers for 2 specific construct(s) and **returns no value for any other sequence**.

> This head saw a single value of: `na`, `li_nh4`, `mg`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_fold_genomic`

**Claim.** GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. Calibrated score for separating canonical i-motif motifs found in reproducible iMab CUT&Tag peaks (live HEK293T, GSE220882) from composition-matched motifs taken from dinucleotide shuffles of the same peak windows. The positive label is antibody occupancy at a locus, not a measured folding event, and no buffer was measured -- the attached condition is a nominal guess with every field flagged imputed. Read it as 'does this look like the motifs enriched under iMab peaks', and never as P(folds). Its value over im_fold is sequence diversity: thousands of independent loci instead of variants of a handful of designed constructs.

- kind: iM; task: binary; target column: `folded`
- training rows: 16,207; features: 160
- evidence tiers: derived, experimental
- sources: `gse220882_hek293t_im_cutandtag`, `shuffled::gse220882_hek293t_im_cutandtag`
- sequence groups: 3,907 from 16,207 sequences (largest cluster 50, 2,216 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| auroc | 0.7790 |
| auprc | 0.7903 |
| brier_uncalibrated | 0.1914 |
| ece_uncalibrated | 0.0230 |
| calibration_method | platt |
| ece_isotonic | 0.0056 |
| ece_platt | 0.0068 |
| n_distinct_calibrated_values | 7476 |
| brier | 0.1908 |
| ece | 0.0068 |
| mce | 0.0183 |
| accuracy_at_0.5 | 0.7036 |
| positive_rate | 0.5234 |
| calibration_scope | ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 140 – 140 | 1 | yes |
| na | 10 – 10 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0.5 – 0.5 | 1 | yes |
| ph | 7.2 – 7.2 | 1 | yes |
| temperature | 37 – 37 | 1 | yes |
| length | 15 – 48 | 34 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `g4_fold_genomic`

**Claim.** GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. The BG4 counterpart of im_fold_genomic, on the same peaks and the same matched-shuffle control. It exists mainly as a control: the two heads share a pipeline, so a gap between them says something about the antibodies and the underlying motif grammar rather than about the model.

- kind: G4; task: binary; target column: `folded`
- training rows: 11,841; features: 160
- evidence tiers: derived, experimental
- sources: `gse220882_hek293t_g4_cutandtag`, `shuffled::gse220882_hek293t_g4_cutandtag`
- sequence groups: 2,929 from 11,841 sequences (largest cluster 37, 1,192 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| auroc | 0.7645 |
| auprc | 0.7973 |
| brier_uncalibrated | 0.1916 |
| ece_uncalibrated | 0.0383 |
| calibration_method | isotonic |
| ece_isotonic | 0.0083 |
| ece_platt | 0.0214 |
| n_distinct_calibrated_values | 257 |
| brier | 0.1903 |
| ece | 0.0083 |
| mce | 0.0140 |
| accuracy_at_0.5 | 0.7062 |
| positive_rate | 0.5663 |
| calibration_scope | ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 140 – 140 | 1 | yes |
| na | 10 – 10 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0.5 – 0.5 | 1 | yes |
| ph | 7.2 – 7.2 | 1 | yes |
| temperature | 37 – 37 | 1 | yes |
| length | 15 – 33 | 19 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `im_architecture`

**Claim.** ARCHITECTURE ONLY, AND NEARLY DETERMINISTIC BY CONSTRUCTION. The positive label is 'carries a canonical four-C-tract motif' and the negatives are shuffles from which that motif was rejected, so a near-perfect AUROC here means the features can recompute the regex -- it is NOT evidence of biological predictive power. This head is a placeholder that exists so the iM pipeline is exercisable end to end; im_fold and im_fold_genomic are the heads that now do the work it was standing in for.

- kind: iM; task: binary; target column: `folded`
- training rows: 16,000; features: 160
- evidence tiers: derived
- sources: `im_architecture_negatives`, `im_architecture_positives`
- sequence groups: 3,849 from 16,000 sequences (largest cluster 74, 1,976 singletons)
- CV: 5 grouped folds x 3 seeds

| metric | value |
|---|---|
| auroc | 0.9966 |
| auprc | 0.9973 |
| brier_uncalibrated | 0.0131 |
| ece_uncalibrated | 0.0077 |
| calibration_method | platt |
| ece_isotonic | 0.0052 |
| ece_platt | 0.0067 |
| n_distinct_calibrated_values | 1505 |
| brier | 0.0130 |
| ece | 0.0067 |
| mce | 0.0270 |
| accuracy_at_0.5 | 0.9844 |
| positive_rate | 0.5009 |
| calibration_scope | ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability. |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 100 – 100 | 1 | yes |
| na | 0 – 0 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0 – 0 | 1 | yes |
| ph | 5.8 – 5.8 | 1 | yes |
| temperature | 25 – 25 | 1 | yes |
| length | 15 – 50 | 36 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

### `locus_peak_overlap_state`

**Claim.** GENOMIC PEAK-OVERLAP PROXY FOR A 201-NT WINDOW. NOT CO-OCCUPANCY AND NOT COMPETITION. Calibrated posterior over how a 201-nt genomic window was classified by two CUT&Tag experiments run as PARALLEL REACTIONS on separate aliquots of the same HEK293T population: reproducible peaks from both antibodies, from BG4 only, from iMab only, or from neither (a composition-matched shuffle of the same windows -- generated, not an observed unoccupied locus). Because the two antibodies never touched the same aliquot, a 'both' call is population-level co-localisation of peaks; nothing here observes two structures on one DNA molecule, and nothing here is a free energy. Balanced accuracy 0.450 against a 0.250 four-class floor, iM-only recall 0.156 -- the class most users are asking about is the one it recovers worst. Trained exclusively on 201-nt windows, so a query of oligonucleotide length is outside the applicability domain by construction. Read it as an exploratory genomic-resemblance score.

- kind: locus; task: multiclass; target column: `topology`
- training rows: 70,628; features: 182
- evidence tiers: derived, experimental
- sources: `gse220882_hek293t_locus_state`, `shuffled::gse220882_hek293t_locus_state`
- sequence groups: 3,940 from 70,628 sequences (largest cluster 624, 3,076 singletons)
- CV: 5 grouped folds x 5 seeds

| metric | value |
|---|---|
| balanced_accuracy | 0.4504 |
| accuracy | 0.6143 |
| auroc_ovr_macro | 0.8000 |
| brier_uncalibrated | 0.5008 |
| ece_uncalibrated | 0.0162 |
| temperature | 0.9369 |
| brier | 0.5004 |
| ece | 0.0101 |
| calibration_scope | ECE is measured over the three topology classes on sequences already known to form a G4, in K+ buffer. It says nothing about whether a sequence folds, and nothing about topology in Na+. |

| class | n | recall |
|---|---|---|
| both | 13156 | 0.392 |
| g4_only | 15088 | 0.353 |
| im_only | 7070 | 0.156 |
| neither | 35314 | 0.901 |

**Applicability domain** — the model saw only this much variation, so queries outside it are flagged at prediction time:

| variable | range | distinct values seen | enforced |
|---|---|---|---|
| k | 140 – 140 | 1 | yes |
| na | 10 – 10 | 1 | yes |
| li_nh4 | 0 – 0 | 1 | yes |
| mg | 0.5 – 0.5 | 1 | yes |
| ph | 7.2 – 7.2 | 1 | yes |
| temperature | 37 – 37 | 1 | yes |
| length | 201 – 201 | 1 | yes |
| crowder_pct | 0 – 0 | 1 | yes |
| strand_conc | 5 – 5 | 1 | yes |

> This head saw a single value of: `k`, `na`, `li_nh4`, `mg`, `ph`, `temperature`, `length`, `crowder_pct`, `strand_conc`. It cannot resolve those variables, and predictions that vary them are extrapolation, not inference.

## Intended use

- Prioritising candidate G4/i-motif elements for experimental follow-up, *under a stated buffer*, with a calibrated probability rather than a rank.
- Asking how a prediction changes across a cation or pH gradient (`quadcond sweep`), including whether the model is entitled to answer.
- Retrieving the nearest real measurements to any query (`quadcond atlas neighbours`), so a score can be checked against evidence.

## Out of scope

- **Structure prediction.** No 3D coordinates, no ligand docking.
- **RNA.** The bundled data is DNA; RNA sequences are accepted (U→T) but no RNA-specific head has been trained.
- **In-cell occupancy.** Nothing here models chromatin, transcription, or protein binding. A high in-vitro folding probability is not a claim about a structure existing in a nucleus.
- **Clinical or diagnostic use.**

## Known biases

- The topology data is predominantly human and predominantly 20–25 nt; long-looped, bulged and multimeric species are under-represented.
- The topology classes are imbalanced (parallel is the plurality), so balanced accuracy and per-class recall are reported alongside accuracy.
- Distillation heads inherit every bias of the model they distil, including its training-buffer coverage.

## Reproducing

```bash
python scripts/01_build_atlas.py
python scripts/02_train.py
python scripts/03_model_card.py
quadcond report
```

<details><summary>Raw model summary (JSON)</summary>

```json
{
  "version": "0.4.7",
  "dataset_fingerprint_sha256": "4c5cca5a98e1259452824cc3f5aaad151a11671a4391b99f5abaa67e716d061b",
  "artifact_sha256": "",
  "atlas_snapshot": {
    "path": "data/atlas.db",
    "n_records": 398375,
    "summary": [
      {
        "source": "g4stab_webdb_predicted",
        "kind": "G4",
        "evidence_tier": "predicted",
        "label_class": "biophysical",
        "n": 254761,
        "n_folded": 0,
        "n_topology": 0,
        "n_tm": 254761,
        "n_pht": 0,
        "n_unique_seq": 40000,
        "n_conditions": 5
      },
      {
        "source": "gse220882_hek293t_locus_state",
        "kind": "locus",
        "evidence_tier": "experimental",
        "label_class": "genomic_proxy",
        "n": 35314,
        "n_folded": 0,
        "n_topology": 35314,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 35312,
        "n_conditions": 1
      },
      {
        "source": "shuffled::gse220882_hek293t_locus_state",
        "kind": "locus",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 35314,
        "n_folded": 0,
        "n_topology": 35314,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 35314,
        "n_conditions": 1
      },
      {
        "source": "complementary_strand_im_candidates",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 12000,
        "n_folded": 0,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 11715,
        "n_conditions": 1
      },
      {
        "source": "im_architecture_negatives",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 12000,
        "n_folded": 12000,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 11991,
        "n_conditions": 1
      },
      {
        "source": "im_architecture_positives",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 12000,
        "n_folded": 12000,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 11715,
        "n_conditions": 1
      },
      {
        "source": "gse220882_hek293t_im_cutandtag",
        "kind": "iM",
        "evidence_tier": "experimental",
        "label_class": "genomic_proxy",
        "n": 8482,
        "n_folded": 8482,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 8123,
        "n_conditions": 1
      },
      {
        "source": "shuffled::gse220882_hek293t_im_cutandtag",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 7725,
        "n_folded": 7725,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 7724,
        "n_conditions": 1
      },
      {
        "source": "gse220882_hek293t_g4_cutandtag",
        "kind": "G4",
        "evidence_tier": "experimental",
        "label_class": "genomic_proxy",
        "n": 6706,
        "n_folded": 6706,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 6340,
        "n_conditions": 1
      },
      {
        "source": "shuffled::gse220882_hek293t_g4_cutandtag",
        "kind": "G4",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 5135,
        "n_folded": 5135,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 5121,
        "n_conditions": 1
      },
      {
        "source": "g4stab_experimental_tm",
        "kind": "G4",
        "evidence_tier": "experimental",
        "label_class": "biophysical",
        "n": 2367,
        "n_folded": 2367,
        "n_topology": 0,
        "n_tm": 2274,
        "n_pht": 0,
        "n_unique_seq": 1062,
        "n_conditions": 261
      },
      {
        "source": "shuffled::g4stab_experimental_tm",
        "kind": "G4",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 2367,
        "n_folded": 2367,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 2346,
        "n_conditions": 261
      },
      {
        "source": "g4sp_topology",
        "kind": "G4",
        "evidence_tier": "experimental",
        "label_class": "biophysical",
        "n": 1005,
        "n_folded": 1005,
        "n_topology": 1005,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 1005,
        "n_conditions": 1
      },
      {
        "source": "shuffled::g4sp_topology",
        "kind": "G4",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 1005,
        "n_folded": 1005,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 1005,
        "n_conditions": 1
      },
      {
        "source": "duvma_im_stability_landscape",
        "kind": "iM",
        "evidence_tier": "experimental",
        "label_class": "biophysical",
        "n": 937,
        "n_folded": 937,
        "n_topology": 0,
        "n_tm": 379,
        "n_pht": 558,
        "n_unique_seq": 2,
        "n_conditions": 192
      },
      {
        "source": "shuffled::duvma_im_stability_landscape",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 937,
        "n_folded": 937,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 912,
        "n_conditions": 192
      },
      {
        "source": "im_pht_biophysical",
        "kind": "iM",
        "evidence_tier": "experimental",
        "label_class": "biophysical",
        "n": 160,
        "n_folded": 160,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 160,
        "n_unique_seq": 147,
        "n_conditions": 1
      },
      {
        "source": "shuffled::im_pht_biophysical",
        "kind": "iM",
        "evidence_tier": "derived",
        "label_class": "catalog",
        "n": 160,
        "n_folded": 160,
        "n_topology": 0,
        "n_tm": 0,
        "n_pht": 0,
        "n_unique_seq": 160,
        "n_conditions": 1
      }
    ]
  },
  "heads": {
    "g4_fold": {
      "name": "g4_fold",
      "kind": "G4",
      "task": "binary",
      "target": "folded",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "auroc": 0.9701901888273958,
        "auprc": 0.9706190086646269,
        "brier_uncalibrated": 0.06906744285185676,
        "ece_uncalibrated": 0.020550326985585624,
        "calibration_method": "platt",
        "ece_isotonic": 0.004137041876253732,
        "ece_platt": 0.008075821404432828,
        "n_distinct_calibrated_values": 2959,
        "brier": 0.06784587816580487,
        "ece": 0.008075821404432828,
        "mce": 0.026834819493397832,
        "accuracy_at_0.5": 0.9045077105575327,
        "reliability": {
          "confidence": [
            0.0029596087951287062,
            0.009602966946567568,
            0.027991091533521688,
            0.09245420660925915,
            0.32757343392770616,
            0.672235412965208,
            0.9027712522120207,
            0.9750784853419635,
            0.9916270591565026,
            0.9973708541814641
          ],
          "accuracy": [
            0.0,
            0.002967359050445104,
            0.03115727002967359,
            0.10518518518518519,
            0.33827893175074186,
            0.6454005934718101,
            0.914074074074074,
            0.973293768545994,
            0.9896142433234422,
            1.0
          ],
          "count": [
            675,
            674,
            674,
            675,
            674,
            674,
            675,
            674,
            674,
            675
          ]
        },
        "positive_rate": 0.5,
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability."
      },
      "training_meta": {
        "n_rows": 6744,
        "n_features": 160,
        "tiers": [
          "derived",
          "experimental"
        ],
        "sources": [
          "g4sp_topology",
          "g4stab_experimental_tm",
          "shuffled::g4sp_topology",
          "shuffled::g4stab_experimental_tm"
        ],
        "grouping": {
          "n_items": 6744,
          "n_groups": 1619,
          "largest_group": 71,
          "singleton_groups": 511,
          "mean_group_size": 4.165534280420013
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 83.2,
        "claim": "P(sequence adopts a G-quadruplex) against a composition-matched dinucleotide-shuffled background. Pinned to its four biophysical sources by name. Filtering on label_class alone was not enough: the GSE220882 negatives are label_class=catalog like every other shuffle, so they passed the filter while their positives did not, and the head silently gained 5,135 unpaired negatives and an AUROC to match. g4_fold_genomic is the head that carries the peak-derived rows.",
        "label_classes": [
          "biophysical",
          "catalog"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": true
      },
      "applicability": {
        "k": {
          "min": 0.0,
          "max": 1015.0,
          "p05": 0.0,
          "p95": 110.0,
          "n_unique": 57
        },
        "na": {
          "min": 0.0,
          "max": 1020.0,
          "p05": 0.0,
          "p95": 100.0,
          "n_unique": 39
        },
        "li_nh4": {
          "min": 0.0,
          "max": 300.0,
          "p05": 0.0,
          "p95": 20.0,
          "n_unique": 19
        },
        "mg": {
          "min": 0.0,
          "max": 10.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 4
        },
        "ph": {
          "min": 4.0,
          "max": 8.0,
          "p05": 7.0,
          "p95": 7.5,
          "n_unique": 21
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 5.0,
          "max": 99.0,
          "p05": 16.0,
          "p95": 37.0,
          "n_unique": 59
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 0.2,
          "max": 300.0,
          "p05": 0.25,
          "p95": 6.0,
          "n_unique": 19
        }
      }
    },
    "g4_topology": {
      "name": "g4_topology",
      "kind": "G4",
      "task": "multiclass",
      "target": "topology",
      "classes": [
        "parallel",
        "antiparallel",
        "hybrid"
      ],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "balanced_accuracy": 0.7238345373127864,
        "accuracy": 0.7522388059701492,
        "auroc_ovr_macro": 0.9048679527696867,
        "brier_uncalibrated": 0.33802205527709006,
        "ece_uncalibrated": 0.07221034029823038,
        "temperature": 1.3457478771154054,
        "brier": 0.32784474673071606,
        "ece": 0.03587319431330903,
        "reliability": {
          "confidence": [
            0.45904413220251683,
            0.5461590004568619,
            0.619395954843809,
            0.6962883114520467,
            0.7784073422143305,
            0.8451655546637233,
            0.8985826596810991,
            0.9315641178084508,
            0.9557227661532725,
            0.9786459459755749
          ],
          "accuracy": [
            0.43564356435643564,
            0.5,
            0.5643564356435643,
            0.65,
            0.6831683168316832,
            0.85,
            0.93,
            0.9702970297029703,
            0.95,
            0.9900990099009901
          ],
          "count": [
            101,
            100,
            101,
            100,
            101,
            100,
            100,
            101,
            100,
            101
          ]
        },
        "per_class": {
          "parallel": {
            "n": 487,
            "recall": 0.864476386036961
          },
          "antiparallel": {
            "n": 216,
            "recall": 0.6944444444444444
          },
          "hybrid": {
            "n": 302,
            "recall": 0.6125827814569537
          }
        },
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE is measured over the three topology classes on sequences already known to form a G4, in K+ buffer. It says nothing about whether a sequence folds, and nothing about topology in Na+."
      },
      "training_meta": {
        "n_rows": 1005,
        "n_features": 160,
        "tiers": [
          "experimental"
        ],
        "sources": [
          "g4sp_topology"
        ],
        "grouping": {
          "n_items": 1005,
          "n_groups": 617,
          "largest_group": 126,
          "singleton_groups": 513,
          "mean_group_size": 1.6288492706645057
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 31.4,
        "claim": "Calibrated posterior over folding topology, given that the sequence forms a G4. Trained on K+ buffer data only.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 100.0,
          "max": 100.0,
          "p05": 100.0,
          "p95": 100.0,
          "n_unique": 1
        },
        "na": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 7.0,
          "max": 7.0,
          "p05": 7.0,
          "p95": 7.0,
          "n_unique": 1
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 13.0,
          "max": 72.0,
          "p05": 18.0,
          "p95": 36.0,
          "n_unique": 41
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "g4_tm": {
      "name": "g4_tm",
      "kind": "G4",
      "task": "regression",
      "target": "tm",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "r2": 0.6700589834856598,
        "rmse": 7.7564071253868825,
        "mae": 5.579288732809505,
        "spearman": 0.8163358638761452,
        "oof_residual_halfwidth": 12.361515045166016,
        "oof_residual_alpha": 0.1,
        "in_sample_coverage_of_residual_pool": 0.9001759014951627,
        "split_conformal_coverage": 0.8887951095753515,
        "split_conformal_coverage_sd": 0.026688967392108,
        "split_conformal_halfwidth_mean": 12.207139984130858,
        "split_conformal_n_splits": 20,
        "split_conformal_note": "half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90%",
        "interval_note": "interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool",
        "target_range": [
          7.7,
          94.0
        ],
        "calibration_scope": "Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale."
      },
      "training_meta": {
        "n_rows": 2274,
        "n_features": 160,
        "tiers": [
          "experimental"
        ],
        "sources": [
          "g4stab_experimental_tm"
        ],
        "grouping": {
          "n_items": 2274,
          "n_groups": 400,
          "largest_group": 305,
          "singleton_groups": 185,
          "mean_group_size": 5.685
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 36.2,
        "claim": "Melting temperature under the supplied buffer. Requires experimental Tm records (G4STAB Supplementary Table 1) in the atlas.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 0.0,
          "max": 1015.0,
          "p05": 0.0,
          "p95": 120.0,
          "n_unique": 56
        },
        "na": {
          "min": 0.0,
          "max": 1020.0,
          "p05": 0.0,
          "p95": 110.0,
          "n_unique": 39
        },
        "li_nh4": {
          "min": 0.0,
          "max": 300.0,
          "p05": 0.0,
          "p95": 30.0,
          "n_unique": 19
        },
        "mg": {
          "min": 0.0,
          "max": 10.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 4
        },
        "ph": {
          "min": 4.0,
          "max": 8.0,
          "p05": 7.0,
          "p95": 7.5,
          "n_unique": 21
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 5.0,
          "max": 99.0,
          "p05": 15.0,
          "p95": 37.34999999999991,
          "n_unique": 56
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 0.2,
          "max": 300.0,
          "p05": 0.25,
          "p95": 10.0,
          "n_unique": 19
        },
        "temperature_is_enforced": false,
        "temperature_note": "NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask."
      }
    },
    "g4_tm_distilled": {
      "name": "g4_tm_distilled",
      "kind": "G4",
      "task": "regression",
      "target": "tm",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 3,
      "metrics": {
        "r2": 0.9631506722490957,
        "rmse": 1.9067025378416678,
        "mae": 1.4765320090611775,
        "spearman": 0.9785367398600405,
        "oof_residual_halfwidth": 3.0373191833496094,
        "oof_residual_alpha": 0.1,
        "in_sample_coverage_of_residual_pool": 0.9000111111111111,
        "split_conformal_coverage": 0.8961118222343449,
        "split_conformal_coverage_sd": 0.01512553543979076,
        "split_conformal_halfwidth_mean": 3.029021759033204,
        "split_conformal_n_splits": 20,
        "split_conformal_note": "half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90%",
        "interval_note": "interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool",
        "target_range": [
          28.5,
          89.1
        ],
        "calibration_scope": "Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale."
      },
      "training_meta": {
        "n_rows": 90000,
        "n_features": 160,
        "tiers": [
          "predicted"
        ],
        "sources": [
          "g4stab_webdb_predicted"
        ],
        "grouping": {
          "n_items": 90000,
          "n_groups": 3552,
          "largest_group": 5355,
          "singleton_groups": 74,
          "mean_group_size": 25.33783783783784
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 3,
        "seconds": 308.6,
        "claim": "DISTILLED: reproduces the G4STAB ensemble's Tm response to cation composition. Useful as a condition-response prior and as a benchmark; it is NOT trained on measurements and inherits every G4STAB bias.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": false,
        "target_semantics": "predicted",
        "target_semantics_label": "prior-model output",
        "biophysically_grounded": false,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 0.0,
          "max": 100.0,
          "p05": 0.0,
          "p95": 100.0,
          "n_unique": 4
        },
        "na": {
          "min": 0.0,
          "max": 100.0,
          "p05": 0.0,
          "p95": 100.0,
          "n_unique": 4
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 7.0,
          "max": 7.0,
          "p05": 7.0,
          "p95": 7.0,
          "n_unique": 1
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 14.0,
          "max": 50.0,
          "p05": 17.0,
          "p95": 47.0,
          "n_unique": 37
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        },
        "temperature_is_enforced": false,
        "temperature_note": "NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask."
      }
    },
    "im_fold": {
      "name": "im_fold",
      "kind": "iM",
      "task": "binary",
      "target": "folded",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "auroc": 0.8704296874999998,
        "auprc": 0.8838150075413367,
        "brier_uncalibrated": 0.15321962330867572,
        "ece_uncalibrated": 0.0896309376868885,
        "calibration_method": "platt",
        "ece_isotonic": 0.0782731065761882,
        "ece_platt": 0.05524021333043322,
        "n_distinct_calibrated_values": 293,
        "brier": 0.14832809236733616,
        "ece": 0.05524021333043322,
        "mce": 0.12672396213581605,
        "accuracy_at_0.5": 0.7875,
        "reliability": {
          "confidence": [
            0.04647772341756462,
            0.09986555027139085,
            0.17686852052644292,
            0.29011583044054456,
            0.4200416615967406,
            0.5614410345406691,
            0.720473962135816,
            0.8248474000696289,
            0.9056554845331379,
            0.948504888363771
          ],
          "accuracy": [
            0.09375,
            0.09375,
            0.15625,
            0.40625,
            0.3125,
            0.59375,
            0.59375,
            0.8125,
            0.9375,
            1.0
          ],
          "count": [
            32,
            32,
            32,
            32,
            32,
            32,
            32,
            32,
            32,
            32
          ]
        },
        "positive_rate": 0.5,
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability."
      },
      "training_meta": {
        "n_rows": 320,
        "n_features": 160,
        "tiers": [
          "derived",
          "experimental"
        ],
        "sources": [
          "im_pht_biophysical",
          "shuffled::im_pht_biophysical"
        ],
        "grouping": {
          "n_items": 320,
          "n_groups": 175,
          "largest_group": 43,
          "singleton_groups": 137,
          "mean_group_size": 1.8285714285714285
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 4.8,
        "claim": "P(sequence forms an i-motif) against a composition-matched dinucleotide-shuffled background. Positives are sequences with a measured transitional pH, so this is grounded in folding measurements -- but the panel is small and heavily focused on designed C-tract variants.",
        "label_classes": [
          "biophysical",
          "catalog"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": true
      },
      "applicability": {
        "k": {
          "min": 100.0,
          "max": 100.0,
          "p05": 100.0,
          "p95": 100.0,
          "n_unique": 1
        },
        "na": {
          "min": 10.0,
          "max": 10.0,
          "p05": 10.0,
          "p95": 10.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 7.0,
          "max": 7.0,
          "p05": 7.0,
          "p95": 7.0,
          "n_unique": 1
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 15.0,
          "max": 116.0,
          "p05": 18.950000000000003,
          "p95": 57.0,
          "n_unique": 43
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "im_pht": {
      "name": "im_pht",
      "kind": "iM",
      "task": "regression",
      "target": "ph_t",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "r2": 0.5915665061450279,
        "rmse": 0.3494439080174955,
        "mae": 0.25486612737178804,
        "spearman": 0.787842313426595,
        "oof_residual_halfwidth": 0.6317110061645508,
        "oof_residual_alpha": 0.1,
        "in_sample_coverage_of_residual_pool": 0.9125,
        "split_conformal_coverage": 0.9002130237501564,
        "split_conformal_coverage_sd": 0.05987825828267437,
        "split_conformal_halfwidth_mean": 0.6139051628112793,
        "split_conformal_n_splits": 20,
        "split_conformal_note": "half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90%",
        "interval_note": "interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool",
        "target_range": [
          5.1,
          7.9
        ],
        "calibration_scope": "Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale."
      },
      "training_meta": {
        "n_rows": 160,
        "n_features": 160,
        "tiers": [
          "experimental"
        ],
        "sources": [
          "im_pht_biophysical"
        ],
        "grouping": {
          "n_items": 160,
          "n_groups": 85,
          "largest_group": 19,
          "singleton_groups": 58,
          "mean_group_size": 1.8823529411764706
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 4.3,
        "claim": "Transitional pH from sequence, at the iM-Seeker reference buffer (100 mM KCl + 10 mM Na cacodylate). Answers 'which of these sequences folds at higher pH?'. Grouped by sequence similarity, so the number is what to expect on a sequence you have not measured. It carries no salt response, because every row was measured in one buffer.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 100.0,
          "max": 100.0,
          "p05": 100.0,
          "p95": 100.0,
          "n_unique": 1
        },
        "na": {
          "min": 10.0,
          "max": 10.0,
          "p05": 10.0,
          "p95": 10.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 7.0,
          "max": 7.0,
          "p05": 7.0,
          "p95": 7.0,
          "n_unique": 1
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 15.0,
          "max": 116.0,
          "p05": 18.95,
          "p95": 57.0,
          "n_unique": 43
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "im_pht_condition": {
      "name": "im_pht_condition",
      "kind": "iM",
      "task": "regression",
      "target": "ph_t",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "r2": 0.8345546834800903,
        "rmse": 0.18905261347485044,
        "mae": 0.13982849982477003,
        "spearman": 0.8752572669536914,
        "oof_residual_halfwidth": 0.29415418624877887,
        "oof_residual_alpha": 0.1,
        "in_sample_coverage_of_residual_pool": 0.9014336917562724,
        "split_conformal_coverage": 0.880494374958045,
        "split_conformal_coverage_sd": 0.07157546679808792,
        "split_conformal_halfwidth_mean": 0.2958197174072265,
        "split_conformal_n_splits": 20,
        "split_conformal_note": "half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90%",
        "interval_note": "interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool",
        "target_range": [
          5.13,
          7.0
        ],
        "calibration_scope": "Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale."
      },
      "training_meta": {
        "n_rows": 558,
        "n_features": 160,
        "tiers": [
          "experimental"
        ],
        "sources": [
          "duvma_im_stability_landscape"
        ],
        "grouping": {
          "n_items": 558,
          "n_groups": 10,
          "largest_group": 68,
          "singleton_groups": 0,
          "mean_group_size": 55.8
        },
        "group_by": "condition",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 4.2,
        "claim": "Transitional pH as a function of buffer, for the two constructs 5DUVMA measured across ten ionic strengths. GROUPED BY BUFFER, not by sequence: folds hold out whole ionic strengths, so the number answers 'does this extrapolate to a salt concentration we did not measure?'. It says nothing about a new sequence -- there are only two here. REFUSES ANY SEQUENCE OUTSIDE ITS TWO-CONSTRUCT ALLOWLIST: no value is returned, because this head gives essentially the same answer whatever sequence it is given, so there is no prediction to report rather than an uncertain one.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 10.0,
          "max": 1000.0,
          "p05": 10.0,
          "p95": 500.0,
          "n_unique": 10
        },
        "na": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 7.0,
          "max": 7.0,
          "p05": 7.0,
          "p95": 7.0,
          "n_unique": 1
        },
        "temperature": {
          "min": 30.0,
          "max": 60.0,
          "p05": 33.0,
          "p95": 58.0,
          "n_unique": 31
        },
        "length": {
          "min": 21.0,
          "max": 45.0,
          "p05": 21.0,
          "p95": 45.0,
          "n_unique": 2
        },
        "sequence_allowlist": [
          "CCCCCCCCCTTTCCCCCCCCCTTTCCCCCCCCCTTTCCCCCCCCC",
          "CCCTAACCCTAACCCTAACCC"
        ],
        "sequence_allowlist_note": "This head was trained on 2 distinct sequences. It learned how those constructs respond to conditions, not how sequence affects the response. Any other sequence is out of domain and is refused.",
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "im_tm_condition": {
      "name": "im_tm_condition",
      "kind": "iM",
      "task": "regression",
      "target": "tm",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "r2": 0.8531538765417034,
        "rmse": 6.458634837633418,
        "mae": 4.776799930381272,
        "spearman": 0.9264658352406328,
        "oof_residual_halfwidth": 11.180231170654295,
        "oof_residual_alpha": 0.1,
        "in_sample_coverage_of_residual_pool": 0.9023746701846965,
        "split_conformal_coverage": 0.9043043951728252,
        "split_conformal_coverage_sd": 0.04289727453470009,
        "split_conformal_halfwidth_mean": 11.392365058898926,
        "split_conformal_n_splits": 20,
        "split_conformal_note": "half the groups calibrate the half-width, the other half measure coverage; mean over 20 random group partitions at nominal 90%",
        "interval_note": "interval half-width = oof_residual_halfwidth; coverage to quote = split_conformal_coverage (held-out groups), not in_sample_coverage_of_residual_pool",
        "target_range": [
          8.05,
          80.37
        ],
        "calibration_scope": "Regression heads carry an OOF-residual interval, not a conformal guarantee: the half-width and its reported coverage come from the same pool of out-of-fold residuals, so the coverage figure is in-sample to that pool. Treat the width as a typical error scale."
      },
      "training_meta": {
        "n_rows": 379,
        "n_features": 160,
        "tiers": [
          "experimental"
        ],
        "sources": [
          "duvma_im_stability_landscape"
        ],
        "grouping": {
          "n_items": 379,
          "n_groups": 10,
          "largest_group": 40,
          "singleton_groups": 0,
          "mean_group_size": 37.9
        },
        "group_by": "condition",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 4.9,
        "claim": "i-motif melting temperature as a function of pH and ionic strength, for the two 5DUVMA constructs. Grouped by buffer, same caveat as im_pht_condition: a condition-response model, not a sequence model. REFUSES ANY SEQUENCE OUTSIDE ITS TWO-CONSTRUCT ALLOWLIST: no value is returned, because this head gives essentially the same answer whatever sequence it is given, so there is no prediction to report rather than an uncertain one.",
        "label_classes": [
          "biophysical"
        ],
        "has_experimental_observation": true,
        "target_semantics": "biophysical",
        "target_semantics_label": "biophysically anchored",
        "biophysically_grounded": true,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 10.0,
          "max": 1000.0,
          "p05": 10.0,
          "p95": 1000.0,
          "n_unique": 10
        },
        "na": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 4.35,
          "max": 7.05,
          "p05": 4.65,
          "p95": 7.05,
          "n_unique": 19
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 21.0,
          "max": 45.0,
          "p05": 21.0,
          "p95": 45.0,
          "n_unique": 2
        },
        "sequence_allowlist": [
          "CCCCCCCCCTTTCCCCCCCCCTTTCCCCCCCCCTTTCCCCCCCCC",
          "CCCTAACCCTAACCCTAACCC"
        ],
        "sequence_allowlist_note": "This head was trained on 2 distinct sequences. It learned how those constructs respond to conditions, not how sequence affects the response. Any other sequence is out of domain and is refused.",
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        },
        "temperature_is_enforced": false,
        "temperature_note": "NOT a domain gate for this head: it PREDICTS a melting temperature and does not take one as an input. Every training row records temperature=25 as a metadata default on a Tm measurement, so this range describes the records rather than a limit on what you may ask."
      }
    },
    "im_fold_genomic": {
      "name": "im_fold_genomic",
      "kind": "iM",
      "task": "binary",
      "target": "folded",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "auroc": 0.7789806779099697,
        "auprc": 0.7903197702092515,
        "brier_uncalibrated": 0.19140602442221702,
        "ece_uncalibrated": 0.023003698081160134,
        "calibration_method": "platt",
        "ece_isotonic": 0.005556679489452704,
        "ece_platt": 0.006835794297730338,
        "n_distinct_calibrated_values": 7476,
        "brier": 0.1908482612007066,
        "ece": 0.006835794297730338,
        "mce": 0.018309717299423656,
        "accuracy_at_0.5": 0.7035848707348676,
        "reliability": {
          "confidence": [
            0.09690198555563032,
            0.2283794157117488,
            0.33802025900521243,
            0.4284266287074031,
            0.504075145821945,
            0.5771103229535514,
            0.6473451876236961,
            0.7189270012500409,
            0.7970469576800072,
            0.897510067879713
          ],
          "accuracy": [
            0.09870450339296731,
            0.243059839605182,
            0.3395061728395062,
            0.42751388032078963,
            0.4959901295496607,
            0.5753086419753086,
            0.6465144972239358,
            0.7006172839506173,
            0.7927205428747687,
            0.9136335595311537
          ],
          "count": [
            1621,
            1621,
            1620,
            1621,
            1621,
            1620,
            1621,
            1620,
            1621,
            1621
          ]
        },
        "positive_rate": 0.5233541062503856,
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability."
      },
      "training_meta": {
        "n_rows": 16207,
        "n_features": 160,
        "tiers": [
          "derived",
          "experimental"
        ],
        "sources": [
          "gse220882_hek293t_im_cutandtag",
          "shuffled::gse220882_hek293t_im_cutandtag"
        ],
        "grouping": {
          "n_items": 16207,
          "n_groups": 3907,
          "largest_group": 50,
          "singleton_groups": 2216,
          "mean_group_size": 4.148195546455081
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 451.7,
        "claim": "GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. Calibrated score for separating canonical i-motif motifs found in reproducible iMab CUT&Tag peaks (live HEK293T, GSE220882) from composition-matched motifs taken from dinucleotide shuffles of the same peak windows. The positive label is antibody occupancy at a locus, not a measured folding event, and no buffer was measured -- the attached condition is a nominal guess with every field flagged imputed. Read it as 'does this look like the motifs enriched under iMab peaks', and never as P(folds). Its value over im_fold is sequence diversity: thousands of independent loci instead of variants of a handful of designed constructs.",
        "label_classes": [
          "catalog",
          "genomic_proxy"
        ],
        "has_experimental_observation": true,
        "target_semantics": "genomic_proxy",
        "target_semantics_label": "genomic proxy",
        "biophysically_grounded": false,
        "mixed_label_classes": true
      },
      "applicability": {
        "k": {
          "min": 140.0,
          "max": 140.0,
          "p05": 140.0,
          "p95": 140.0,
          "n_unique": 1
        },
        "na": {
          "min": 10.0,
          "max": 10.0,
          "p05": 10.0,
          "p95": 10.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.5,
          "max": 0.5,
          "p05": 0.5,
          "p95": 0.5,
          "n_unique": 1
        },
        "ph": {
          "min": 7.2,
          "max": 7.2,
          "p05": 7.2,
          "p95": 7.2,
          "n_unique": 1
        },
        "temperature": {
          "min": 37.0,
          "max": 37.0,
          "p05": 37.0,
          "p95": 37.0,
          "n_unique": 1
        },
        "length": {
          "min": 15.0,
          "max": 48.0,
          "p05": 22.0,
          "p95": 44.0,
          "n_unique": 34
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "g4_fold_genomic": {
      "name": "g4_fold_genomic",
      "kind": "G4",
      "task": "binary",
      "target": "folded",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "auroc": 0.7645305501823564,
        "auprc": 0.7972703234006233,
        "brier_uncalibrated": 0.19163715820047061,
        "ece_uncalibrated": 0.03826803122091427,
        "calibration_method": "isotonic",
        "ece_isotonic": 0.008255371042776504,
        "ece_platt": 0.021369224020754383,
        "n_distinct_calibrated_values": 257,
        "brier": 0.19032557865914607,
        "ece": 0.008255371042776504,
        "mce": 0.014027548564128245,
        "accuracy_at_0.5": 0.7061903555442952,
        "reliability": {
          "confidence": [
            0.08112173563447969,
            0.25475340624313186,
            0.4689469425982132,
            0.5585316802245776,
            0.5971298724485355,
            0.6166464843137504,
            0.6708279431669092,
            0.744796779333359,
            0.7939356441602885,
            0.9218879832313092
          ],
          "accuracy": [
            0.08666666666666667,
            0.2606516290726817,
            0.48006644518272423,
            0.5713159968479118,
            0.5908761766835626,
            0.6069114470842333,
            0.6741479634247715,
            0.7307692307692307,
            0.796718322698268,
            0.9110091743119266
          ],
          "count": [
            1200,
            1197,
            1204,
            1269,
            1381,
            926,
            1203,
            1274,
            1097,
            1090
          ]
        },
        "positive_rate": 0.5663373025926864,
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability."
      },
      "training_meta": {
        "n_rows": 11841,
        "n_features": 160,
        "tiers": [
          "derived",
          "experimental"
        ],
        "sources": [
          "gse220882_hek293t_g4_cutandtag",
          "shuffled::gse220882_hek293t_g4_cutandtag"
        ],
        "grouping": {
          "n_items": 11841,
          "n_groups": 2929,
          "largest_group": 37,
          "singleton_groups": 1192,
          "mean_group_size": 4.042676681461249
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 127.4,
        "claim": "GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. The BG4 counterpart of im_fold_genomic, on the same peaks and the same matched-shuffle control. It exists mainly as a control: the two heads share a pipeline, so a gap between them says something about the antibodies and the underlying motif grammar rather than about the model.",
        "label_classes": [
          "catalog",
          "genomic_proxy"
        ],
        "has_experimental_observation": true,
        "target_semantics": "genomic_proxy",
        "target_semantics_label": "genomic proxy",
        "biophysically_grounded": false,
        "mixed_label_classes": true
      },
      "applicability": {
        "k": {
          "min": 140.0,
          "max": 140.0,
          "p05": 140.0,
          "p95": 140.0,
          "n_unique": 1
        },
        "na": {
          "min": 10.0,
          "max": 10.0,
          "p05": 10.0,
          "p95": 10.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.5,
          "max": 0.5,
          "p05": 0.5,
          "p95": 0.5,
          "n_unique": 1
        },
        "ph": {
          "min": 7.2,
          "max": 7.2,
          "p05": 7.2,
          "p95": 7.2,
          "n_unique": 1
        },
        "temperature": {
          "min": 37.0,
          "max": 37.0,
          "p05": 37.0,
          "p95": 37.0,
          "n_unique": 1
        },
        "length": {
          "min": 15.0,
          "max": 33.0,
          "p05": 19.0,
          "p95": 31.0,
          "n_unique": 19
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "im_architecture": {
      "name": "im_architecture",
      "kind": "iM",
      "task": "binary",
      "target": "folded",
      "classes": [],
      "use_conditions": true,
      "n_estimators": 3,
      "metrics": {
        "auroc": 0.9966078880491571,
        "auprc": 0.9972920081585112,
        "brier_uncalibrated": 0.013121124587354881,
        "ece_uncalibrated": 0.007679089951485591,
        "calibration_method": "platt",
        "ece_isotonic": 0.005224470823338512,
        "ece_platt": 0.006687820450115305,
        "n_distinct_calibrated_values": 1505,
        "brier": 0.013027352402325096,
        "ece": 0.006687820450115305,
        "mce": 0.027001530923803127,
        "accuracy_at_0.5": 0.9844375,
        "reliability": {
          "confidence": [
            0.0015475397615119125,
            0.0031314132768432905,
            0.005498469076196875,
            0.012177555820529748,
            0.06072087906595998,
            0.9312730727560541,
            0.9969559149534402,
            0.9980278369105062,
            0.9985962612256413,
            0.9991451139164078
          ],
          "accuracy": [
            0.0,
            0.0025,
            0.0325,
            0.00875,
            0.04625,
            0.91875,
            1.0,
            1.0,
            1.0,
            1.0
          ],
          "count": [
            1600,
            1600,
            1600,
            1600,
            1600,
            1600,
            1602,
            1598,
            1600,
            1600
          ]
        },
        "positive_rate": 0.500875,
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE and Brier are measured on the constructed, approximately balanced positive-versus-matched-shuffle task this head was trained on. The negatives are generated, not measured non-folders, so these are not absolute probabilities for an arbitrary genomic, transcriptomic or aptamer population. Recalibrate against the prevalence of your own candidate set before treating a number as a probability."
      },
      "training_meta": {
        "n_rows": 16000,
        "n_features": 160,
        "tiers": [
          "derived"
        ],
        "sources": [
          "im_architecture_negatives",
          "im_architecture_positives"
        ],
        "grouping": {
          "n_items": 16000,
          "n_groups": 3849,
          "largest_group": 74,
          "singleton_groups": 1976,
          "mean_group_size": 4.156923876331515
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 3,
        "seconds": 122.3,
        "claim": "ARCHITECTURE ONLY, AND NEARLY DETERMINISTIC BY CONSTRUCTION. The positive label is 'carries a canonical four-C-tract motif' and the negatives are shuffles from which that motif was rejected, so a near-perfect AUROC here means the features can recompute the regex -- it is NOT evidence of biological predictive power. This head is a placeholder that exists so the iM pipeline is exercisable end to end; im_fold and im_fold_genomic are the heads that now do the work it was standing in for.",
        "label_classes": [
          "catalog"
        ],
        "has_experimental_observation": false,
        "target_semantics": "derived",
        "target_semantics_label": "derived / control",
        "biophysically_grounded": false,
        "mixed_label_classes": false
      },
      "applicability": {
        "k": {
          "min": 100.0,
          "max": 100.0,
          "p05": 100.0,
          "p95": 100.0,
          "n_unique": 1
        },
        "na": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "ph": {
          "min": 5.8,
          "max": 5.8,
          "p05": 5.8,
          "p95": 5.8,
          "n_unique": 1
        },
        "temperature": {
          "min": 25.0,
          "max": 25.0,
          "p05": 25.0,
          "p95": 25.0,
          "n_unique": 1
        },
        "length": {
          "min": 15.0,
          "max": 50.0,
          "p05": 19.0,
          "p95": 45.0,
          "n_unique": 36
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    },
    "locus_peak_overlap_state": {
      "name": "locus_peak_overlap_state",
      "kind": "locus",
      "task": "multiclass",
      "target": "topology",
      "classes": [
        "both",
        "g4_only",
        "im_only",
        "neither"
      ],
      "use_conditions": true,
      "n_estimators": 5,
      "metrics": {
        "balanced_accuracy": 0.45035592771429633,
        "accuracy": 0.6143455853202696,
        "auroc_ovr_macro": 0.800048554902052,
        "brier_uncalibrated": 0.5007924640220067,
        "ece_uncalibrated": 0.01622351984809746,
        "temperature": 0.9369423542506078,
        "brier": 0.5004025621747314,
        "ece": 0.010126140874070206,
        "reliability": {
          "confidence": [
            0.3525297931190548,
            0.40527687264829093,
            0.44833235187596654,
            0.4941951054626685,
            0.5482845859107549,
            0.6148253302057357,
            0.6960070327656452,
            0.7831615815007273,
            0.8706228116103385,
            0.9518534818539682
          ],
          "accuracy": [
            0.342347444428713,
            0.39331728727169757,
            0.4454197932889707,
            0.4872557349192863,
            0.5254141299730992,
            0.6082401245929492,
            0.700792976493911,
            0.7948463825569871,
            0.8847515220161405,
            0.9610647033838312
          ],
          "count": [
            7063,
            7063,
            7063,
            7062,
            7063,
            7063,
            7062,
            7063,
            7063,
            7063
          ]
        },
        "per_class": {
          "both": {
            "n": 13156,
            "recall": 0.3920644572818486
          },
          "g4_only": {
            "n": 15088,
            "recall": 0.35319459172852596
          },
          "im_only": {
            "n": 7070,
            "recall": 0.15558698727015557
          },
          "neither": {
            "n": 35314,
            "recall": 0.9005776745766552
          }
        },
        "calibration_note": "calibrated metrics are cross-fitted across sequence groups",
        "calibration_scope": "ECE is measured over the three topology classes on sequences already known to form a G4, in K+ buffer. It says nothing about whether a sequence folds, and nothing about topology in Na+."
      },
      "training_meta": {
        "n_rows": 70628,
        "n_features": 182,
        "tiers": [
          "derived",
          "experimental"
        ],
        "sources": [
          "gse220882_hek293t_locus_state",
          "shuffled::gse220882_hek293t_locus_state"
        ],
        "grouping": {
          "n_items": 70628,
          "n_groups": 3940,
          "largest_group": 624,
          "singleton_groups": 3076,
          "mean_group_size": 17.925888324873096
        },
        "group_by": "sequence",
        "cluster_threshold": 0.9,
        "n_folds": 5,
        "n_seeds": 5,
        "seconds": 1073.2,
        "claim": "GENOMIC PEAK-OVERLAP PROXY FOR A 201-NT WINDOW. NOT CO-OCCUPANCY AND NOT COMPETITION. Calibrated posterior over how a 201-nt genomic window was classified by two CUT&Tag experiments run as PARALLEL REACTIONS on separate aliquots of the same HEK293T population: reproducible peaks from both antibodies, from BG4 only, from iMab only, or from neither (a composition-matched shuffle of the same windows -- generated, not an observed unoccupied locus). Because the two antibodies never touched the same aliquot, a 'both' call is population-level co-localisation of peaks; nothing here observes two structures on one DNA molecule, and nothing here is a free energy. Balanced accuracy 0.450 against a 0.250 four-class floor, iM-only recall 0.156 -- the class most users are asking about is the one it recovers worst. Trained exclusively on 201-nt windows, so a query of oligonucleotide length is outside the applicability domain by construction. Read it as an exploratory genomic-resemblance score.",
        "label_classes": [
          "catalog",
          "genomic_proxy"
        ],
        "has_experimental_observation": true,
        "target_semantics": "genomic_proxy",
        "target_semantics_label": "genomic proxy",
        "biophysically_grounded": false,
        "mixed_label_classes": true,
        "renamed_from": "locus_state",
        "window_nt": 201
      },
      "applicability": {
        "k": {
          "min": 140.0,
          "max": 140.0,
          "p05": 140.0,
          "p95": 140.0,
          "n_unique": 1
        },
        "na": {
          "min": 10.0,
          "max": 10.0,
          "p05": 10.0,
          "p95": 10.0,
          "n_unique": 1
        },
        "li_nh4": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "mg": {
          "min": 0.5,
          "max": 0.5,
          "p05": 0.5,
          "p95": 0.5,
          "n_unique": 1
        },
        "ph": {
          "min": 7.2,
          "max": 7.2,
          "p05": 7.2,
          "p95": 7.2,
          "n_unique": 1
        },
        "temperature": {
          "min": 37.0,
          "max": 37.0,
          "p05": 37.0,
          "p95": 37.0,
          "n_unique": 1
        },
        "length": {
          "min": 201.0,
          "max": 201.0,
          "p05": 201.0,
          "p95": 201.0,
          "n_unique": 1
        },
        "crowder_pct": {
          "min": 0.0,
          "max": 0.0,
          "p05": 0.0,
          "p95": 0.0,
          "n_unique": 1
        },
        "strand_conc": {
          "min": 5.0,
          "max": 5.0,
          "p05": 5.0,
          "p95": 5.0,
          "n_unique": 1
        }
      }
    }
  }
}
```
</details>