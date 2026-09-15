"""Training: grouped cross-validation, out-of-fold calibration, conformal intervals.

The procedure for every head is the same and is deliberately conservative:

1. Pull labelled rows from the atlas at the requested evidence tiers.
2. Cluster sequences (``grouping.cluster_sequences``) and build grouped folds so
   near-duplicates cannot straddle a split.
3. Produce out-of-fold predictions.  All reported metrics come from these.
4. Fit the calibrator on the out-of-fold predictions -- never on training-fold
   predictions, which are over-confident by construction.
5. Refit the base ensemble on all data; keep the OOF-fitted calibrator.
6. For regression, derive a split-conformal half-width from OOF residuals, so
   the reported interval has (marginal) finite-sample coverage.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

from .. import claims
from ..conditions import Condition
from ..features import featurize, feature_names
from .base import Head, MultiTaskModel, rows_fingerprint
from .calibration import (
    BinaryCalibrator, TemperatureScaler, brier_multiclass,
    expected_calibration_error, maximum_calibration_error, multiclass_ece,
    reliability_curve,
)
from .grouping import cluster_sequences, group_report


@dataclass
class TaskSpec:
    """One head to train.

    ``tiers`` is the contract that keeps the model card honest: a head whose
    tiers do not include ``experimental`` is reported everywhere as not
    grounded in measurements.
    """

    name: str
    kind: str
    task: str
    target: str
    classes: Sequence[str] = ()
    tiers: Sequence[str] = ("experimental",)
    min_rows: int = 40
    max_rows: int | None = None
    n_seeds: int | None = None
    claim: str = ""
    sources: Sequence[str] | None = None
    label_classes: Sequence[str] | None = None
    group_by: str = "sequence"


DEFAULT_TASKS = [
    TaskSpec(
        "g4_fold", "G4", "binary", "folded",
        tiers=("experimental", "derived"),
        label_classes=("biophysical", "catalog"),
        sources=("g4stab_experimental_tm", "shuffled::g4stab_experimental_tm",
                 "g4sp_topology", "shuffled::g4sp_topology"),
        claim="P(sequence adopts a G-quadruplex) against a composition-matched "
              "dinucleotide-shuffled background. Pinned to its four biophysical sources "
              "by name. Filtering on label_class alone was not enough: the GSE220882 "
              "negatives are label_class=catalog like every other shuffle, so they "
              "passed the filter while their positives did not, and the head silently "
              "gained 5,135 unpaired negatives and an AUROC to match. "
              "g4_fold_genomic is the head that carries the peak-derived rows.",
    ),
    TaskSpec(
        "g4_topology", "G4", "multiclass", "topology",
        classes=("parallel", "antiparallel", "hybrid"),
        claim="Calibrated posterior over folding topology, given that the sequence "
              "forms a G4. Trained on K+ buffer data only.",
    ),
    TaskSpec(
        "g4_tm", "G4", "regression", "tm",
        claim="Melting temperature under the supplied buffer. Requires experimental "
              "Tm records (G4STAB Supplementary Table 1) in the atlas.",
    ),
    TaskSpec(
        "g4_tm_distilled", "G4", "regression", "tm",
        tiers=("predicted",), max_rows=90_000, n_seeds=3,
        claim="DISTILLED: reproduces the G4STAB ensemble's Tm response to cation "
              "composition. Useful as a condition-response prior and as a benchmark; "
              "it is NOT trained on measurements and inherits every G4STAB bias.",
    ),
    TaskSpec(
        "im_fold", "iM", "binary", "folded",
        tiers=("experimental", "derived"),
        sources=("im_pht_biophysical", "shuffled::im_pht_biophysical"),
        label_classes=("biophysical", "catalog"),
        claim="P(sequence forms an i-motif) against a composition-matched "
              "dinucleotide-shuffled background. Positives are sequences with a measured "
              "transitional pH, so this is grounded in folding measurements -- but the "
              "panel is small and heavily focused on designed C-tract variants.",
    ),
    TaskSpec(
        "im_pht", "iM", "regression", "ph_t",
        label_classes=("biophysical",), sources=("im_pht_biophysical",),
        claim="Transitional pH from sequence, at the iM-Seeker reference buffer "
              "(100 mM KCl + 10 mM Na cacodylate). Answers 'which of these sequences "
              "folds at higher pH?'. Grouped by sequence similarity, so the number is "
              "what to expect on a sequence you have not measured. It carries no salt "
              "response, because every row was measured in one buffer.",
    ),
    TaskSpec(
        "im_pht_condition", "iM", "regression", "ph_t",
        label_classes=("biophysical",), sources=("duvma_im_stability_landscape",),
        group_by="condition",
        claim="Transitional pH as a function of buffer, for the two constructs 5DUVMA "
              "measured across ten ionic strengths. GROUPED BY BUFFER, not by sequence: "
              "folds hold out whole ionic strengths, so the number answers 'does this "
              "extrapolate to a salt concentration we did not measure?'. It says nothing "
              "about a new sequence -- there are only two here.",
    ),
    TaskSpec(
        "im_tm_condition", "iM", "regression", "tm",
        label_classes=("biophysical",), sources=("duvma_im_stability_landscape",),
        group_by="condition",
        claim="i-motif melting temperature as a function of pH and ionic strength, for "
              "the two 5DUVMA constructs. Grouped by buffer, same caveat as "
              "im_pht_condition: a condition-response model, not a sequence model.",
    ),
    TaskSpec(
        "im_fold_genomic", "iM", "binary", "folded",
        tiers=("experimental", "derived"),
        sources=("gse220882_hek293t_im_cutandtag",
                 "shuffled::gse220882_hek293t_im_cutandtag"),
        label_classes=("genomic_proxy", "catalog"),
        claim="P(sequence forms an i-motif in cells), trained on canonical motifs drawn "
              "from reproducible iMab CUT&Tag peaks in live HEK293T cells (GSE220882) "
              "against composition-matched motifs from dinucleotide shuffles of the same "
              "peak windows. GENOMIC PROXY: the positive label is antibody occupancy at a "
              "locus, not a folding measurement, and no condition was measured -- this "
              "head must not be read as biophysics. What it adds over im_fold is sequence "
              "diversity: thousands of independent loci instead of variants of a handful "
              "of designed constructs.",
    ),
    TaskSpec(
        "g4_fold_genomic", "G4", "binary", "folded",
        tiers=("experimental", "derived"),
        sources=("gse220882_hek293t_g4_cutandtag",
                 "shuffled::gse220882_hek293t_g4_cutandtag"),
        label_classes=("genomic_proxy", "catalog"),
        claim="The BG4 counterpart of im_fold_genomic, on the same peaks and the same "
              "matched-shuffle control. It exists mainly as a control: the two heads "
              "share a pipeline, so a gap between them says something about the "
              "antibodies and the underlying motif grammar rather than about the model.",
    ),
    TaskSpec(
        "locus_peak_overlap_state", "locus", "multiclass", "topology",
        classes=("both", "g4_only", "im_only", "neither"),
        tiers=("experimental", "derived"),
        sources=("gse220882_hek293t_locus_state",
                 "shuffled::gse220882_hek293t_locus_state"),
        label_classes=("genomic_proxy", "catalog"),
        claim="GENOMIC PEAK-OVERLAP PROXY FOR A 201-NT WINDOW. NOT CO-OCCUPANCY "
              "AND NOT COMPETITION. BG4 and iMab CUT&Tag were run as PARALLEL "
              "REACTIONS on separate aliquots of one HEK293T population, so a "
              "'both' call means both antibodies gave reproducible peaks over "
              "this window in that population -- nothing in the design observes "
              "two structures on the same DNA molecule. The 'neither' class is a "
              "composition-matched shuffle, not an observed unoccupied locus. "
              "Balanced accuracy 0.450 against a 0.250 floor; iM-only recall "
              "0.156. Every training window is 201 nt, so an oligo-length query "
              "is out of domain by construction. An exploratory genomic-"
              "resemblance score, and no more than that.",
    ),
    TaskSpec(
        "im_architecture", "iM", "binary", "folded",
        tiers=("derived",), max_rows=16_000, n_seeds=3,
        sources=("im_architecture_positives", "im_architecture_negatives"),
        claim="ARCHITECTURE ONLY, AND NEARLY DETERMINISTIC BY CONSTRUCTION. The positive "
              "label is 'carries a canonical four-C-tract motif' and the negatives are "
              "shuffles from which that motif was rejected, so a near-perfect AUROC here "
              "means the features can recompute the regex -- it is NOT evidence of "
              "biological predictive power. This head is a placeholder that exists so the "
              "iM pipeline is exercisable end to end; im_fold and im_fold_genomic "
              "are the heads that now do the work it was standing in for.",
    ),
]


def _estimator(task: str, seed: int, n_rows: int):
    from xgboost import XGBClassifier, XGBRegressor

    common = dict(
        n_estimators=400 if n_rows > 300 else 250,
        learning_rate=0.05,
        max_depth=4 if n_rows < 1500 else 6,
        subsample=0.85,
        colsample_bytree=0.7,
        min_child_weight=2,
        reg_lambda=2.0,
        random_state=seed,
        n_jobs=2,
        tree_method="hist",
    )
    if task == "regression":
        return XGBRegressor(objective="reg:squarederror", **common)
    return XGBClassifier(eval_metric="mlogloss", **common)


def _rows_to_arrays(rows, spec: TaskSpec, use_conditions: bool):
    seqs, conds, ys, meta = [], [], [], []
    for r in rows:
        val = r[spec.target]
        if val is None:
            continue
        if spec.task == "multiclass":
            if val not in spec.classes:
                continue
            y = list(spec.classes).index(val)
        elif spec.task == "binary":
            y = int(val)
        else:
            y = float(val)
        seqs.append(r["sequence"])
        conds.append(
            Condition.from_mapping(
                {c: r[c] for c in ("k", "na", "li_nh4", "mg", "ph", "temperature",
                                   "crowder_pct", "strand_conc")},
                track_imputed=False,
            )
        )
        ys.append(y)
        meta.append({"source": r["source"], "tier": r["evidence_tier"], "id": r["record_id"]})
    X = featurize(seqs, conds, kind=spec.kind, use_conditions=use_conditions)
    return seqs, conds, X, np.asarray(ys), meta


def _applicability(conds: list[Condition], seqs: list[str]) -> dict:
    def rng(vals):
        a = np.asarray(vals, dtype=float)
        return {"min": float(a.min()), "max": float(a.max()),
                "p05": float(np.percentile(a, 5)), "p95": float(np.percentile(a, 95)),
                "n_unique": int(len(np.unique(np.round(a, 3))))}

    lens = [len(s) for s in seqs]
    return {
        "k": rng([c.k for c in conds]),
        "na": rng([c.na for c in conds]),
        "li_nh4": rng([c.li_nh4 for c in conds]),
        "mg": rng([c.mg for c in conds]),
        "ph": rng([c.ph for c in conds]),
        "temperature": rng([c.temperature for c in conds]),
        # Present on every Condition, and omitted here until v0.4.3 -- so the
        # domain check had nothing to compare a query against and 40% PEG or
        # 500 uM strand passed as in-domain against training data that never
        # left 0% and 5 uM. A field the model carries but the applicability
        # table does not is a field nobody is checking.
        "crowder_pct": rng([c.crowder_pct for c in conds]),
        "strand_conc": rng([c.strand_conc for c in conds]),
        "length": rng(lens),
    }



def _crossfit_calibrated(oof, y, groups, task, n_classes, n_folds=5, method="isotonic"):
    """Honest calibrated predictions.

    Fitting a calibrator on the out-of-fold predictions and then scoring those
    same predictions understates ECE, sometimes to zero.  So we cross-fit the
    calibrator too: for each held-out group fold, fit on the rest and transform
    the held-out part.  The returned array is what the calibration metrics are
    computed on; the deployed calibrator is separately fitted on all OOF data.
    """
    from sklearn.model_selection import GroupKFold

    n = len(y)
    folds = min(n_folds, max(2, len(np.unique(groups))))
    out = np.empty((n, n_classes), dtype=float)
    for tr, te in GroupKFold(n_splits=folds).split(np.zeros(n), y, groups):
        if task == "binary":
            cal = BinaryCalibrator(method).fit(oof[tr, 1], y[tr])
            p = cal.transform(oof[te, 1])
            out[te, 1] = p
            out[te, 0] = 1 - p
        else:
            cal = TemperatureScaler().fit(oof[tr], y[tr])
            out[te] = cal.transform(oof[te])
    return out


def _condition_groups(conds) -> np.ndarray:
    """Fold groups defined by buffer rather than by sequence.

    For a head whose job is the *condition response* of a small number of
    constructs, grouping by sequence is meaningless -- there are two groups, and
    holding one out asks the model to predict a sequence it has never seen,
    which is a different question. Grouping by buffer instead holds out whole
    ionic strengths, so the reported number answers the question the head is
    actually for: does this extrapolate to a salt concentration we did not
    measure?
    """
    keys = [
        (round(c.k, 3), round(c.na, 3), round(c.li_nh4, 3), round(c.mg, 3))
        for c in conds
    ]
    uniq = {k: i for i, k in enumerate(sorted(set(keys)))}
    return np.array([uniq[k] for k in keys])



def _grouped_conformal(
    resid: np.ndarray,
    groups: np.ndarray,
    *,
    alpha: float = 0.1,
    n_splits: int = 20,
    seed: int = 0,
) -> dict:
    """Held-out coverage of an OOF-residual interval, split by group.

    Half the groups calibrate the half-width, the other half measure how often
    it actually contains the residual; repeated over random group partitions and
    averaged, because a single split of ten buffers is far too noisy to quote.

    Returns the mean held-out coverage, its spread, and the mean half-width the
    calibration halves produced -- so a half-width that only covers because it
    is enormous is visible rather than flattering.
    """
    uniq = np.unique(groups)
    if len(uniq) < 4:
        return {
            "split_conformal_coverage": None,
            "split_conformal_note": (
                f"only {len(uniq)} groups; a group-split coverage estimate would be "
                "noise. Quote the half-width as an error scale and nothing more."
            ),
        }
    rng = np.random.default_rng(seed)
    covs, widths = [], []
    for _ in range(n_splits):
        perm = rng.permutation(uniq)
        cal_groups = set(perm[: len(perm) // 2].tolist())
        cal = np.isin(groups, list(cal_groups))
        if cal.all() or not cal.any():
            continue
        q = float(np.quantile(resid[cal], 1 - alpha, method="higher"))
        covs.append(float((resid[~cal] <= q).mean()))
        widths.append(q)
    if not covs:
        return {"split_conformal_coverage": None,
                "split_conformal_note": "no usable group split"}
    return {
        "split_conformal_coverage": float(np.mean(covs)),
        "split_conformal_coverage_sd": float(np.std(covs)),
        "split_conformal_halfwidth_mean": float(np.mean(widths)),
        "split_conformal_n_splits": len(covs),
        "split_conformal_note": (
            f"half the groups calibrate the half-width, the other half measure "
            f"coverage; mean over {len(covs)} random group partitions at "
            f"nominal {1 - alpha:.0%}"
        ),
    }


def train_head(
    rows,
    spec: TaskSpec,
    *,
    use_conditions: bool = True,
    n_seeds: int = 5,
    n_folds: int = 5,
    cluster_threshold: float = 0.90,
    conformal_alpha: float = 0.1,
    seed: int = 0,
    verbose: bool = True,
) -> Head | None:
    t0 = time.time()
    seqs, conds, X, y, meta = _rows_to_arrays(rows, spec, use_conditions)
    if len(y) < spec.min_rows:
        if verbose:
            print(f"  [skip] {spec.name}: only {len(y)} labelled rows (need {spec.min_rows})")
        return None

    if spec.group_by == "condition":
        groups = _condition_groups(conds)
    else:
        groups = cluster_sequences(seqs, threshold=cluster_threshold)
    n_groups = len(np.unique(groups))
    folds = min(n_folds, max(2, n_groups))

    if spec.task == "multiclass" or spec.task == "binary":
        # StratifiedGroupKFold needs each class present; fall back if not
        try:
            splitter = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
            split = list(splitter.split(X, y, groups))
        except ValueError:
            split = list(GroupKFold(n_splits=folds).split(X, y, groups))
    else:
        split = list(GroupKFold(n_splits=folds).split(X, y, groups))

    n_classes = len(spec.classes) if spec.task == "multiclass" else 2
    if spec.task == "regression":
        oof = np.full(len(y), np.nan)
    else:
        oof = np.full((len(y), n_classes), np.nan)

    for tr, te in split:
        preds = []
        for s in range(n_seeds):
            m = _estimator(spec.task, seed + s, len(tr))
            m.fit(X[tr], y[tr])
            preds.append(m.predict(X[te]) if spec.task == "regression" else m.predict_proba(X[te]))
        p = np.mean(preds, axis=0)
        if spec.task == "regression":
            oof[te] = p
        else:
            if p.shape[1] != n_classes:  # a fold missed a class
                full = np.zeros((len(te), n_classes))
                full[:, : p.shape[1]] = p
                p = full
            oof[te] = p

    # ------------------------------------------------------------- calibration
    calibrator = None
    metrics: dict = {}
    if spec.task == "binary":
        from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

        p_raw = oof[:, 1]
        metrics["auroc"] = float(roc_auc_score(y, p_raw))
        metrics["auprc"] = float(average_precision_score(y, p_raw))
        metrics["brier_uncalibrated"] = float(brier_score_loss(y, p_raw))
        metrics["ece_uncalibrated"] = expected_calibration_error(y, p_raw)
        # Choose the calibrator by cross-fitted ECE rather than assuming one.
        # Isotonic is the more flexible fit but on small datasets it produces
        # wide flat plateaus, so two very different sequences can come back with
        # the same probability. Platt is smooth and monotone and often wins here.
        # A small margin favours Platt to avoid trading resolution for noise.
        cal_scores = {}
        for method in ("isotonic", "platt"):
            pc = _crossfit_calibrated(oof, y, groups, "binary", 2, method=method)[:, 1]
            cal_scores[method] = (expected_calibration_error(y, pc), pc)
        best = "platt" if cal_scores["platt"][0] <= cal_scores["isotonic"][0] + 0.005 else "isotonic"
        metrics["calibration_method"] = best
        metrics["ece_isotonic"] = cal_scores["isotonic"][0]
        metrics["ece_platt"] = cal_scores["platt"][0]
        calibrator = BinaryCalibrator(best).fit(p_raw, y)
        p_cal = cal_scores[best][1]
        metrics["n_distinct_calibrated_values"] = int(len(np.unique(np.round(p_cal, 4))))
        metrics["brier"] = float(brier_score_loss(y, p_cal))
        metrics["ece"] = expected_calibration_error(y, p_cal)
        metrics["mce"] = maximum_calibration_error(y, p_cal)
        metrics["accuracy_at_0.5"] = float(((p_cal >= 0.5).astype(int) == y).mean())
        conf, acc, cnt = reliability_curve(y, p_cal)
        metrics["reliability"] = {"confidence": conf.tolist(), "accuracy": acc.tolist(),
                                  "count": cnt.tolist()}
        metrics["positive_rate"] = float(np.mean(y))
        metrics["calibration_note"] = "calibrated metrics are cross-fitted across sequence groups"

    elif spec.task == "multiclass":
        from sklearn.metrics import balanced_accuracy_score, roc_auc_score

        metrics["balanced_accuracy"] = float(balanced_accuracy_score(y, oof.argmax(1)))
        metrics["accuracy"] = float((oof.argmax(1) == y).mean())
        try:
            metrics["auroc_ovr_macro"] = float(
                roc_auc_score(y, oof, multi_class="ovr", average="macro")
            )
        except ValueError:
            metrics["auroc_ovr_macro"] = float("nan")
        metrics["brier_uncalibrated"] = brier_multiclass(y, oof)
        metrics["ece_uncalibrated"] = multiclass_ece(y, oof)
        calibrator = TemperatureScaler().fit(oof, y)
        cal = _crossfit_calibrated(oof, y, groups, "multiclass", n_classes)
        metrics["temperature"] = calibrator.temperature
        metrics["brier"] = brier_multiclass(y, cal)
        metrics["ece"] = multiclass_ece(y, cal)
        correct = (cal.argmax(1) == y).astype(float)
        conf, acc, cnt = reliability_curve(correct, cal.max(1))
        metrics["reliability"] = {"confidence": conf.tolist(), "accuracy": acc.tolist(),
                                  "count": cnt.tolist()}
        per_class = {}
        for i, c in enumerate(spec.classes):
            m = y == i
            per_class[c] = {"n": int(m.sum()),
                            "recall": float((cal.argmax(1)[m] == i).mean()) if m.any() else None}
        metrics["per_class"] = per_class
        metrics["calibration_note"] = "calibrated metrics are cross-fitted across sequence groups"

    else:
        from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

        metrics["r2"] = float(r2_score(y, oof))
        metrics["rmse"] = float(np.sqrt(mean_squared_error(y, oof)))
        metrics["mae"] = float(mean_absolute_error(y, oof))
        metrics["spearman"] = float(
            __import__("scipy.stats", fromlist=["spearmanr"]).spearmanr(y, oof).statistic
        )
        resid = np.abs(y - oof)
        q = float(np.quantile(resid, 1 - conformal_alpha, method="higher"))
        metrics["oof_residual_halfwidth"] = q
        metrics["oof_residual_alpha"] = conformal_alpha
        metrics["in_sample_coverage_of_residual_pool"] = float((resid <= q).mean())

        # Grouped split-conformal, which the number above is not.
        #
        # Taking the quantile of a residual pool and then reporting the fraction
        # of that same pool it covers is circular: the coverage is ~1-alpha by
        # construction, whatever the model does. A real coverage estimate needs
        # the quantile fitted on one set of groups and measured on another, and
        # the split has to be by group rather than by row -- splitting rows would
        # put near-duplicate sequences on both sides and inflate coverage the
        # same way a random CV split inflates R^2.
        split = _grouped_conformal(resid, groups, alpha=conformal_alpha, seed=seed)
        metrics.update(split)
        metrics["interval_note"] = (
            "interval half-width = oof_residual_halfwidth; "
            "coverage to quote = split_conformal_coverage (held-out groups), "
            "not in_sample_coverage_of_residual_pool"
        )
        metrics["target_range"] = [float(np.min(y)), float(np.max(y))]

    # ------------------------------------------------------------- final refit
    estimators = []
    for s in range(n_seeds):
        m = _estimator(spec.task, seed + 100 + s, len(y))
        m.fit(X, y)
        estimators.append(m)

    head = Head(
        name=spec.name,
        kind=spec.kind,
        task=spec.task,
        target=spec.target,
        feature_names=feature_names(spec.kind, use_conditions=use_conditions),
        classes=list(spec.classes),
        use_conditions=use_conditions,
        estimators=estimators,
        calibrator=calibrator,
        conformal_q=metrics.get("oof_residual_halfwidth"),
        conformal_alpha=conformal_alpha,
        metrics=metrics,
        training_meta={
            "n_rows": int(len(y)),
            "n_features": int(X.shape[1]),
            "tiers": sorted({m["tier"] for m in meta}),
            "sources": sorted({m["source"] for m in meta}),
            "grouping": group_report(groups),
            "group_by": spec.group_by,
            "cluster_threshold": cluster_threshold,
            "n_folds": folds,
            "n_seeds": n_seeds,
            "seconds": round(time.time() - t0, 1),
        },
        applicability=_applicability(conds, seqs),
    )
    if verbose:
        key = {"binary": "auroc", "multiclass": "balanced_accuracy", "regression": "r2"}[spec.task]
        print(f"  [ok]   {spec.name}: n={len(y)} groups={n_groups} "
              f"{key}={metrics.get(key):.3f} ece={metrics.get('ece', float('nan')):.3f} "
              f"({head.training_meta['seconds']}s)")
    return head


def train_all(
    atlas,
    *,
    tasks: Sequence[TaskSpec] = tuple(DEFAULT_TASKS),
    use_conditions: bool = True,
    allow_predicted: bool = True,
    **kw,
) -> MultiTaskModel:
    model = MultiTaskModel()
    all_rows = []
    for spec in tasks:
        tiers = list(spec.tiers)
        if "predicted" in tiers and not allow_predicted:
            print(f"  [skip] {spec.name}: needs the predicted tier "
                  f"(pass allow_predicted=True to train it)")
            continue
        rows = atlas.query(kind=spec.kind, tiers=tiers, label=spec.target,
                           sources=spec.sources, label_classes=spec.label_classes)
        if spec.max_rows and len(rows) > spec.max_rows:
            rs = np.random.RandomState(0)
            keep = rs.choice(len(rows), spec.max_rows, replace=False)
            rows = [rows[i] for i in sorted(keep)]
        kw2 = dict(kw)
        if spec.n_seeds:
            kw2["n_seeds"] = spec.n_seeds
        all_rows.extend(rows)
        head = train_head(rows, spec, use_conditions=use_conditions, **kw2)
        if head is not None:
            head.training_meta["claim"] = claims.CORRECTED_CLAIMS.get(
                spec.name, spec.claim)
            head.training_meta["label_classes"] = sorted(
                {r["label_class"] for r in rows} if rows else set()
            )
            # Three fields, not one. "experimental" in tiers answers "was
            # anything measured", which is not "was folding measured" -- the
            # collapse of those two is what let v0.4.0's model card call a
            # CUT&Tag proxy head measurement-grounded. quadcond/claims.py
            # derives all three from the head's own sources and tiers.
            head.training_meta.update({
                k: v for k, v in claims.semantics(
                    spec.name, head.training_meta, head.task).items()
                if k != "calibration_scope"
            })
            head.metrics["calibration_scope"] = claims.semantics(
                spec.name, head.training_meta, head.task)["calibration_scope"]
            head.training_meta["mixed_label_classes"] = (
                len(head.training_meta["label_classes"]) > 1
            )
        if head is not None:
            model.add(head)
    model.dataset_fingerprint = rows_fingerprint(all_rows)
    model.atlas_snapshot = {
        "path": str(atlas.path),
        "n_records": atlas.count(),
        "summary": atlas.summary().to_dict(orient="records"),
    }
    return model
