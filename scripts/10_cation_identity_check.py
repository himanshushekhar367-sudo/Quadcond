#!/usr/bin/env python3
"""Does g4_tm know K+ from Na+ -- on sequences it was never shown?

Potassium stabilises a G-quadruplex substantially more than sodium at the same
concentration; it is the single best-known fact about G4 buffer chemistry. A head
reporting R^2 0.670 could still be reading total monovalent concentration and
ignoring which cation is present, and the aggregate metric would not notice. The
test: every sequence in the atlas measured *both* K+-only and Na+-only at
comparable concentration and pH, predicted under both, measured delta against
predicted delta.

**Why this file was rewritten.** v0.4.1 ran that comparison with the shipped
model -- which was trained on the very records being compared. It is a
mechanistic consistency check (does the fitted function respond to cation
identity at all?) and it is not validation, because the answer is partly
memorised. v0.4.2 scores every pair with a model that never saw its sequence or
anything within the 0.90 k-mer-similarity cluster around it, by refitting g4_tm
under the same grouped folds the head is trained with and predicting each pair
from the fold that held its cluster out. Both numbers are reported: the in-sample
one is what memorisation buys, the out-of-fold one is what the head knows.

**Confidence intervals resample sequence clusters, not pairs.** The 411 pairs
come from ~151 sequences, several of which contribute a dozen pairs each, and
near-duplicates contribute the same information twice. A pair-level bootstrap
would treat those as independent draws and report an interval several times too
tight.

**Monotonicity.** Tm rises with [K+] over the range these measurements cover. A
gradient-boosted ensemble has no such constraint, and on sparse coverage it can
produce a titration curve that goes back down. The ladder walks the *shipped*
model, deliberately: it is a property of the artefact users will run, not a
validation statistic.

    python scripts/10_cation_identity_check.py
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from quadcond.atlas import Atlas
from quadcond.conditions import Condition
from quadcond.features import featurize
from quadcond.models.grouping import cluster_sequences
from quadcond.models.predict import Predictor
from quadcond.models.train import DEFAULT_TASKS, _estimator, _rows_to_arrays

# Constructs whose K+ response is documented well enough to eyeball.
LADDER_SEQUENCES = {
    "hTelo 22AG": "AGGGTTAGGGTTAGGGTTAGGG",
    "c-MYC Pu22": "TGAGGGTGGGTAGGGTGGGTAA",
    "c-KIT1": "AGGGAGGGCGCTGGGAGGAGGG",
    "BCL2 Pu39": "AGGGGCGGGCGCGGGAGGAAGGGGGCGGGAGCGGGGCTG",
}
LADDER_K = (0.0, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 150.0, 200.0)


def matched_pairs(db: str, *, conc_tol_frac: float = 0.15, ph_tol: float = 0.3):
    """Sequences measured under K+-only and Na+-only at comparable conditions."""
    con = sqlite3.connect(f"file:{Path(db)}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "SELECT sequence, k, na, ph, tm FROM records "
            "WHERE evidence_tier='experimental' AND label_class='biophysical' "
            "AND tm IS NOT NULL AND li_nh4=0 AND mg=0"
        ).fetchall()
    finally:
        con.close()

    k_only: dict[str, list] = {}
    na_only: dict[str, list] = {}
    for seq, k, na, ph, tm in rows:
        if na == 0 and k > 0:
            k_only.setdefault(seq, []).append((k, ph, tm))
        elif k == 0 and na > 0:
            na_only.setdefault(seq, []).append((na, ph, tm))

    pairs = []
    for seq in set(k_only) & set(na_only):
        for k, phk, tmk in k_only[seq]:
            for na, phn, tmn in na_only[seq]:
                if abs(k - na) <= max(5.0, conc_tol_frac * max(k, na)) \
                        and abs(phk - phn) <= ph_tol:
                    pairs.append({"sequence": seq, "k": k, "na": na,
                                  "ph": (phk + phn) / 2,
                                  "tm_k": tmk, "tm_na": tmn})
    return pairs



# --------------------------------------------------------------- out-of-fold
class OutOfFold:
    """g4_tm refit under the head's own grouped folds, queryable by sequence.

    ``predict(seq, cond)`` returns the prediction from the fold in which ``seq``
    was held out, so no pair is ever scored by a model that saw it. Sequences the
    training set does not contain have no held-out fold and are reported as
    uncovered rather than silently scored with the full model.
    """

    def __init__(self, atlas_path: str, *, n_folds: int = 5, n_seeds: int = 3,
                 cluster_threshold: float = 0.90, seed: int = 0, verbose: bool = True):
        from sklearn.model_selection import GroupKFold

        spec = next(s for s in DEFAULT_TASKS if s.name == "g4_tm")
        atlas = Atlas(atlas_path)
        rows = atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target,
                           sources=spec.sources, label_classes=spec.label_classes)
        seqs, conds, X, y, _ = _rows_to_arrays(rows, spec, True)
        groups = cluster_sequences(seqs, threshold=cluster_threshold)

        # Every row of one sequence carries one group id -- identical sequences
        # have identical k-mer profiles, so the clusterer cannot split them.
        self.seq_group: dict[str, int] = {}
        for s, g in zip(seqs, groups):
            self.seq_group.setdefault(s, int(g))

        folds = min(n_folds, max(2, len(np.unique(groups))))
        self.models: list[list] = []
        self.group_fold: dict[int, int] = {}
        for f, (tr, te) in enumerate(GroupKFold(n_splits=folds).split(X, y, groups)):
            for g in np.unique(groups[te]):
                self.group_fold[int(g)] = f
            fitted = []
            for s in range(n_seeds):
                m = _estimator(spec.task, seed + s, len(tr))
                m.fit(X[tr], y[tr])
                fitted.append(m)
            self.models.append(fitted)
            if verbose:
                print(f"    fold {f + 1}/{folds}: fit on {len(tr)} rows, "
                      f"{len(np.unique(groups[te]))} clusters held out", flush=True)
        self.n_rows, self.n_groups, self.n_folds = len(y), int(len(np.unique(groups))), folds

    def covers(self, seq: str) -> bool:
        return seq in self.seq_group

    def predict(self, seq: str, cond: Condition) -> float:
        f = self.group_fold[self.seq_group[seq]]
        X = featurize([seq], [cond], kind="G4", use_conditions=True)
        return float(np.mean([m.predict(X)[0] for m in self.models[f]]))


def _cond(k: float, na: float, ph: float) -> Condition:
    return Condition.from_mapping(
        {"k": k, "na": na, "ph": ph, "temperature": 25.0}, track_imputed=False)


def _delta_stats(meas: np.ndarray, pred: np.ndarray) -> dict:
    return {
        "measured_delta_mean": round(float(meas.mean()), 2),
        "measured_delta_median": round(float(np.median(meas)), 2),
        "predicted_delta_mean": round(float(pred.mean()), 2),
        "predicted_delta_median": round(float(np.median(pred)), 2),
        "sign_agreement": round(float((np.sign(meas) == np.sign(pred)).mean()), 4),
        "pearson_r": round(float(np.corrcoef(meas, pred)[0, 1]), 4),
        "mae_of_delta": round(float(np.abs(meas - pred).mean()), 2),
    }


def cluster_bootstrap(pairs: list[dict], key: str, cluster_of: dict[str, int],
                      *, n: int = 2000, seed: int = 0) -> dict:
    """Percentile CIs from resampling sequence-similarity clusters.

    The unit of resampling is the cluster, because that is the unit the folds
    hold out and the unit at which these measurements are actually independent:
    hTelo 22AG with one extra flanking base is not a second observation of the
    K+/Na+ effect.
    """
    by: dict[int, list[dict]] = {}
    for p in pairs:
        by.setdefault(cluster_of[p["sequence"]], []).append(p)
    keys = list(by)
    rs = np.random.RandomState(seed)
    draws: dict[str, list[float]] = {}
    for _ in range(n):
        sample = [q for c in rs.choice(len(keys), len(keys), replace=True)
                  for q in by[keys[c]]]
        m = np.array([q["measured_delta"] for q in sample])
        d = np.array([q[key] for q in sample])
        if len(np.unique(m)) < 2 or len(np.unique(d)) < 2:
            continue
        for stat, val in _delta_stats(m, d).items():
            draws.setdefault(stat, []).append(val)
    out = {}
    for stat, vals in draws.items():
        a = np.asarray(vals, dtype=float)
        out[stat] = {"lo": round(float(np.percentile(a, 2.5)), 3),
                     "hi": round(float(np.percentile(a, 97.5)), 3)}
    out["n_clusters"] = len(keys)
    out["n_resamples"] = n
    return out


def cation_identity(pred: Predictor, pairs: list[dict], oof: "OutOfFold | None") -> dict:
    """Predicted vs measured dTm(K - Na), in-sample and out-of-fold.

    The in-sample column exists to be compared against, not to be quoted. If the
    two are close, the head has learned cation identity as a transferable rule;
    if the out-of-fold column collapses, it memorised these constructs.
    """
    meas, in_samp, out_fold, uncovered = [], [], [], []
    for p in pairs:
        ck, cn = _cond(p["k"], 0.0, p["ph"]), _cond(0.0, p["na"], p["ph"])
        vk = pred.predict(p["sequence"], ck, heads=["g4_tm"],
                          n_neighbours=0)[0]["predictions"]["g4_tm"]["value"]
        vn = pred.predict(p["sequence"], cn, heads=["g4_tm"],
                          n_neighbours=0)[0]["predictions"]["g4_tm"]["value"]
        p["in_sample_delta"] = round(vk - vn, 2)
        p["measured_delta"] = round(p["tm_k"] - p["tm_na"], 2)
        if oof is not None and oof.covers(p["sequence"]):
            p["predicted_delta"] = round(
                oof.predict(p["sequence"], ck) - oof.predict(p["sequence"], cn), 2)
            out_fold.append(p["predicted_delta"])
        else:
            uncovered.append(p["sequence"])
            p["predicted_delta"] = None
        meas.append(p["measured_delta"])
        in_samp.append(p["in_sample_delta"])

    meas = np.asarray(meas)
    covered = [p for p in pairs if p["predicted_delta"] is not None]
    res = {
        "n_pairs": len(pairs),
        "n_sequences": len({p["sequence"] for p in pairs}),
        "n_pairs_out_of_fold": len(covered),
        "n_sequences_not_in_training_rows": len(set(uncovered)),
        "in_sample": _delta_stats(meas, np.asarray(in_samp)),
        "worst_misses": _worst(pairs),
        "sequences_with_contradictory_measurements": _contradictory(pairs),
    }
    res["in_sample"]["note"] = (
        "scored with the shipped model, which was trained on these records. A "
        "mechanistic consistency check, not validation.")
    if covered:
        mo = np.array([p["measured_delta"] for p in covered])
        po = np.array([p["predicted_delta"] for p in covered])
        res["out_of_fold"] = _delta_stats(mo, po)
        res["out_of_fold"]["note"] = (
            f"each pair scored by the fold that held its sequence-similarity "
            f"cluster out; {oof.n_folds} grouped folds over {oof.n_groups} "
            f"clusters of {oof.n_rows} rows.")
        res["out_of_fold_ci95"] = cluster_bootstrap(
            covered, "predicted_delta", oof.seq_group)
        res["in_sample_ci95"] = cluster_bootstrap(
            covered, "in_sample_delta", oof.seq_group)
    return res


def _worst(pairs: list[dict], n: int = 8) -> list[dict]:
    """Largest disagreements, one row per sequence.

    Deduplicated deliberately: a sequence with several irreconcilable literature
    values fills the whole list otherwise, and the interesting question is how
    many *distinct* constructs the head gets badly wrong.
    """
    best: dict[str, dict] = {}
    for p in pairs:
        if p.get("predicted_delta") is None:
            continue
        err = abs(p["measured_delta"] - p["predicted_delta"])
        cur = best.get(p["sequence"])
        if cur is None or err > cur["_err"]:
            best[p["sequence"]] = {
                "sequence": p["sequence"], "k": p["k"], "ph": p["ph"],
                "measured_delta": p["measured_delta"],
                "predicted_delta": p["predicted_delta"], "_err": err,
            }
    out = sorted(best.values(), key=lambda d: -d["_err"])[:n]
    for d in out:
        d.pop("_err")
    return out


def _contradictory(pairs: list[dict], spread: float = 10.0) -> list[dict]:
    """Sequences whose own literature values disagree by more than ``spread``.

    Worth separating from model error. Some constructs carry measurements from
    several papers that differ by tens of degrees at the same nominal buffer --
    ``TCGGCGGGAGGGCGGGGCGGGGCGGA`` has matched pairs implying anything from -23
    to +40 C for the same K/Na comparison. No single prediction can agree with
    all of them, and counting that as a modelling failure would be dishonest in
    the flattering direction: it makes the head look worse than the evidence can
    establish, which is its own kind of misreporting.
    """
    by: dict[str, list[float]] = {}
    for p in pairs:
        by.setdefault(p["sequence"], []).append(p["measured_delta"])
    out = []
    for seq, vals in by.items():
        if len(vals) > 1 and (max(vals) - min(vals)) > spread:
            out.append({"sequence": seq, "n_pairs": len(vals),
                        "measured_delta_min": min(vals),
                        "measured_delta_max": max(vals)})
    return sorted(out, key=lambda d: -(d["measured_delta_max"] - d["measured_delta_min"]))


def monotonicity(pred: Predictor, ph: float = 7.0) -> dict:
    out = {}
    total_rev = 0
    for name, seq in LADDER_SEQUENCES.items():
        vals = []
        for k in LADDER_K:
            v = pred.predict(seq,
                             Condition.from_mapping({"k": k, "ph": ph, "temperature": 25.0}),
                             heads=["g4_tm"], n_neighbours=0)[0]["predictions"]["g4_tm"]["value"]
            vals.append(round(v, 2))
        d = np.diff(vals)
        rev = int((d < -0.05).sum())
        total_rev += rev
        out[name] = {
            "k_mM": list(LADDER_K), "tm": vals,
            "span": round(vals[-1] - vals[0], 2),
            "reversals": rev,
            "largest_reversal": round(float(-d.min()), 2) if (d < 0).any() else 0.0,
        }
    out["_summary"] = {
        "sequences": len(LADDER_SEQUENCES),
        "steps_per_sequence": len(LADDER_K) - 1,
        "total_reversals": total_rev,
        "note": (
            "A reversal is a step where predicted Tm falls as [K+] rises. The "
            "ensemble has no monotonicity constraint, so these are artefacts of "
            "sparse coverage rather than a physical claim. The overall span is "
            "the number that carries the condition response."
        ),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="artifacts/quadcond_model.joblib")
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--out", default="artifacts/cation_identity_check.json")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--in-sample-only", action="store_true",
                    help="skip the refit; reports the v0.4.1 number and nothing else")
    a = ap.parse_args()

    pred = Predictor.load(a.model)
    pairs = matched_pairs(a.db)
    if not pairs:
        sys.exit("no sequence in this atlas was measured under both K+-only and "
                 "Na+-only at comparable concentration; nothing to check")

    oof = None
    if not a.in_sample_only:
        print("refitting g4_tm under grouped folds so no pair is scored by a model "
              "that saw it:", flush=True)
        oof = OutOfFold(a.db, n_folds=a.folds, n_seeds=a.seeds)
        print()

    ident = cation_identity(pred, pairs, oof)
    mono = monotonicity(pred)

    print(f"K+ / Na+ matched pairs: {ident['n_pairs']} over "
          f"{ident['n_sequences']} sequences")
    if oof is not None:
        print(f"  {ident['n_pairs_out_of_fold']} of them have a held-out prediction "
              f"({ident['n_sequences_not_in_training_rows']} sequences absent from "
              "the g4_tm training rows)")
    print()

    def line(tag: str, st: dict, ci: dict | None) -> None:
        def band(stat: str) -> str:
            if not ci or stat not in ci:
                return ""
            return f"  [95% CI {ci[stat]['lo']:+.1f}, {ci[stat]['hi']:+.1f}]"
        print(f"  {tag}")
        print(f"    predicted dTm(K - Na)  mean {st['predicted_delta_mean']:+.1f}"
              f"{band('predicted_delta_mean')}")
        agr = f"{st['sign_agreement']:.0%}"
        if ci and "sign_agreement" in ci:
            agr += f" [{ci['sign_agreement']['lo']:.0%}, {ci['sign_agreement']['hi']:.0%}]"
        r = f"{st['pearson_r']:.3f}"
        if ci and "pearson_r" in ci:
            r += f" [{ci['pearson_r']['lo']:.2f}, {ci['pearson_r']['hi']:.2f}]"
        print(f"    sign agreement {agr}   pearson r {r}   "
              f"MAE {st['mae_of_delta']:.1f} C")

    print(f"  measured   dTm(K - Na)  mean {ident['in_sample']['measured_delta_mean']:+.1f}"
          f"  median {ident['in_sample']['measured_delta_median']:+.1f}\n")
    if "out_of_fold" in ident:
        line("OUT-OF-FOLD  (the number to quote)",
             ident["out_of_fold"], ident.get("out_of_fold_ci95"))
        print()
    line("in-sample  (v0.4.1 reported this; it is a consistency check, not validation)",
         ident["in_sample"], ident.get("in_sample_ci95"))

    if "out_of_fold" in ident:
        o, i = ident["out_of_fold"], ident["in_sample"]
        print(f"\n  memorisation accounts for {i['pearson_r'] - o['pearson_r']:+.3f} of r "
              f"and {i['sign_agreement'] - o['sign_agreement']:+.1%} of sign agreement.")
        ci = ident.get("out_of_fold_ci95", {})
        lo = ci.get("sign_agreement", {}).get("lo")
        if lo is not None and lo > 0.5:
            print("  The held-out interval excludes chance sign agreement: on a sequence "
                  "\n  it has never seen, the head still separates potassium from sodium "
                  "\n  in the right direction.")
        elif lo is not None:
            print("  The held-out interval includes chance sign agreement. Cation identity "
                  "\n  is not established as a transferable rule by this evidence.")

    print("\n  largest disagreements, out-of-fold (one row per sequence):")
    for w in ident["worst_misses"][:5]:
        print(f"    {w['sequence'][:26]:<26} K/Na {w['k']:g} pH {w['ph']:g}   "
              f"measured {w['measured_delta']:+6.1f}   predicted {w['predicted_delta']:+6.1f}")
    contra = ident["sequences_with_contradictory_measurements"]
    if contra:
        print(f"\n  {len(contra)} of the {ident['n_sequences']} sequences carry literature "
              "values that disagree with\n  each other by more than 10 C at the same nominal "
              "buffer. No prediction can\n  match all of them, so they cap what this check "
              "can establish:")
        for c in contra[:3]:
            print(f"    {c['sequence'][:26]:<26} {c['n_pairs']} pairs spanning "
                  f"{c['measured_delta_min']:+.1f} to {c['measured_delta_max']:+.1f} C")

    print(f"\nK+ ladder on the SHIPPED model, pH 7, {LADDER_K[0]:g}-{LADDER_K[-1]:g} mM")
    print("  (a property of the artefact users will run, not a validation statistic):")
    for name in LADDER_SEQUENCES:
        m = mono[name]
        print(f"  {name:<14} span {m['span']:+6.1f} C   reversals "
              f"{m['reversals']}/{len(LADDER_K) - 1}"
              + (f"   largest {m['largest_reversal']:.1f} C" if m["reversals"] else ""))
    print(f"\n  {mono['_summary']['total_reversals']} non-monotonic steps in total. "
          "The ensemble carries no\n  monotonicity constraint; on the thin part of the "
          "salt range that shows.")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(
        {"cation_identity": ident, "monotonicity": mono}, indent=2))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
