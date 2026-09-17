"""Re-derive QuadCond grouped out-of-fold predictions for the biophysical heads,
exactly as quadcond.models.train.train_head does, and save per-row tables."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/qc")
from sklearn.model_selection import StratifiedGroupKFold, GroupKFold
from quadcond.atlas import Atlas
from quadcond.models import train as T
from quadcond.models.grouping import cluster_sequences

atlas = Atlas("/home/claude/qc/data/atlas_core.db")
want = ["g4_fold", "g4_topology", "g4_tm", "im_fold", "im_pht"]
specs = {s.name: s for s in T.DEFAULT_TASKS}
N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 3
for name in want:
    spec = specs[name]
    rows = atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target,
                       sources=spec.sources, label_classes=spec.label_classes)
    for use_c in (True, False):
        seqs, conds, X, y, meta = T._rows_to_arrays(rows, spec, use_c)
        groups = cluster_sequences(seqs, threshold=0.90)
        folds = min(5, max(2, len(np.unique(groups))))
        if spec.task in ("binary", "multiclass"):
            split = list(StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=0).split(X, y, groups))
        else:
            split = list(GroupKFold(n_splits=folds).split(X, y, groups))
        nc = len(spec.classes) if spec.task == "multiclass" else 2
        oof = np.full(len(y), np.nan) if spec.task == "regression" else np.full((len(y), nc), np.nan)
        fold_id = np.full(len(y), -1)
        for k, (tr, te) in enumerate(split):
            fold_id[te] = k
            ps = []
            for s in range(N_SEEDS):
                m = T._estimator(spec.task, s, len(tr)); m.fit(X[tr], y[tr])
                ps.append(m.predict(X[te]) if spec.task == "regression" else m.predict_proba(X[te]))
            p = np.mean(ps, axis=0)
            if spec.task != "regression" and p.shape[1] != nc:
                f = np.zeros((len(te), nc)); f[:, :p.shape[1]] = p; p = f
            oof[te] = p
        df = pd.DataFrame({"sequence": seqs, "y": y, "group": groups, "fold": fold_id,
                           "source": [m_["source"] for m_ in meta]})
        for c in ("k", "na", "li_nh4", "mg", "ph", "temperature", "strand_conc"):
            df[c] = [getattr(cd, c, None) for cd in conds]
        if spec.task == "regression":
            df["pred"] = oof
        elif spec.task == "binary":
            df["pred"] = oof[:, 1]
        else:
            for i, c in enumerate(spec.classes):
                df[f"p_{c}"] = oof[:, i]
        tag = "full" if use_c else "seqonly"
        df.to_csv(f"/home/claude/bench/oof_{name}_{tag}.csv", index=False)
        print(name, tag, len(df), flush=True)
