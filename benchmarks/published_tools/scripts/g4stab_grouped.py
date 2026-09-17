"""Retrain the published G4STAB architecture (cloned from trained_model_00, fresh weights,
published optimiser settings) under QuadCond's grouped folds, and also under random 5-fold."""
import os, sys, numpy as np, pandas as pd
os.chdir("/home/claude/tools/donn-liew_G4STAB"); sys.path.insert(0, ".")
sys.argv = ["x"]
import tensorflow as tf
import g4stab_predictor as G
from sklearn.model_selection import KFold, GroupShuffleSplit
d = pd.read_csv("/home/claude/bench/oof_g4_tm_full.csv")
seqs = list(d.sequence.str.upper())
conc = d[["k", "na", "li_nh4"]].fillna(0).values.tolist()
ph = d.ph.fillna(7.0).tolist()
Xo, Xk, Xs, Xp = G.seqmap(seqs, conc, ph)
base = tf.keras.models.load_model("trained_models/trained_model_00.keras",
                                  custom_objects={"CastFloat32": G.CastFloat32}, compile=False)
names = [i.name for i in base.inputs]; shapes = [tuple(i.shape[1:]) for i in base.inputs]
print(names, shapes)
feats = {"ohe": Xo, "kmer": Xk, "salt": Xs, "ph": Xp}
def assign():
    out = {}
    pool = dict(feats)
    for n, sh in zip(names, shapes):
        for k, v in list(pool.items()):
            if tuple(v.shape[1:]) == sh:
                out[n] = v; pool.pop(k); break
    return out
X = assign(); assert len(X) == len(names), X.keys()
y = d.y.values.astype("float32")
def fit_predict(tr, te, seed):
    tf.keras.utils.set_random_seed(seed)
    m = tf.keras.models.clone_model(base)
    m.compile(optimizer=tf.keras.optimizers.Nadam(learning_rate=1.55e-4, beta_1=0.995, beta_2=0.9915, weight_decay=4e-5), loss="mse")
    gss = GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=seed)
    a, b = next(gss.split(tr, groups=d.group.values[tr])); trn, val = tr[a], tr[b]
    es = tf.keras.callbacks.EarlyStopping(patience=15, restore_best_weights=True)
    m.fit({k: v[trn] for k, v in X.items()}, y[trn], validation_data=({k: v[val] for k, v in X.items()}, y[val]),
          epochs=250, batch_size=32, verbose=0, callbacks=[es])
    return m.predict({k: v[te] for k, v in X.items()}, verbose=0).ravel()
res = {}
for scheme in ["grouped"]:
    oof = np.zeros(len(y))
    folds = d.fold.values if scheme == "grouped" else np.zeros(len(y), int)
    if scheme == "random":
        for k, (_, te) in enumerate(KFold(5, shuffle=True, random_state=0).split(y)): folds[te] = k
    for k in range(5):
        tr = np.where(folds != k)[0]; te = np.where(folds == k)[0]
        oof[te] = np.mean([fit_predict(tr, te, s) for s in (0,)], axis=0)
        print(scheme, k, flush=True)
    res[f"g4stab_retrained_{scheme}"] = oof
    res[f"fold_{scheme}"] = folds
out = pd.DataFrame(res); out["sequence"] = seqs; out["y"] = y
out.to_csv("/home/claude/bench/scores_g4stab_retrained.csv", index=False)
