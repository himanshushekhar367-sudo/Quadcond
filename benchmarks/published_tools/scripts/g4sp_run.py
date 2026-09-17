import sys, pickle, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
A = "/home/claude/tools/donn-liew_G4ShapePredictor/g4sp application code/"
sys.argv = ["x"]
def conv(s, pad=100):
    m = {'A':1,'T':2,'C':3,'G':4,'N':0}; l=[m.get(c,0) for c in s.upper()]
    while len(l) < pad: l.insert(0,0); l.append(0)
    return np.array(l[:pad])
d = pd.read_csv("oof_g4_topology_full.csv")
X = np.array([conv(s) for s in d.sequence]); y = d.y.values  # QuadCond order: parallel, antiparallel, hybrid == G4SP 0,1,2
out = pd.DataFrame({"sequence": d.sequence, "y": y, "fold": d.fold})
from sklearn.base import clone
for name in ["LightGBMClassifier", "CatBoostClassifier", "XGBClassifier"]:
    m = pickle.load(open(A + name + ".pkl", "rb"))
    P = m.predict_proba(X)
    for i, c in enumerate(["parallel","antiparallel","hybrid"]): out[f"{name}_pretrained_p_{c}"] = P[:, i]
    # retrain same estimator/hyper-parameters under QuadCond's grouped folds
    oof = np.zeros((len(y), 3))
    for k in sorted(d.fold.unique()):
        tr, te = d.fold.values != k, d.fold.values == k
        if name == "CatBoostClassifier":
            from catboost import CatBoostClassifier
            mm = CatBoostClassifier(**{kk: vv for kk, vv in m.get_params().items()}); mm.set_params(verbose=0)
        elif name == "XGBClassifier":
            break
        else:
            mm = clone(m)
        mm.fit(X[tr], y[tr]); oof[te] = mm.predict_proba(X[te])
    else:
        for i, c in enumerate(["parallel","antiparallel","hybrid"]): out[f"{name}_grouped_p_{c}"] = oof[:, i]
    print(name, flush=True)
out.to_csv("scores_g4sp.csv", index=False)
