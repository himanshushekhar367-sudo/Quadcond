import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/qc")
from sklearn.model_selection import KFold
from quadcond.atlas import Atlas
from quadcond.models import train as T
atlas = Atlas("/home/claude/qc/data/atlas_core.db")
spec = {s.name: s for s in T.DEFAULT_TASKS}["g4_tm"]
rows = atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target, sources=spec.sources, label_classes=spec.label_classes)
seqs, conds, X, y, meta = T._rows_to_arrays(rows, spec, True)
oof = np.zeros(len(y))
for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
    oof[te] = np.mean([T._estimator("regression", s, len(tr)).fit(X[tr], y[tr]).predict(X[te]) for s in range(3)], axis=0)
pd.DataFrame({"sequence": seqs, "y": y, "pred": oof}).to_csv("oof_g4_tm_random.csv", index=False)
from sklearn.metrics import r2_score, mean_absolute_error
print(r2_score(y, oof), mean_absolute_error(y, oof))
