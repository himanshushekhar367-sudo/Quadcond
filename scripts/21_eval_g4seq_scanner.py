import pandas as pd, numpy as np, json, sys
from sklearn.metrics import roc_auc_score, average_precision_score
d = pd.read_csv("test_scores.csv.gz")
d["id"] = ["w%d" % i for i in range(len(d))]
for f, cols in [("test_pqsfinder.csv", None), ("test_keras.csv", None)]:
    try:
        x = pd.read_csv(f); d = d.merge(x, on="id", how="left")
    except Exception as e:
        print("missing", f, e)
TOOLS = [("QuadCond genome-scan (G4-seq)", "score"), ("G4Hunter (max |25-nt window|)", "g4h_absmax"),
         ("G4Hunter (|whole-window mean|)", "g4h_absmean"), ("pqsfinder", "pqsfinder"),
         ("canonical G3 regex count", "g4_regex"), ("DeepG4", "deepg4"), ("G4mismatch (K)", "g4mismatch"),
         ("G4detector (K, PQ-neg model)", "g4detector_K_pq"), ("G4detector (K, random-neg model)", "g4detector_K_random"),
         ("G4detector (K, dishuffle model)", "g4detector_K_dishuffle")]
rows = []
for split in ["test_heldout_chr", "test_mouse"]:
    s = d[d.split == split]
    for name, col in TOOLS:
        if col not in s or s[col].isna().all():
            continue
        for neg in ["random", "pq"]:
            t = s[s.negtype.isin(["pos", neg])].dropna(subset=[col])
            rows.append(dict(split=split, negatives=neg, tool=name,
                             auroc=round(roc_auc_score(t.y, t[col]), 4),
                             auprc=round(average_precision_score(t.y, t[col]), 4), n=len(t)))
R = pd.DataFrame(rows)
R.to_csv("scanner_benchmark.csv", index=False)
print(R.pivot_table(index="tool", columns=["split", "negatives"], values="auroc").to_string())
