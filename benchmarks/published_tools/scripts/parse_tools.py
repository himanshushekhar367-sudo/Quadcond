import pandas as pd, numpy as np, re
a = pd.read_csv("allseqs.csv")
# G4Catchall: max |score| per id (0 if no hit)
for tag in ["default", "g2"]:
    try:
        b = pd.read_csv(f"g4catchall/{tag}.bed", sep="\t", header=None)
        m = b.groupby(0)[7].apply(lambda v: np.abs(v).max())
        a[f"g4catchall_{tag}"] = a.id.map(m).fillna(0.0)
    except Exception as e: print("catchall", tag, e)
# G4Boost: max g4_prob and min mfe per id
try:
    g = pd.read_csv("g4boost/allseqs.fa.g4scores.csv", sep="\t")
    a["g4boost_prob"] = a.id.map(g.groupby("seq").g4_prob.max()).fillna(0.0)
    a["g4boost_mfe"] = a.id.map(g.groupby("seq").mfe_pred.min()).fillna(0.0)
except Exception as e: print("g4boost", e)
try:
    g = pd.read_csv("g4boost_g3/allseqs.fa.g4scores.csv", sep="\t")
    a["g4boost_g3_prob"] = a.id.map(g.groupby("seq").g4_prob.max()).fillna(0.0)
    a["g4boost_g3_mfe"] = a.id.map(g.groupby("seq").mfe_pred.min()).fillna(0.0)
except Exception as e: print("g4boost g3", e)
# iM-Seeker: max probability and max strength per id
try:
    t = pd.read_csv("imseeker/out3/iM-seeker_final_prediction.txt", sep="\t")
    t["id"] = t["Putative i-motif id"].str.extract(r"^(s\d+)")
    t["prob"] = t["Folding probability"].str.extract(r",\s*([0-9.eE-]+)\]").astype(float)
    a["imseeker_prob"] = a.id.map(t.groupby("id").prob.max())
    a["imseeker_strength"] = a.id.map(t.groupby("id")["Folding strength"].max())
    a["imseeker_found"] = a.id.isin(t.id)
except Exception as e: print("imseeker", e)
for f in ["scores_rules.csv", "scores_pqsfinder.csv", "scores_keras.csv"]:
    try:
        x = pd.read_csv(f); x = x.drop(columns=[c for c in ["sequence"] if c in x.columns])
        a = a.merge(x, on="id", how="left")
    except Exception as e: print(f, e)
a.to_csv("tool_scores_all.csv", index=False)
print(a.columns.tolist(), len(a))
