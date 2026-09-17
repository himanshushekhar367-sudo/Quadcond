import sys, pickle, importlib.util, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.argv = ["x"]
spec = importlib.util.spec_from_file_location("ims", "/home/claude/tools/YANGB1_iM-Seeker/iM-Seeker.py")
ims = importlib.util.module_from_spec(spec); spec.loader.exec_module(ims)
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
t = pd.read_csv("imseeker/out3/iM-seeker_final_prediction.txt", sep="\t")
parts = t["Putative i-motif id"].str.split("|")
t["id"] = parts.str[0].str.replace(r"[+-]$", "", regex=True)
t["comp"] = parts.apply(lambda p: p[-7:])
rep = t.sort_values("Folding strength", ascending=False).drop_duplicates("id")
a = pd.read_csv("allseqs.csv").set_index("sequence")
d = pd.read_csv("oof_im_pht_full.csv"); d["sequence"] = d.sequence.str.upper(); d["id"] = a.reindex(d.sequence).id.values
d = d.merge(rep[["id", "comp", "Folding strength"]], on="id", how="left")
has = d.comp.notna().values
X = np.array([ims.getFeatures(*c) if isinstance(c, list) else [np.nan] * 1 for c in d.comp], dtype=object)
Xh = np.array([list(x) for x in X[has]], float)
reg = pickle.load(open("/mnt/user-data/uploads/Quadcond/24587160/pickle_model_regression.pkl", "rb"))
import xgboost as xgb
params = reg.get_xgb_params()
oof = np.full(len(d), np.nan)
yh = d.y.values[has]; fh = d.fold.values[has]; pred = np.full(has.sum(), np.nan)
for k in np.unique(fh):
    tr, te = fh != k, fh == k
    m = xgb.XGBRegressor(**{kk: v for kk, v in params.items() if v is not None}, n_estimators=reg.n_estimators)
    m.fit(Xh[tr], yh[tr]); pred[te] = m.predict(Xh[te])
oof[has] = pred
d["imseeker_grouped"] = oof
d[["sequence", "y", "group", "fold", "Folding strength", "imseeker_grouped"]].to_csv("scores_imseeker_grouped.csv", index=False)
print("covered", has.sum(), "of", len(d), "spearman grouped", spearmanr(yh, pred)[0], "r2 (pH units)", r2_score(yh, pred),
      "qc on same rows", spearmanr(yh, d.pred.values[has])[0], r2_score(yh, d.pred.values[has]))
