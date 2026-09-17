import sys, json, pandas as pd, numpy as np
sys.path.insert(0, "/home/claude/qc")
from quadcond.models.predict import Predictor
from quadcond.conditions import Condition
from quadcond.scans import mutation_scan
p = Predictor.load("/home/claude/qc/artifacts/quadcond_model.joblib")
cond = Condition(k=140, na=10, mg=1, ph=7.4, temperature=37)
ctrls = {"MYC Pu27": "TGGGGAGGGTGGGGAGGGTGGGGAAGG", "KIT c-kit1": "AGGGAGGGCGCTGGGAGGAGGG",
         "KIT c-kit2": "CGGGCGGGCGCGAGGGAGGGG", "VEGFA Pu22": "GGGGCGGGCCGGGGGCGGGGTCC",
         "BCL2 Pu39": "AGGGGCGGGCGCGGGAGGAAGGGGGCGGGAGCGGGGCTG", "hTERT hT21": "GGGGAGGGGCTGGGAGGGCCCG",
         "chr22 negative-control window": "AGGAGGGCAGAGAGCTGGGGCCTCGGACTCACCCGACGCTTGTGATGAGCTGCACCCAGGA"}
rows = []
for name, s in ctrls.items():
    r = mutation_scan(p, s, cond, heads=["g4_tm"])
    txt = json.dumps(r, default=str)
    muts = r["substitutions"]
    deltas = [m["heads"]["g4_tm"].get("delta") for m in muts]
    deltas = [x for x in deltas if x is not None]
    lost = sum(m["heads"]["g4_tm"].get("motif", {}).get("state") == "motif_lost" for m in muts)
    rows.append(dict(control=name, length=len(s), n_substitutions=len(muts), n_with_delta=len(deltas),
                     max_abs_delta_Tm=round(max(map(abs, deltas)), 2) if deltas else None,
                     n_abs_delta_ge_5C=sum(abs(x) >= 5 for x in deltas), n_motif_lost=lost, most_destabilising=min(deltas) if deltas else None))
    if name == "MYC Pu27": open("posctrl_myc_sample.json", "w").write(txt[:4000])
print(pd.DataFrame(rows).to_string(index=False))
pd.DataFrame(rows).to_csv("posctrl_summary.csv", index=False)
