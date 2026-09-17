import sys, pandas as pd, numpy as np
sys.path.insert(0,"/home/claude/qc")
from quadcond.models.predict import Predictor
from quadcond.conditions import Condition
from rules import revcomp
p = Predictor.load("/home/claude/qc/artifacts/quadcond_model.joblib")
cond = Condition(k=100)
ext = pd.read_csv("ext_g4seq.csv")
out = []
B = 200
seqs = list(ext.sequence)
for i in range(0, len(seqs), B):
    chunk = seqs[i:i+B]
    fw = p.predict(chunk, cond, heads=["g4_fold","g4_fold_genomic"], n_neighbours=0)
    rv = p.predict([revcomp(s) for s in chunk], cond, heads=["g4_fold","g4_fold_genomic"], n_neighbours=0)
    for s, a, b in zip(chunk, fw, rv):
        el = p.scan(s, cond, kinds=("G4",), both_strands=True)
        sc = [e["predictions"]["g4_fold"].get("probability") for e in el]
        sc = [x for x in sc if x is not None]
        out.append(dict(sequence=s,
            qc_whole=max(a["predictions"]["g4_fold"]["probability"], b["predictions"]["g4_fold"]["probability"]),
            qc_genomic=max(a["predictions"]["g4_fold_genomic"]["probability"], b["predictions"]["g4_fold_genomic"]["probability"]),
            qc_scan=max(sc) if sc else 0.0, qc_n_elements=len(el)))
    print(i, flush=True)
pd.DataFrame(out).to_csv("scores_qc_ext.csv", index=False)
