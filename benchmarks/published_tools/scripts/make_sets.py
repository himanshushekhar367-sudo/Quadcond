import pandas as pd, numpy as np, random
from rules import *
# unique sequences across all sets -> one FASTA for external tools
seqs = set()
for f in ["oof_g4_fold_full.csv","oof_g4_topology_full.csv","oof_g4_tm_full.csv","oof_im_fold_full.csv","oof_im_pht_full.csv"]:
    seqs |= set(pd.read_csv(f).sequence.str.upper())
# external G4-seq (G4detector benchmark, human, K+): positives vs PQ-regex negatives and random negatives
import zipfile, io
def read_fa_zip(z, name, n, seed):
    with zipfile.ZipFile(z) as zz:
        txt = zz.read(name).decode()
    recs = [r.split("\n",1) for r in txt.split(">")[1:]]
    recs = [(h.strip(), s.replace("\n","").upper()) for h,s in recs]
    random.Random(seed).shuffle(recs); return recs[:n]
B="/home/claude/tools/OrensteinLab_G4detector/benchmarks/"
pos = read_fa_zip(B+"K_positives.zip","pos_ex_K_125.fa",3000,1)
npq = read_fa_zip(B+"human_K_negatives.zip","neg_ex_K_pq_125.fa",3000,2)
nrd = read_fa_zip(B+"human_K_negatives.zip","neg_ex_K_random_125.fa",3000,3)
rows=[(s,1,"pos") for _,s in pos]+[(s,0,"pq") for _,s in npq]+[(s,0,"random") for _,s in nrd]
ext=pd.DataFrame(rows,columns=["sequence","y","negtype"]); ext=ext[~ext.sequence.str.contains("[^ACGT]")]
ext.to_csv("ext_g4seq.csv",index=False); print(ext.negtype.value_counts(), ext.sequence.str.len().describe())
seqs |= set(ext.sequence)
seqs = sorted(seqs)
ids = {s:f"s{i}" for i,s in enumerate(seqs)}
pd.DataFrame({"id":[ids[s] for s in seqs],"sequence":seqs}).to_csv("allseqs.csv",index=False)
with open("allseqs.fa","w") as fh:
    for s in seqs: fh.write(f">{ids[s]}\n{s}\n")
print(len(seqs))
# rule scores
R=[]
for s in seqs:
    m,mx,mn = g4hunter(s)
    R.append(dict(id=ids[s],sequence=s,g4h_mean=m,g4h_max=mx,g4h_min=mn,
        qgrs=qgrs_gscore(s) if len(s)<=130 else np.nan,
        g4_regex=regex_hits(s,G4_CANON), g4_regex_long=regex_hits(s,G4_LONG), g2_regex=regex_hits(s,G2_REL),
        im_regex=regex_hits(s,IM_CANON), gc=gc(s), length=len(s)))
pd.DataFrame(R).to_csv("scores_rules.csv",index=False)
