import pandas as pd, numpy as np, itertools, collections
def pairs_for(df, keycols, target):
    out=[]
    df = df.copy(); df["sequence"]=df.sequence.str.upper()
    for key, g in df.groupby(keycols, dropna=False):
        # average replicate measurements of the same sequence in the same buffer
        g2 = g.groupby("sequence").agg(y=(target,"mean"), pred=("pred","mean")).reset_index()
        seqs = list(g2.sequence)
        bylen = collections.defaultdict(list)
        for i,s in enumerate(seqs): bylen[len(s)].append(i)
        for L, idx in bylen.items():
            for a,b in itertools.combinations(idx,2):
                s1,s2 = seqs[a],seqs[b]
                diff=[k for k in range(L) if s1[k]!=s2[k]]
                if len(diff)==1:
                    k=diff[0]
                    out.append(dict(ref=s1,alt=s2,pos=k,ref_base=s1[k],alt_base=s2[k],
                        y_ref=g2.y[a],y_alt=g2.y[b],d_true=g2.y[b]-g2.y[a],
                        qc_ref=g2.pred[a],qc_alt=g2.pred[b],d_qc=g2.pred[b]-g2.pred[a],
                        **{c:v for c,v in zip(keycols, key if isinstance(key,tuple) else (key,))}))
    return pd.DataFrame(out)
d=pd.read_csv("oof_g4_tm_full.csv"); d["yv"]=d.y
p=pairs_for(d.rename(columns={"y":"tmv"}), ["k","na","li_nh4","mg","ph","strand_conc"], "tmv")
p.to_csv("pairs_g4_tm.csv",index=False); print("G4 Tm single-substitution pairs:",len(p), "distinct refs", p.ref.nunique())
print(p.d_true.describe())
d=pd.read_csv("oof_im_pht_full.csv")
p=pairs_for(d.rename(columns={"y":"phv"}), ["k"], "phv")
p.to_csv("pairs_im_pht.csv",index=False); print("iM pHT pairs:",len(p)); print(p.d_true.describe())
