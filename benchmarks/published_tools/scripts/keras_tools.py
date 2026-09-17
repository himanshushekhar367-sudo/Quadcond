import os, sys, numpy as np, pandas as pd
os.environ["TF_USE_LEGACY_KERAS"] = "1"; os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
import tensorflow as tf
import tf_keras as keras
T = "/home/claude/tools/"
seqs = pd.read_csv("allseqs.csv")
S = list(seqs.sequence)
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
out = pd.DataFrame({"id": seqs.id})
# ---- DeepG4 (classic, sequence-only), one-hot A,C,G,T, right zero-pad to 201, crop at 201
m = keras.models.load_model(T+"raphaelmourad_DeepG4/inst/extdata/DeepG4_classic_rescale_BW_sampling_02_03_2021/2021-03-02T16-17-28Z/best_model.h5", compile=False)
print("DeepG4 input", m.input_shape)
L = m.input_shape[1]
def oh_deepg4(s):
    a = np.zeros((L,4), np.float32); idx={"A":0,"C":1,"G":2,"T":3}
    for i,c in enumerate(s[:L]):
        if c in idx: a[i,idx[c]] = 1
    return a
p1 = m.predict(np.stack([oh_deepg4(s) for s in S]), batch_size=512, verbose=0).ravel()
p2 = m.predict(np.stack([oh_deepg4(rc(s)) for s in S]), batch_size=512, verbose=0).ravel()
out["deepg4"] = np.maximum(p1, p2); print("deepg4 done", flush=True)
# ---- G4detector (K+, sequence-only base models; 124-nt one-hot ACGT centred, N padding = zeros)
def oh_g4det(s, win=124):
    if len(s) > win:
        c = round(len(s)/2); s = s[c-62:c+62]
    elif len(s) < win:
        z1 = round((win-len(s))/2); s = "N"*z1 + s + "N"*(win-len(s)-z1)
    a = np.zeros((win,4), np.float32); idx={"A":0,"C":1,"G":2,"T":3}
    for i,c in enumerate(s):
        if c in idx: a[i,idx[c]] = 1
    return a
for name in ["K_pq_125nt_base_model", "K_random_125nt_base_model", "K_dishuffle_125nt_base_model"]:
    sm = tf.saved_model.load(T+"OrensteinLab_G4detector/models/"+name)
    fn = sm.signatures["serving_default"]
    kw = list(fn.structured_input_signature[1].keys())[0]
    lay = lambda x: fn(**{kw: tf.constant(x, dtype=fn.structured_input_signature[1][kw].dtype)})
    X = np.stack([oh_g4det(s) for s in S])
    preds = []
    for i in range(0, len(X), 512):
        r = lay(X[i:i+512]); r = list(r.values())[0] if isinstance(r, dict) else r
        preds.append(np.asarray(r).ravel())
    out["g4detector_"+name.split("_125")[0]] = np.concatenate(preds); print(name, flush=True)
# ---- G4mismatch whole-genome K model: 15-nt windows with 100-nt flanks, N=0.25; max over windows, both strands
sys.path.insert(0, T+"OrensteinLab_G4mismatch")
from models import pearson
mm = keras.models.load_model(T+"OrensteinLab_G4mismatch/models/model_K.h5", custom_objects={"pearson": pearson}, compile=False)
print("G4mismatch input", mm.input_shape)
W = mm.input_shape[1]; flank = (W-15)//2
def oh_mm(s):
    tr = {"A":0,"C":1,"G":2,"T":3}
    a = np.full((len(s),4), 0.25, np.float32)
    for i,c in enumerate(s):
        if c in tr: a[i] = 0; a[i,tr[c]] = 1
    return a
out.to_csv("scores_keras_partial.csv", index=False)
res = np.full(len(S), -np.inf)
buf, owner = [], []
def flush():
    global buf, owner
    if not buf: return
    pr = mm.predict(np.stack(buf), batch_size=1024, verbose=0).ravel()
    np.maximum.at(res, np.array(owner), pr)
    buf, owner = [], []
for k, s in enumerate(S):
    for q in (s, rc(s)):
        c = oh_mm("N"*flank + q + "N"*max(flank, W - flank - len(q)))
        for i in range(0, max(1, len(c) - W + 1), 2):
            buf.append(c[i:i+W]); owner.append(k)
    if len(buf) > 20000: flush()
    if k % 2000 == 0: print("mm", k, flush=True)
flush()
out["g4mismatch"] = res
out.to_csv("scores_keras.csv", index=False)
