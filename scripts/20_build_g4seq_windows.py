"""Assemble the G4-seq (K+) window dataset with a chromosome hold-out and a mouse external set."""
import zipfile, random, re, pandas as pd
B = "/home/claude/tools/OrensteinLab_G4detector/benchmarks/"
def read(z, name):
    txt = zipfile.ZipFile(z).read(name).decode()
    for rec in txt.split(">")[1:]:
        h, s = rec.split("\n", 1)
        s = s.replace("\n", "").upper()
        if len(s) >= 100 and not re.search("[^ACGT]", s):
            yield h.strip(), s
HOLD = {"chr2", "chr8", "chr17"}
rng = random.Random(7)
rows = []
for z, n, lab, neg in [("K_positives.zip", "pos_ex_K_125.fa", 1, "pos"),
                       ("human_K_negatives.zip", "neg_ex_K_pq_125.fa", 0, "pq"),
                       ("human_K_negatives.zip", "neg_ex_K_random_125.fa", 0, "random")]:
    recs = list(read(B + z, n)); rng.shuffle(recs)
    tr = [r for r in recs if r[0].split(":")[0] not in HOLD]
    te = [r for r in recs if r[0].split(":")[0] in HOLD]
    ntr = {"pos": 60000, "pq": 30000, "random": 30000}[neg]
    nte = {"pos": 6000, "pq": 3000, "random": 3000}[neg]
    rows += [dict(id=h, sequence=s, y=lab, negtype=neg, split="train") for h, s in tr[:ntr]]
    rows += [dict(id=h, sequence=s, y=lab, negtype=neg, split="test_heldout_chr") for h, s in te[:nte]]
    print(n, len(tr), len(te))
for z, n, lab, neg in [("K_positives.zip", "pos_ex_K_Mouse_125.fa", 1, "pos"),
                       ("mouse_K_negatives.zip", "neg_ex_K_Mouse_pq_125.fa", 0, "pq"),
                       ("mouse_K_negatives.zip", "neg_ex_K_Mouse_rand_125.fa", 0, "random")]:
    recs = list(read(B + z, n)); rng.shuffle(recs)
    k = {"pos": 6000, "pq": 3000, "random": 3000}[neg]
    rows += [dict(id="mm:" + h, sequence=s, y=lab, negtype=neg, split="test_mouse") for h, s in recs[:k]]
d = pd.DataFrame(rows).drop_duplicates("sequence")
d.to_csv("g4seq_windows.csv.gz", index=False)
print(d.groupby(["split", "negtype"]).size())
