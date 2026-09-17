# The genomic path: `quadcond genome-scan`

## Why a second model

`g4_fold` is an oligonucleotide model: purified sequences of 5-99 nt against
dinucleotide shuffles of themselves, which is the task its 0.97 AUROC describes.
Benchmarked on human G4-seq windows it does not transfer
(`benchmarks/published_tools`, section 2):

| | vs random genomic windows | vs motif-matching windows G4-seq did not observe |
|---|---|---|
| `g4_fold`, whole 124-nt window | 0.657 | **0.289** |
| `g4_fold`, best motif in the window | 0.687 | 0.476 |
| G4Hunter (max 25-nt window) | 0.966 | 0.685 |

At or below chance on the hard negatives is not a number to hand a user behind
an applicability warning, so as of 0.5.1 the folding classifiers **refuse** a
window longer than their training range and name this path instead. The genomic
question now has its own model.

## What the scanner is

`quadcond/genome_scan.py`, asset `quadcond_g4seq_scanner.joblib` (1.8 MB,
optional), trained by `scripts/19_train_g4seq_scanner.py`.

- **Label**: an observed G4-seq window (polymerase stalling on purified human
  genomic DNA in K+ buffer; Marsico et al. 2019, GEO GSE110582), taken from the
  G4detector benchmark release. Negatives: random genomic windows, and
  canonical-motif windows G4-seq did not observe.
- **Training**: 119,865 windows of 124 nt (59,988 positive). Chromosomes 2, 8
  and 17 are held out entirely; a mouse G4-seq set is kept as a cross-species
  test. Calibration is isotonic, fitted on chromosome-grouped out-of-fold
  predictions.
- **Features**: the package's own sequence features computed on both strands
  (G-richer strand first, so the score is strand-symmetric) plus ten window
  statistics (windowed G4Hunter, motif counts, G-run structure): 266 in total.
- **Claim**: P(window resembles an observed K+ G4-seq window), calibrated at the
  50 % prevalence of the training set. It is a **genomic in-vitro observation**,
  not a folding measurement in a buffer you choose, so the condition vector is
  not an input. `target_semantics: genomic_in_vitro`,
  `biophysically_grounded: false`.

## Results

AUROC. Bold = best per column. *in-sample* marks a tool trained on these very
files (G4detector's own split is random, so held-out chromosomes are not held
out for it).

| | held-out chr2/8/17: vs random | vs unobserved PQS | mouse: vs random | mouse: vs unobserved PQS |
|---|---|---|---|---|
| **QuadCond genome-scan** | 0.942 | **0.928** | 0.928 | **0.775** |
| G4Hunter (max 25-nt window) | 0.947 | 0.672 | 0.930 | 0.554 |
| G4mismatch (K) | 0.965 | 0.849 | **0.942** | 0.636 |
| G4detector (K, random-neg) *in-sample* | **0.970** | 0.818 | 0.956 | 0.623 |
| G4detector (K, PQ-neg) *in-sample* | 0.441 | 0.928 | 0.470 | 0.761 |
| G4detector (K, dishuffle) *in-sample* | 0.855 | 0.768 | 0.829 | 0.589 |
| DeepG4 | 0.725 | 0.541 | 0.754 | 0.495 |
| canonical G3 regex count | 0.691 | 0.532 | 0.689 | 0.473 |
| `g4_fold` (oligo head, for reference) | 0.657 | 0.289 | — | — |

Internal chromosome-grouped CV: AUROC 0.939. Held-out chromosomes: Brier 0.099,
ECE 0.009. Mouse: Brier 0.166, ECE 0.100 — the calibration does not transfer
across species even though the ranking largely does, which is why the score is
labelled human-genomic.

The point of the table is the second column. Telling observed G4-seq windows
from random genome is easy and G4Hunter already does it; telling them from
*other canonical motifs that were not observed* is the discrimination a genomic
scanner has to have, and there the rule-based scores are near chance while this
scanner holds up, including on a genome it never saw.

## Known failure mode: isolated motifs in G-poor context

The scanner learned what G4-seq observes, which is largely G-rich
neighbourhoods. A single canonical motif implanted in an AT-rich background
scores low: a human telomeric repeat in random sequence scores 0.04. In a
matched experiment, implanting `GGGTTAGGGTTAGGGTTAGGG` into 300 random genomic
windows moved the median score from 0.067 to 0.064, while implanting the MYC
NHE III1 motif moved it to 0.98.

So `scan()` does not rely on the scanner alone. Any canonical motif with
|G4Hunter| >= 1.2 that no scanner region covers is reported as its own region
with `call_basis: "canonical_motif_rescue"` and the scanner score of the window
centred on it. Regions found by the model carry `call_basis: "g4seq_scanner"`.
The distinction is in the output; nothing is silently merged.

## Use

```bash
quadcond assets fetch --name g4seq_scanner
quadcond genome-scan --fasta promoters.fa --out regions.csv --k 140 --na 10 --ph 7.4 --temperature 37
```

Each row is a region, its `g4seq_score`, its `call_basis`, one canonical motif,
and that motif's `g4_tm` and `g4_topology` under the buffer you asked for — the
division of labour the benchmark argues for: G4-seq-trained model to find the
windows, biophysical heads to say what happens to them in a buffer.

Over HTTP: `POST /scan/genome` with `{"sequence": ..., "step": 25,
"threshold": 0.5, "k": 140, ...}`, capped at 20 kb per request.

For whole chromosomes, use `genomewide/` instead: it enumerates motifs directly
and is built for resumable multi-core runs.
