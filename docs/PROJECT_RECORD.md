# QuadCond and AENNA-3D — project record

**Software release:** v0.5.1
**Record compiled:** 28 September 2026
**Author:** Auroni Deep · Kusuma School of Biological Sciences, IIT Delhi

This is a working record of what has been built, what has been measured, and what
is genuinely unresolved. It is written to be read by someone who has not been in
the room, including a reviewer.

---

## 1. Read this first: the fold-provenance problem

The most serious open problem is not a missing dataset. It is this:

> The exact 398,375-record atlas the shipped v0.5.1 models were trained on is not
> published, and the fold assignments used at training time were never persisted.

Only `atlas_core.db` (5.4 MB, a subset) ships. The full `atlas.db` (354 MB) has a
recorded content fingerprint but no release URL. The model sidecars record
*summary* metadata about the grouping — number of groups, grouping variable, fold
count, seed count — but never the assignment itself. Grouped cross-validation was
treated as a procedure rather than as an artifact, and a procedure is reproducible
only if its inputs and its determinism are both guaranteed.

### This has now been measured, and replay does not reproduce the split

`quadcond/models/foldmanifest.py` reconstructs fold assignments by replaying
`build_folds` on rows read today, and records an explicit
`assignment_provenance` flag distinguishing `reconstructed` from `historical`.
It currently reports `reconstructed`, and its own documentation is careful to say
that a matching row fingerprint is "evidence, not proof, of historical fold
identity".

Comparing the replayed manifest against the group counts recorded in the shipped
sidecars settles it:

| Head | Rows (sidecar / replay) | Groups (sidecar / replay) | Reproduced |
| --- | --- | --- | --- |
| g4_fold | 6,744 / 6,744 | 1,619 / 1,624 | no |
| g4_topology | 1,005 / 1,005 | 617 / 618 | no |
| g4_tm | 2,274 / 2,274 | 400 / 407 | no |
| im_fold | 320 / 320 | 175 / 171 | no |
| im_pht | 160 / 160 | 85 / 85 | yes |
| im_pht_condition | 558 / 558 | 10 / 10 | yes |
| im_tm_condition | 379 / 379 | 10 / 10 | yes |

Row counts match exactly everywhere, so the data being read is right. Four of the
five **sequence**-grouped heads fail to reproduce; `im_pht`, the smallest at 160
rows, reproduces, and both **condition**-grouped heads reproduce exactly. That
pattern is what tie-breaking instability would look like — fewer rows and fewer
tied lengths mean fewer opportunities to diverge — and it locates the problem in
the sequence clustering path rather than in the data or in condition grouping.

The likely mechanism is visible in `_greedy_clusters`: it orders sequences by
`np.argsort([-len(...)])`, and NumPy's default sort is not stable. Many sequences
share a length, so tied rows are permuted differently across runs, platforms or
NumPy versions — and because the algorithm is greedy, a different traversal order
picks different centroids and yields a different number of groups. Worth
confirming, and worth fixing with `kind="stable"` plus a canonical row sort before
clustering, so that v0.6 does not inherit the problem.

Both of these are now scripted under `benchmarks/fold_provenance/`:
`compare_sidecar_vs_manifest.py` regenerates the table above from the sidecar and
the manifest, and `check_determinism.py` tests the tie-ordering mechanism directly
by clustering the same sequences twice and then shuffled. The second has not been
run — it needs numpy, so it belongs in `qcgw` or the checkout's `.venv`.

### What this does not invalidate

Being precise matters here, because the instinct is to assume everything is
affected and most of it is not.

- **The head-to-head comparison stands.** Within the benchmark run every
  competitor architecture was retrained on *the same* folds as QuadCond. The claim
  is relative ordering, and that does not depend on which partition the folds were,
  only on it being identical across tools.
- **The genome-wide result stands.** Held-out chromosomes, sampled motifs, a
  within-locus flank control. It never touches atlas folds.
- **The scanner evaluation stands.** Chromosomes 2, 8 and 17 excluded entirely from
  scanner training; the mouse set independent.
- **The positive controls stand.** Six promoter G4s scanned in full against an
  external regulatory score.

### What it does affect

**Table 1 is a record, not a reproducible claim.** Its per-head metrics are copied
from build-time sidecars, computed under a partition that no longer exists and
cannot be regenerated. They remain a truthful account of the training run. They
are not something a third party can verify.

**One sentence in the manuscript now overstates.** The draft says the recomputed
out-of-fold figures "agree with the frozen record to two decimal places" (0.670
and 0.672 for G4 melting temperature), which reads as confirmation that the
original evaluation was reproduced. It was not. Those are two *different* fold
partitions giving nearly the same answer. That is genuinely reassuring — it says
the metric is stable to the partition, which is a better thing to report than a
false reproduction — but it has to be said that way.

### What to do, in order

1. **Fix the determinism** — stable sort, canonical row order — and confirm that
   two consecutive replays now agree with each other. Until that holds, no
   partition is worth publishing.
2. **Publish the full atlas** with a fixed row order and its fingerprint. Zenodo
   accepts 354 MB without difficulty. Required regardless of anything else.
3. **Adopt the fold manifest as canonical from v0.5.2** — sequence hash, group id,
   fold, seed — shipped as an artifact with its own SHA-256, and recompute the
   reported per-head metrics under it. Table 1's numbers will shift slightly. That
   is the price of every number in the paper being regenerable from published
   files, and it is worth paying.
4. **Disclose the artifact/metric relationship plainly**: the released model was
   trained under a partition that could not be reconstructed exactly, the
   difference is documented in group counts, and all reported metrics are computed
   under the published partition. A reviewer who is told this will accept it. A
   reviewer who discovers it will not.
5. **Reword the two-decimal sentence** as described above.

### Judgement

This does not block submission to *Bioinformatics* or *Briefings in
Bioinformatics*, provided the atlas is published, the canonical partition ships,
and the manuscript is explicit about which numbers are records and which are
reproducible. It does block a clean claim of end-to-end reproducibility today, and
for a web-server venue reviewers will press on exactly this.

The encouraging part is that the hard diagnostic work is already done. The fold
manifest exists, it refuses to overclaim, and the overlap audit fails closed when
leakage cannot be established rather than assuming innocence. What remains is
mostly publishing and one honest paragraph.

---

## 2. What the system is

Three layers. Keeping them distinct is what keeps the claims honest.

**QuadCond** is the predictive core and the only part making quantitative claims.
Condition-aware melting temperature and transitional pH for G-quadruplexes and
i-motifs, with K⁺, Na⁺, Li⁺/NH₄⁺, Mg²⁺, pH, temperature, crowder and strand
concentration as explicit inputs. Twelve heads, calibrated, each governed by a
claims module that fixes what it is licensed to say. Applicability is evaluated
per head, and a query outside a head's domain is *refused*: it returns no score,
no posterior, no interval. A refusal is a different object from a flagged number,
and no display setting can turn one into the other.

**AENNA-3D** is the browser workbench — React, TypeScript, Three.js — with a
geometry builder behind it that places nucleotides at published helical
parameters: 3.3 Å rise and 30° twist between G-tetrads at a 9.7 Å C1′ radius for
parallel quadruplexes; 3.1 Å per intercalated C:C⁺ step with 6.2 Å between
successive pairs of one duplex at a 5.4 Å C1′ radius for i-motifs. It is a model
builder, not a structure predictor, and it should never be described as a third
prediction method. It corrected a real defect: the previous release applied
B-form duplex constants to every fold, which misstates a quadruplex rise by
roughly a factor of two and places the tetrad hydrogen-bond network where it
cannot be.

**AlphaGenome Atlas** supplies the regulatory axis, joined on exact chromosome,
position, reference and alternate allele. The two axes are never combined into one
number; a variant is ranked by the *lower* of its two within-request percentile
ranks, which is high only when both are.

---

## 3. What was measured against other tools

Eighteen published G-quadruplex and i-motif tools, on tasks where a shared target
genuinely exists. One rule throughout: where a released model had been trained on
the evaluation data, its score is reported as an in-sample upper bound *and* the
same architecture is retrained on identical grouped folds. The retrained figure is
the comparable one.

| Task | QuadCond | Best comparable competitor |
| --- | --- | --- |
| G4 melting temperature, R² (grouped) | **0.672** | 0.478 (G4STAB architecture, retrained) |
| Buffer response, median per-sequence ρ | **0.896** | 0.800 (G4STAB retrained); 0.000 for every buffer-blind tool |
| Single-substitution ΔTm, ρ (n = 1,020) | **0.537** | 0.406 (Δpqsfinder); 0.036 for ΔG4Hunter |
| Single-substitution ΔTm, sign accuracy ≥ 2 °C | **0.809** | 0.483 for ΔG4Hunter — chance |
| i-motif folding, AUROC | **0.851** | 0.620 (iM-Seeker) |
| i-motif transitional pH, ρ (shared 137) | **0.730** | 0.465 (iM-Seeker architecture, retrained) |

Two findings generalise beyond this tool.

Released models evaluated on their own training collections reach figures no
out-of-sample protocol reproduces — 0.835 against 0.478 for one architecture on
the same data. And the choice between grouped and random splitting moved
QuadCond's own head by 0.15 R² (0.672 grouped, 0.826 random). Neither observation
criticises the tools concerned; together they are why figures quoted from separate
papers cannot be compared, and why this benchmark retrained rather than quoted.

The single-substitution result carries the most practical weight. Variant
annotation for quadruplexes currently rests on differences in a buffer-blind
motif score, and against 1,020 measured pairs that difference is at chance for
sign and near zero for rank.

---

## 4. The failure that was reported, and what came of it

On genomic windows the folding heads performed badly: AUROC 0.29–0.39 against
sequence-matched unobserved negatives, worse than a plain motif score. A model
trained on short synthetic oligonucleotides does not transfer to 124-nt genomic
windows.

Two responses. The folding heads now **refuse** genome-length input and name the
scanner instead. And a purpose-built window scanner was trained on G4-seq K⁺ data
with chromosomes 2, 8 and 17 held out entirely.

| Evaluation set | QuadCond scanner | G4Hunter | G4mismatch | G4detector (random-neg) |
| --- | --- | --- | --- | --- |
| Held-out chr, random negatives | 0.942 | 0.947 | 0.965 | 0.970 |
| Held-out chr, PQS-matched negatives | **0.928** | 0.672 | 0.849 | 0.818 |
| Mouse, PQS-matched negatives | **0.775** | 0.554 | 0.636 | 0.623 |

The second and third rows are the informative ones: PQS-matched negatives remove
the G-richness contrast, which is most of what a motif score measures. Taking the
lower of each tool's two scores as a summary of robustness to the choice of
negative set, the scanner is the best of the field. Both convolutional G4-seq
models perform well on whichever negative set resembles their own training
negatives and considerably worse on the other, which is a reason to report both
splits rather than one.

---

## 5. Variant triage, and a bug worth recording

Substitutions that abolish a motif have no ΔTm — there is no motif left to
evaluate. They were dropping out of the structural ranking entirely: the most
severe structural class, unranked. They now take the top structural rank with
their basis recorded as motif loss rather than as a magnitude.

The promoter positive controls showed this was not cosmetic. Across six promoter
G4s — MYC Pu27, KIT c-kit1 and c-kit2, VEGFA Pu22, BCL2 Pu39, hTERT hT21 — 633
substitutions joined to a regulatory record with **zero unclassified**: 65 high on
both axes, 177 structural only, 94 regulatory only, 297 neither. 111 abolish the
motif, and two of the six top-ranked variants across the panel are in that class.
A length-matched chr22 window with no quadruplex motif produced 183 substitutions
and no structural delta at all, which is the expected behaviour and the reason it
was included.

---

## 6. The genome-wide integration test

The question the earlier draft flagged as unrun: across the genome, do structural
disruption and predicted regulatory effect agree?

23 chromosomes, 150 sampled motifs each, 109,125 Atlas records with zero failed
queries, 2,707,107 scored substitutions.

**They do not.** Motif-destroying substitutions carry no excess predicted
regulatory effect over substitutions in the flanking window of the same locus.
Adjusted for CpG context, substitution type, GC and chromosome, and clustered by
motif, the G4 coefficient is −0.59 percentile points (95% CI −1.43 to +0.25,
p = 0.17) and the i-motif coefficient −0.97 (−1.54 to −0.39). Rank correlation
between predicted destabilisation magnitude and impact score is 0.018 and 0.013.

Three things make this a bounded null rather than an absence of power.

The **control was fixed**. The distant GC-matched control turned out to sit
*above* the flanking control (+4.42 percentile points, p = 3.6 × 10⁻⁷), at GC
0.703 against 0.769 and CpG 11.4% against 5.9% — because a high-GC motif that
cannot find a matched partner within 2–20 kb contributes no control at all, so the
control set is drawn from the easier loci. Every number computed against it was
reading that selection effect.

The **readout was shown to respond**. Measured on flanking substitutions alone and
on the same percentile scale, a substitution at a CpG shifts the impact score by
+9.67 percentile points for G4 loci and +10.11 for i-motif loci. Substitution type
separates strongly. CpG and substitution type explain 0.69% of the rank; structural
class explains 0.068%, adding 0.056% over composition alone.

The **bound is tight**. Destroying a G4 shifts the impact percentile by at most
1.43 points in either direction; an i-motif, at most 1.54.

### What this licenses saying

The claim is about the model, not about biology: *AlphaGenome's variant-effect
predictions carry no signal about G4 or i-motif structural disruption beyond what
local sequence composition already explains.* This design cannot separate "G4 loss
has no regulatory consequence" from "a model trained on sequence-to-function data
has no representation of a non-canonical secondary structure", and the second is
at least as plausible.

That independence is the premise the 2×2 triage encodes. Ranking by the lower of
two percentile ranks treats the axes as separate lines of evidence rather than as
corroboration, and a near-zero incremental contribution of structural class is
direct evidence for that. The corollary is a negative recommendation with some
force: a predicted regulatory impact score should not be used as a proxy for
structure-mediated regulatory effect, in either direction.

Testing the biology needs a readout that is not another sequence model. Fine-mapped
eQTL credible sets, MPRA measurements over quadruplex variants, allele-specific
quadruplex sequencing, or allele-specific expression from phased long-read
transcriptome data would each support the same paired design with a measured
outcome.

---

## 7. The manuscript, and the verification discipline

The draft was rewritten for v0.5.1 rather than patched. The v0.4.9 version said a
systematic comparator benchmark was outstanding and made no performance claim;
both statements are gone. Four new Results sections — the head-to-head, the
genomic failure and the scanner, the positive controls, the genome-wide test —
three generated tables and four new figures. The previous draft is kept unchanged.

The part worth pointing a supervisor at is `verify_numbers.py`. It holds a list of
quantitative claims appearing in the rendered manuscript and recomputes each one
from frozen result files; a claim that cannot be recomputed fails the run. It
currently verifies 36 claims clean. It caught two rounding errors and one
arithmetic mistake during drafting. The manuscript builder itself fits nothing —
it copies frozen numbers into place.

Three further pieces exist to prevent avoidable mistakes rather than to produce
results. **`benchmarks/ingest_audit`** measures how much of a candidate dataset
overlaps the atlas *and* each benchmark evaluation set before anything is
retrained, because a curated literature compilation could quietly make the
headline figures partly in-sample; tested against a set known to overlap
completely, it reports 100%. **`benchmarks/external_validation`** holds panels of
measurements published *after* the model artifact was frozen and refuses to score
a panel without a publication date or with a mismatched model hash. **The claims
module** keeps evidence semantics attached to every prediction, so a number lifted
into a figure still says whether it came from a melting curve or an antibody peak.

---

## 8. What remains

**Blocking submission**

1. Publish the full atlas and settle the fold question (§1).
2. Create the v0.5.1 GitHub release with its four asset files, then a Zenodo DOI.
3. A public server URL — the workbench currently runs on localhost, and for a
   web-server venue a reviewer must be able to open it without installing
   anything. Institutional hosting already exists for QuaDB.
4. G4detector's author list; the peptide RNA i-motif citation. Journal and DOI are
   verified for the first, not the author list.
5. Corresponding-author contact line, funding, CRediT roles, conflict declaration.

**Highest-value additions, ranked**

1. **Six CD melting curves** on top-ranked promoter variants. This is the single
   change that most moves which journal will take the paper, it is within the
   lab's existing capability, and nothing else on this list is close.
2. **Leave-one-buffer-out validation.** The sharpest question a reviewer can ask
   is whether the buffer-response advantage is physics or fitting of the training
   table's buffer columns. A day of work.
3. **The eQTL test.** Fine-mapped GTEx credible sets, same paired within-locus
   design, a measured outcome. It can distinguish the two readings of the
   genome-wide null.
4. **Release the 1,020-pair ΔTm set** as a standalone benchmark. Benchmark sets
   get cited for years.

**Queued, deliberately not done yet**

G4All is the substantial new dataset, and its value is not the 1,243 melting
temperatures but the 913 *measured* non-G4 sequences, which would let the folding
classifier be evaluated against real negatives instead of dinucleotide shuffles.
It is not ingested, because ingesting a literature compilation before an overlap
audit is the fastest way to contaminate the benchmark the paper rests on. The
audit exists; run it first, and describe the addition by how much of it is new
rather than by its headline row count.

The 50 plasma-derived oligonucleotides measured in a near-physiological
mixed-cation buffer and the FAM230E/F K⁺ titration are declared as external
validation panels with full provenance and no rows: both were published after the
artifact was frozen, which is what makes them usable that way. The sequences still
have to be extracted from the supplements by hand.

---

## 9. Overall assessment

The tool is publishable. Where a shared target exists, the condition-aware heads
lead the comparable field, and the single-substitution result has immediate
practical consequences for how quadruplex variants are annotated. The genomic
failure was found, reported, and answered with a component that holds up on
held-out chromosomes and across species. The integration hypothesis was tested
properly and returned a bounded null that supports the design it was meant to
validate.

The gap between where this is and a strong paper is not analysis. It is a
published atlas, a resolved fold question, a live server, and six melting curves.
