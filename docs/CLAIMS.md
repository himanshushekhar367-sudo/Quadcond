# What each head is allowed to claim

**7 biophysically anchored · 3 genomic-proxy · 2 derived/predicted auxiliary**

QuadCond ships twelve heads. They are not twelve versions of the same thing, and
the difference is not accuracy — it is what was measured. Every prediction
carries three fields that say so:

| field | question it answers |
|---|---|
| `has_experimental_observation` | Was anything measured upstream of this label? |
| `target_semantics` | What is the label — `biophysical`, `genomic_proxy`, `derived`, `predicted`? |
| `biophysically_grounded` | Is the label a physical measurement of the structure this head predicts? |

Only the third licenses a folding or stability claim. v0.4.0 shipped a single
`grounded_in_measurements` field that collapsed the first and third, and the
consequence was a model card describing a CUT&Tag proxy head as
"measurement-grounded: yes" on the same day the results report called it a proxy.
The field was the bug.

---

## Biophysically anchored (7)

The label is a physical measurement of the structure the head predicts. These are
the heads whose numbers may be quoted as evidence about DNA.

### `g4_tm` — melting temperature, R² 0.670, RMSE 7.8 °C
**May claim:** the melting temperature of a G-quadruplex-forming oligonucleotide
under a stated cation composition and pH, for sequences resembling the 1,062 in
G4STAB Supplementary Table 1. This is the strongest result in the release: the
condition ablation is **+0.416 R²** on 2,274 real measurements across 261 buffers
— counting a buffer as a distinct (K⁺, Na⁺, Li⁺/NH₄⁺, Mg²⁺, pH) combination.
Releases before v0.4.3 said 255, which is the same count with magnesium left out
of the key; neither number was wrong, but the definition was never stated and the
figure was a literal in the report template rather than a query against the
atlas. It is now computed at generation time.

**Verified beyond R², out-of-fold** (`scripts/10_cation_identity_check.py`). An
R² of 0.670 would look identical whether the head reads *which* cation is present
or only how much monovalent salt there is, so the atlas was searched for sequences
measured both K⁺-only and Na⁺-only at comparable concentration and pH — 411
matched pairs over 151 sequences. Each pair is scored by a refit of `g4_tm` from
the fold that **held its whole sequence-similarity cluster out**, so no pair is
predicted by a model that saw it:

| | measured | out-of-fold | in-sample |
|---|---|---|---|
| mean ΔTm(K−Na) | +12.9 °C | **+13.6 °C** [12.0, 18.4] | +13.5 °C [11.0, 21.6] |
| sign agreement | — | **94%** [90%, 98%] | 95% [91%, 98%] |
| Pearson r | — | **0.652** [0.48, 0.81] | 0.875 [0.76, 0.96] |
| MAE of ΔTm | — | 5.3 °C | 3.1 °C |

Intervals are 95% percentile bootstrap **over sequence-similarity clusters**, not
over pairs: 411 pairs come from 151 sequences and hTelo 22AG alone contributes
187, so a pair-level interval would be several times too tight. The head
separates potassium from sodium on a sequence it has never seen, in the right
direction, at close to the right magnitude — and the held-out interval on sign
agreement excludes chance.

Memorisation is visible and quantified: it buys **+0.223 of r** and +1.5% of sign
agreement. v0.4.1 reported the in-sample column alone, which was a mechanistic
consistency check described as verification; the out-of-fold column is what the
head knows.

**May not claim:** accuracy on an individual construct from that aggregate. Eight
of those 151 sequences carry published values disagreeing with *each other* by
more than 10 °C at the same nominal buffer — hTelo 22AG alone spans −4.7 to
+32.1 °C across its matched pairs. One cherry-picked comparison can make this head
look excellent or broken; neither is informative.

**May not claim:** that a titration curve is physically shaped. Walking a K⁺
ladder from 0 to 200 mM gives a 32–45 °C span (the real condition response) with
**4 non-monotonic steps across 32** on a four-construct panel. The ensemble has no
monotonicity constraint and the 5–50 mM stretch is thin. Read endpoints, not
slopes.

**May not claim:** anything at pH 5. The full pH span is 4.0–8.0 but the middle
90% is 7.0–7.5, and K⁺ reaches 1,015 mM with a 95th percentile of 120 mM. The
predictor flags a query outside the range; it cannot flag thinness inside it.
Read the applicability table before trusting an extrapolation dressed as an
interpolation.

### `g4_topology` — parallel / antiparallel / hybrid, balanced accuracy 0.724
**May claim:** a calibrated posterior over topology for a sequence already known
to form a G4, **in K⁺**.

**May not claim:** that the sequence folds — it is conditioned on folding. And
nothing about Na⁺: the training set is K⁺-only, which the ablation shows plainly
(+0.003 from conditioning). This head is not condition-aware and does not pretend
to be.

### `im_pht` — transitional pH from sequence, R² 0.592
**May claim:** which of two i-motif-forming sequences folds at higher pH, at the
iM-Seeker reference buffer (100 mM KCl + 10 mM Na cacodylate). Grouped by
sequence similarity, so the figure is what to expect on a sequence you have not
measured.

**May not claim:** any salt response. All 160 measurements are in one buffer;
conditioning adds +0.001 R². And the 160 are variants of a small number of parent
constructs, so the effective sequence diversity is far below 160.

### `im_pht_condition` — transitional pH versus buffer, R² 0.835 (grouped by buffer)
### `im_tm_condition` — i-motif Tm versus pH and ionic strength, R² 0.853 (grouped by buffer)

**These two heads return no value for a sequence outside their two-construct
allowlist.** Asked about an unseen sequence, `im_tm_condition` returned 44.82 °C;
asked about Tel21C itself it returned 44.73 °C. Not a degraded answer on
unfamiliar input — the *same* answer, because there is no sequence discrimination
in the head to degrade. Two sequences cannot teach a model how sequence affects a
condition response; they can only teach it how two molecules behave.

v0.4.2 called this a refusal while still emitting the number beside an
`in_domain: false` flag, and the CLI printed it. Since v0.4.3 the `value`,
`interval` and `folded_fraction_at_condition` keys are **absent** — a consumer
that forgets to check the flag gets a `KeyError`, which is the correct failure.
A flag next to a number relies on every downstream reader remembering which one
wins.

**May claim:** how **Tel21C and C9 specifically** respond to salt and pH.
Folds hold out whole ionic strengths, so the number answers "does this
extrapolate to a salt concentration nobody measured". The ablation on
`im_tm_condition` is **+0.790 R²** — the starkest condition effect in the release.

**May not claim:** anything about a new sequence. There are two sequences in this
panel. Two. A condition-response model, not a sequence model, and the pair are
reported separately from `im_pht` for exactly that reason: pooling them and
grouping by sequence gave **R² −0.257**, worse than predicting the mean.

### `g4_fold` — P(G4) versus matched shuffles, AUROC 0.970, ECE 0.008
### `im_fold` — P(i-motif) versus matched shuffles, AUROC 0.870, ECE 0.055
**May claim:** discrimination between measured G4s / i-motifs and
dinucleotide-preserving shuffles of themselves, calibrated on that task.

**May not claim:** these are absolute probabilities for a genomic, transcriptomic
or aptamer population. The negatives are *generated*, not experimentally tested
non-folders, and the task is roughly balanced by construction — neither is true
of any candidate set a user will actually query. Recalibrate against your own
prevalence. `im_fold` additionally does not transfer: on the cross-assay CD panel
it separates i-motifs from negatives by **+0.024** and cannot distinguish an
i-motif from a G4.

**Where the number is served.** `/evidence` used to print these two under each
structure card as `probabilityOfFolding`, beside a `source`, a claim basis and
an applicability domain that all described `g4_tm` / `im_pht`. Nothing in the
response said the number came from a different head, and every cue around it
invited reading it as P(this sequence folds). They now occupy their own
`discrimination` block carrying their own head name, strand, strand sequence,
claim basis, applicability and the sentence above; the cards carry no
probability of any kind. A number is only as honest as the provenance printed
next to it, and a correct value under the wrong provenance is a wrong answer.

---

## Genomic proxy (3)

The label is **antibody occupancy at a genomic locus**, not folding. Both heads
have `has_experimental_observation = true` and
`biophysically_grounded = false`, and both are printed with a
`[GENOMIC PROXY — not a folding probability]` marker by the CLI.

### `im_fold_genomic` — AUROC 0.779, ECE 0.007
### `g4_fold_genomic` — AUROC 0.765, ECE 0.008

**May claim:** that a sequence resembles the canonical motifs enriched under
reproducible iMab (or BG4) CUT&Tag peaks in live HEK293T cells, as against
composition-matched motifs drawn from dinucleotide shuffles of the same peak
windows. On the cross-assay CD panel `im_fold_genomic` separates the authors'
i-motifs from their negatives by **+0.244** and scores their G4 controls at
0.304 — the discrimination `im_fold` cannot make. That gain came from sequence
diversity (8,482 independent loci), not from a better model.

**May not claim:** P(folds), in cells or anywhere else. Three separate reasons,
each sufficient on its own:

1. **CUT&Tag measures occupancy, not structure.** An antibody bound at a locus is
   not a melting transition.
2. **The antibody is contested.** iMab specificity and the contribution of
   accessible chromatin to targeted CUT&Tag signal are both open questions in the
   literature. A score built on that signal inherits the uncertainty.
3. **No buffer was measured.** The attached condition is a nominal guess at
   intracellular chemistry with **all eight fields flagged imputed**, so these
   heads cannot carry a condition response and none is claimed.

`g4_fold_genomic` exists mostly as a control: the two heads share one pipeline,
so a gap between them says something about the antibodies and the motif grammar
rather than about the model.

### `locus_peak_overlap_state` — four-class peak overlap in a 201-nt window

Renamed in v0.4.3. It was `locus_state`, described as "joint occupancy", and the
experiment does not observe that.

**What the design actually is.** Zanin et al. ran iMab and BG4 CUT&Tag as
**parallel reactions on separate aliquots** of one HEK293T population. A window
called "both" had reproducible peaks from both antibodies *in that population*.
Nothing anywhere in the design observes two structures on the same DNA molecule,
and "joint occupancy" invites precisely that reading.

**May claim:** a calibrated posterior over how a **201-nt genomic window** was
classified by those two experiments — peaks from both, from BG4 only, from iMab
only, or from neither. Trained on 35,314 merged loci: 13,156 both, 15,088
G4-only, 7,070 iM-only.

**May not claim,** and the numbers matter more than the prose:

| | |
|---|---|
| window length | **201 nt in every training row**, so an oligonucleotide query is out of domain by construction, not by degree |
| balanced accuracy | **0.450** against a 0.250 four-class floor |
| iM-only recall | **0.156** — the class users ask about most is the one it recovers worst |
| the "neither" class | a composition-matched **shuffle**, not an observed unoccupied locus: three classes are observations, the fourth is generated |

It is an exploratory genomic-resemblance score. Not co-occupancy, not
competition, not free energy, and not a folding probability. It remains the best
available evidence that *the sequence at a site carries information about which
antibodies bound there*, which is a real and modest result.

## Derived / prior-model auxiliary (2)

No measurement of any kind stands behind these labels.

### `g4_tm_distilled` — R² 0.963, `target_semantics = predicted`
**May claim:** to reproduce the G4STAB ensemble's Tm response to cation
composition. Useful as a condition-response prior and as a benchmark target.

**May not claim:** anything about DNA. It is a model fitted to another model's
output, and it inherits every G4STAB bias. Its R² of 0.963 is agreement with a
predictor, not accuracy against measurements — the head that carries the
measurements is `g4_tm`, at 0.670.

### `im_architecture` — AUROC 0.997, `target_semantics = derived`
**May claim:** essentially nothing. The positive label is "carries a canonical
four-C-tract motif" and the negatives are shuffles from which that motif was
rejected, so a near-perfect AUROC means the feature set can recompute the regex.
It exists so the i-motif pipeline is exercisable end to end.

---

## Not solved, and not claimed

- **Thermodynamic G4/i-motif competition.** `locus_peak_overlap_state` replaced
  the independence proxy with an observation, which is progress and is not the
  same thing as solving this. v0.4.2 went further and served a Boltzmann
  ensemble over {G4, i-motif, unfolded}; that has been **withdrawn**, because the
  state competing with both at physiological conditions is the duplex and it was
  not in the partition function — nor was strand concentration or stoichiometry.
  Percentages over an incomplete state space are not uncertain, they are about
  something else. The endpoint now returns strand-resolved evidence that does not
  sum to one. What is still missing is the *in vitro* experiment: one
  duplex, both strands, a condition grid, and a measurement of which structure
  wins at each point. Until that exists, no head here can say which fold a locus
  adopts in a defined buffer — only which antibodies bound it in one cell line.
  `quadcond competition` remains available for the sequences
  `locus_peak_overlap_state` has no genomic context for, but as of **v0.4.8** it
  returns nothing that combines the two structures. The joint 2x2 table, the
  signed `preference_g4_minus_im` and the interpretation sentence naming a
  favoured structure are **withdrawn**. They were labelled EXPERIMENTAL and the
  label did not change what they claimed: the table multiplied two probabilities
  as if forming a G4 and forming an i-motif at one locus were independent events
  — mutually exclusive by construction, since both need the duplex open on
  opposite strands — and its four entries summed to one, which is what makes a
  table read as a distribution over states. The preference subtracted a G4
  discrimination score from an i-motif *architecture* score: different labels,
  different tasks, different calibration sets, and a difference does not acquire
  a meaning from being computed. What remains is each strand's predictions under
  one buffer, with their own applicability and claim basis, and a
  `combination_withdrawn` block naming what was removed.

- **RNA.** Every measurement in the atlas is DNA, and the service now **refuses**
  RNA queries rather than answering them. rG4s are effectively always parallel,
  more stable at the same ionic strength, and the 2′-OH changes the loop
  energetics most of these features encode. A DNA-trained head returns a
  confident number for an RNA sequence with nothing behind it, and a warning
  beside such a number reads as a caveat on a real answer rather than as the
  absence of one.

- **Free energy from the heuristic.** The viewer's energy function was
  calibrated against 2,274 measured melting temperatures converted to folding
  free energies under a two-state van 't Hoff model. That fixed its *scale* —
  the span dropped from 78 RT to 4.1 RT, so the Boltzmann step has real dynamic
  range instead of saturating its own exponent clamp. It did not make the
  function accurate: **R² 0.074** against measured ΔG. It is therefore used only
  for folds no trained head covers (cruciform, triplex, AC-motif), and every
  such number is badged `heuristic` in the interface.

- **Tertiary structure.** The 3D geometry is an idealised model at published
  helical parameters — 3.3 Å rise and ~30° twist per G-tetrad, 6.2 Å between
  successive C:C⁺ pairs of one i-motif duplex with the partner's pairs
  intercalated at 3.1 Å. It is validated for tetrad planarity, four-fold
  occupancy, backbone step length and hard-sphere clashes, and it is not a
  prediction of atomic coordinates. Exported PDB files say so in their own
  REMARK records.

  Three G4 archetypes are drawn — parallel, antiparallel and **hybrid (3+1)**,
  the last being the form human telomeric DNA adopts in K⁺ and therefore the
  subject of the most common query this tool receives. It was previously absent
  from the geometry layer while the heads went on predicting it, and the viewer
  resolved the missing archetype to the parallel one. Where a predicted topology
  has no archetype the viewer now draws **nothing** and says so: a picture is a
  claim about which fold, and substituting a different fold answers a question
  nobody asked.
- **A general condition-aware i-motif model.** The sequence axis has 160
  measurements in one buffer; the condition axis has two sequences across many.
  Both are real; their product is not yet covered. Closing it needs a matrix of
  roughly 50–100 diverse sequences — canonical, long-loop and bulged motifs plus
  **measured non-folders** — across at least four ionic conditions with full pH
  titrations and replicates, evaluated with sequence-family, buffer-regime and
  independent-source holdouts. That is the decisive v0.5 experiment.
- **RNA, cross-species, non-canonical motifs.** The framework accommodates them;
  no training evidence in this release supports a claim about them.

---

## What the tool does with these heads

The heads are the evidence. The workflow is what a researcher does with it, and
three of its parts make claims of their own that belong here.

**Mutation effects are a different task from absolute prediction.** `/scan/mutations`
reports ΔTm and ΔpH_T from the head that owns the quantity — never a difference
between two structures' energies, which is a property of neither. The
uncertainty is computed **paired**: each ensemble member predicts wild type and
mutant, and the spread of those differences is reported, because the two errors
are strongly correlated and propagating the marginal intervals would overstate
the noise by a wide margin.

*May claim:* that this model's ensemble agrees, or does not, about the direction
and rough size of a substitution's effect, under a stated buffer, for a mutant
that still carries the motif.

*May not claim:* a calibrated effect size. **No measured wild-type/mutant panel
has been held out to validate mutation-effect prediction as its own task.** A
head's error on absolute melting temperatures does not transfer to its error on
differences, and no row here is labelled significant — there is no null
distribution to test against, and a magnitude threshold is a threshold with the
word attached.

**A condition response is measured, not assumed.** A head appears on a response
surface because its training data spanned that axis, which licenses a dependence
without demonstrating one — a composition-matched shuffle inherits its parent
row's buffer, so a classification head can show wide training variation in an
axis its label never depended on. Every series therefore reports the observed
spread against the head's own error scale. On Tel22 across K⁺ 0–150 mM, `g4_tm`
moves 34.4 °C and `g4_fold` moves 0.0008; only one of those is a response.

Heads that never saw the axis vary are named as inert with the reason, not
omitted. A head missing from a chart reads as "not applicable", which conflates
an absence of evidence with a refusal for this sequence.

**Synthetic data can never present as measurement.** A generated corpus is
tagged `evidence_tier="synthetic"`, which `target_semantics` checks before
everything else and which no downstream formatting can undo: such a head reports
`biophysically_grounded: false` and `has_experimental_observation: false`, the
headline leads with "SYNTHETIC MODEL — not for scientific use", and the service
refuses the artifact unless explicitly started for it. Mixed tiers poison rather
than dilute the claim, because nobody can say which rows drove which prediction.

---

## The one-sentence positioning

QuadCond is a provenance-aware G4/i-motif framework that demonstrates
reproducible condition effects on measured stability while keeping genomic
proxies, derived labels and direct biophysical evidence scientifically separate.
It is not yet a validated general folding predictor for either structure.
