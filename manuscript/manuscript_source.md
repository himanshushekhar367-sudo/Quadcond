# QuadCond and AENNA-3D: a benchmarked, condition-aware workbench for G-quadruplex and i-motif evidence

Auroni Deep<sup>1,†</sup>, Himanshu Shekhar<sup>1,†</sup>, Shilpi Minocha<sup>1</sup>, Vivekanandan Perumal<sup>1</sup> and Saran Kumar<sup>1,*</sup>

<sup>1</sup>Kusuma School of Biological Sciences, Indian Institute of Technology Delhi, Hauz Khas, New Delhi 110016, India

<sup>†</sup>These authors contributed equally to this work.

<sup>*</sup>To whom correspondence should be addressed. Saran Kumar. [AUTHOR QUERY: institutional email and telephone for the corresponding-author line.]

Running title: Condition-aware G4 and i-motif evidence

Keywords: G-quadruplex; i-motif; DNA secondary structure; melting temperature; transitional pH; benchmarking; applicability domain

Author working draft, 18 September 2026, describing software release v0.5.1. Internal model results are retained from the frozen v0.4.7 sidecars and were not refitted. The head-to-head benchmark, the genomic scanner and the genome-wide variant analysis are new to this revision and were computed for it. Bracketed deployment fields and citations marked CITATION TO ADD require completion before submission. Public deployment and archival DOI remain outstanding.

## Abstract

G-quadruplexes and i-motifs form in guanine-rich and cytosine-rich DNA, and whether either forms depends on the buffer as much as on the sequence. We present QuadCond and AENNA-3D, a Python prediction service and a browser workbench that hold different kinds of evidence apart while placing them in one strand-aware interface. Every output states the strand it was computed on, the training evidence behind its target, its applicability and the model that produced it; a query outside a head's domain is refused and returns no number at all. We benchmarked the service against eighteen published G-quadruplex and i-motif tools on shared tasks, retraining competitor architectures on identical grouped folds wherever a released model had seen the evaluation data. Under grouped cross-validation, G-quadruplex melting temperature reached an R² of 0.672 against 0.478 for a retrained G4STAB architecture; i-motif folding reached an AUROC of 0.851 against 0.620 for iM-Seeker; and on 1,020 measured single-substitution pairs the predicted change in melting temperature correlated at a Spearman ρ of 0.537 with the correct sign in 81% of pairs of at least 2 °C, where a ΔG4Hunter score was at chance (ρ 0.036, sign accuracy 0.48). Within-sequence buffer response, which no buffer-blind tool can express, reached a median per-sequence ρ of 0.90. The folding heads fail on genomic windows, which we report and then address: a purpose-built window scanner trained on G4-seq reaches an AUROC of 0.928 against sequence-matched unobserved negatives on held-out chromosomes, where G4Hunter reaches 0.672, and it is the only tool above 0.92 on both negative sets. Joining the structural axis to AlphaGenome Atlas variant impact scores across 23 chromosomes and 2.7 million scored substitutions, motif-destroying variants carried no excess predicted regulatory effect over same-locus flanking variants once sequence composition was controlled; the bound is 1.4 percentile points, against a 9.7-point shift for CpG context measured on the same scale. The two axes are therefore non-redundant, which is the premise the workbench's triage is built on. Web server: [VERIFIED PUBLIC HTTPS URL].

## Introduction

Guanine-rich DNA can fold into G-quadruplexes (G4s), and the complementary cytosine-rich strand can fold into i-motifs held together by intercalated C:C⁺ base pairs. Antibody imaging and antibody-based sequencing place both in human cells and link G4s to regulatory chromatin and transcription [@biffi;@zeraati;@hansel]. i-Motifs have been proposed as regulatory switches at promoters, and their folding responds to the cellular environment: pH and cation composition set whether a given C-rich tract folds at all, and cellular polyamines can disrupt folded i-motifs outright [@im_switch;@polyamine]. G4 stability depends on which monovalent cation is present and at what concentration, not merely on the presence of a cation [@cations]. A sequence that folds in one buffer is therefore not a sequence that folds, and the same measurement made under two salt conditions can support opposite conclusions.

That dependence is what most software leaves out. QGRS Mapper, G4Hunter and pqsfinder identify or score potential G4-forming sequences from the sequence alone [@qgrs;@g4hunter;@pqsfinder]. G4-iM Grinder searches G4 and i-motif patterns together, including higher-order arrangements [@grinder]. QuaDB retrieves putative quadruplex and i-motif sequences by transcript or gene identifier for rapid lookup [@quadb]. Machine-learning predictors extend the sequence-only question: DeepG4 predicts cell-type-specific active G4 regions from sequence context, and G4Boost estimates folding state and thermodynamic stability from sequence-intrinsic features [@deepg4;@g4boost]; convolutional models trained directly on G4-seq signal, including G4detector [CITATION TO ADD: G4detector, IEEE/ACM Trans. Comput. Biol. Bioinform., doi:10.1109/TCBB.2021.3073595 — DOI and journal verified, author list not yet confirmed] and G4mismatch [@g4mismatch], predict observed quadruplex propensity genome-wide. Two recent tools do take a condition as input — G4STAB predicts G4 melting temperature from sequence and salt concentration, and iM-Seeker scores i-motif folding status and strength for C-rich sequences [@g4stab;@im_web;@im_method] — and G4ShapePredictor predicts folding topology [@g4shape]. Curated resources record the experimental record itself: ONQUADRO holds experimentally determined quadruplex structures, G4LDB catalogues G4 and i-motif ligands, and a recent high-throughput study measured i-motif stability across a designed condition landscape [@onquadro;@g4ldb;@duvma]. Mutation design and variant annotation are served separately again, by G4Killer [@g4killer] and by G4SNVHunter [@g4snvhunter].

Each of these answers a well-posed question, and for most pairs of them no common score exists. That has had an unfortunate consequence for the field: tools are rarely compared, and where they are, a released model is often evaluated on data it was trained on. Two further distinctions are easy to lose. A G4 and an i-motif at one duplex position sit on opposite strands, so a score computed on the strand the user happened to paste may be a score about the wrong molecule. And a prediction can be outside the range its model was trained on without anything in the output saying so — the applicability-domain problem long recognised in structure–activity modelling [@applicability].

This paper does three things. It describes QuadCond and AENNA-3D, which place condition-dependent estimates, strand assignment, applicability and comparison workflows in one interface. It reports a head-to-head benchmark against eighteen published tools on the tasks where a shared target does exist, retraining competitor architectures on identical folds wherever a released model had seen the evaluation data, and reporting the tasks on which QuadCond loses as well as those on which it wins. And it uses the resulting variant workflow to ask, genome-wide, whether structural disruption and predicted regulatory effect agree.

## Materials and methods

### Data organization and evidence semantics

The frozen training atlas holds 398,375 records: 254,761 labelled predicted, 88,643 derived and 54,971 experimental. These are records, not independent sequences, and not 398,375 physical measurements. The experimental tier also carries genomic assay labels, so the tier alone does not separate a melting curve from an antibody peak. A second field, label class, distinguishes biophysical observations from genomic proxies, catalog annotations and synthetic labels. Both fields have to be read before an output is given a scientific meaning, and both travel with every response and every export.

The principal biophysical sources are the topology collection used by G4ShapePredictor, the measured stability collection used by G4STAB, the i-motif transitional-pH collection associated with iM-Seeker, and condition-dependent i-motif measurements from the high-throughput stability study [@g4shape;@g4stab;@im_method;@duvma]. Their training subsets contain 1,005 topology records, 2,274 usable G4 melting-temperature records, 160 sequence-oriented i-motif transitional-pH records, 558 condition-dependent transitional-pH records and 379 condition-dependent melting-temperature records. The last two subsets cover two constructs, C9 and Tel21C. The wider atlas also holds predicted stability values, dinucleotide-preserving shuffled controls [@shuffle] and derived architecture labels, none of which is treated as experimental validation.

Genomic proxy models use sequence windows from the HEK293T i-motif mapping data in GSE220882 [@zanin]. Their outputs are associations with a training label. Reports on iMab specificity have reached different conclusions under different experimental conditions [@imab_critique;@imab_response], and the contribution of accessible chromatin to G4-targeted CUT&Tag signal has been contested [@cuttag;@cuttag_response]. These are unsettled assay-level questions, and they are the reason a proxy score is never presented as a probability that a molecule is folded. A four-class locus model describes overlap of G4 and i-motif peak labels in 201-nt windows; because the two assays used separate aliquots, the class labelled "both" reports population-level co-localisation and observes nothing about two structures on one molecule.

### Model construction and internal evaluation

QuadCond uses task-specific gradient-boosted tree models built with XGBoost [@xgb]. Sequence features cover nucleotide composition, short k-mer frequencies, G- and C-tract and loop descriptors, and sequence propensity summaries. Where the training source supports it, models also take buffer variables, their transformations, and sequence–condition interactions. The frozen package contains 12 prediction heads: seven biophysical, three genomic proxy and two auxiliary heads built on predicted or derived labels. Sharing one model artifact does not make their targets interchangeable, and the per-head documentation follows the model-card practice of stating intended use and evaluation context alongside the numbers [@modelcards].

Sequences are grouped for evaluation on normalized 4-mer representations. Collections below 6,000 items use greedy cosine-similarity grouping at a threshold of 0.90; larger collections use MiniBatchKMeans with at most 4,000 clusters. Grouping reduces sequence similarity across folds; it does not guarantee that every related sequence lands in one fold. Five-fold evaluation is recorded for each head, with five fitting seeds for most heads and three for the two auxiliary heads. The condition-dependent i-motif heads instead hold out 10 condition groups across five folds, which measures transfer across those condition groups for the same two constructs and not transfer to unseen sequence families.

Binary classification uses cross-fitted probability calibration where recorded, and multiclass classification uses temperature scaling; both follow standard post-hoc calibration practice for tree ensembles and multiclass scores [@calibration;@isotonic]. Calibration is evaluated against the supplied labels, and it cannot turn shuffled negatives or antibody labels into a physical folding probability — a calibrated score on a constructed task is calibrated for that task. Regression intervals use out-of-fold residuals as an empirical error scale. The metadata separately reports coverage from 20 partitions in which half the groups set a width and the other half measure coverage; we report that held-out-group coverage rather than coverage measured on the same residual pool that set the width. These intervals do not establish conditional coverage at a particular sequence or buffer, and they are not validated intervals for differences between mutations.

### Benchmark against published tools

Comparisons were run only on tasks where the tools involved answer the same question, and each comparison uses one fixed evaluation set with one target. QuadCond contributes grouped out-of-fold predictions throughout, so no QuadCond number in the benchmark is an in-sample number.

Competitor tools were obtained from their published distributions and run on the identical inputs. Where a released model had been trained on the evaluation data — which is the case for G4STAB on the melting-temperature collection, for iM-Seeker on the transitional-pH collection, and for the three G4ShapePredictor estimators on the topology collection — we report the released model's score marked as in-sample and additionally retrain the same estimator and hyper-parameters on QuadCond's grouped folds. The retrained figure is the comparable one; the in-sample figure is retained as an upper bound on what that architecture can express, and both are reported so that the gap between them is visible. To make the protocol difference itself visible, we also evaluated QuadCond's melting-temperature head under random rather than grouped five-fold splitting, matching the protocol under which G4STAB was published.

Tasks were: G4 folding against dinucleotide-preserving shuffles [@shuffle], n = 5,214; G4 topology, three classes, n = 1,005; buffer-resolved G4 melting temperature, n = 2,274; within-sequence buffer response, 128 sequences with 956 measurements across buffers; single-substitution ΔTm for 1,020 pairs measured in the same buffer; i-motif folding against dinucleotide shuffles, n = 304; and i-motif transitional pH, n = 160, with a 137-sequence subset on which iM-Seeker returns a motif call. Confidence intervals are 1,000-replicate bootstrap intervals resampled over sequence groups rather than over rows, so a group contributing many measurements cannot narrow an interval by itself. Buffer-blind tools produce a single score per sequence and are constant across buffers by construction; this is recorded as a zero rather than omitted, because it is the quantity the comparison is about.

Genomic performance was evaluated separately on G4-seq K⁺ data for the human genome, with two negative sets: random genomic windows matched on length, and unobserved windows matched on predicted quadruplex-sequence content, which is the harder and more informative control because it removes the trivial G-richness signal.

### Genomic window scanner

The fold heads are trained on oligonucleotides and fail on genomic windows (Results). Release v0.5.1 therefore refuses those heads on genome-length input and routes the question to a separate scanner. The scanner represents a 124-nt window by the QuadCond feature set computed on both strands with the G-richer strand placed first, plus ten window-level summary statistics, giving 266 features; strand symmetry means the score does not depend on which strand was supplied. It was trained on G4-seq K⁺ windows with sequence-matched unobserved negatives, holding out chromosomes 2, 8 and 17 entirely, and is evaluated on those held-out chromosomes and, separately, on mouse G4-seq windows as a cross-species test. A canonical-motif rescue pass reports any strong G4Hunter motif not covered by a called region, so a high-confidence motif is never silently dropped because the window model missed it; rescued regions are labelled with their call basis.

### Variant workflow and the regulatory axis

The mutation workflow evaluates single-base substitutions against a frozen reference and reports the difference together with the spread across paired estimators — each estimator's own reference-minus-mutant difference, rather than a subtraction of two marginal intervals. Substitutions that abolish the motif have no ΔTm to report, because there is no motif left to evaluate; these are assigned the highest structural rank with their basis recorded as motif loss rather than as a magnitude, so the most severe structural class is not excluded from the ranking for want of a number.

The regulatory axis uses AlphaGenome Atlas variant impact scores [@alphagenome_atlas], joined on chromosome, position, reference and alternate allele, so the join is exact rather than overlap-based. Three interchangeable sources are supported: the published Tabix bundle, a user-supplied score table, and the Atlas API. The two axes are never combined into one number; a variant is ranked by the smaller of its two within-request percentile ranks, which is high only when both are high, and its quadrant against stated thresholds is reported alongside.

### Genome-wide test of structural and regulatory concordance

For each of chromosomes 1–22 and X we identified canonical G4 and i-motif motifs, sampled 150 per chromosome, and enumerated every single-nucleotide substitution in each motif and in a 100-nt window flanking it on each side. Structural consequence was computed for every substitution and classified as motif loss, destabilising, moderate, neutral or stabilising. Regulatory effect was the absolute Atlas variant impact score, queried live for the motif and its flanks in a single call so that both sets share the query batch.

Two control sets were used. A distant control matched each motif on GC content 2–20 kb away, and a flanking control took the 100-nt window on each side of the motif itself. The flanking control is the primary one: it holds the locus, the regulatory neighbourhood and the regional score distribution fixed, and differs from the motif only in motif membership. The distant control is reported because it is what a naive design would use and because its imbalance is itself a result.

The primary test is a within-locus paired Wilcoxon signed-rank comparison of motif-destroying against flanking substitutions, one pair per locus. The primary adjusted model regresses the within-chromosome percentile rank of the absolute impact score on structural class, substitution type, CpG context, GC content and chromosome, with standard errors clustered by motif and the flanking class as reference; a rank outcome is used because top-percentile membership is rare among motif substitutions and a logistic fit on that many events separates rather than estimates. Covariate balance across classes is reported in full.

Because a null result is only interpretable if the readout can respond to anything at this scale, we added assay-sensitivity controls computed on the flanking substitutions alone, so that no control depends on a motif call: the effect of CpG context and of substitution type on the same percentile scale, and the share of rank variance explained by sequence composition against that explained by structural class.

### Prediction service and applicability

The Python service exposes sequence evidence, condition scans, mutation scans, batch processing and genome scanning. The input strand and its reverse complement are both represented 5′ to 3′. G4 and i-motif heads are routed to the strand carrying the motif each head was trained on, and every result names its strand; the reverse complement is not intrinsically the i-motif strand, and on a C-rich input the assignment reverses. Conditions supplied by the caller are recorded separately from imputed defaults, so a defaulted 100 mM K⁺ cannot later be read back as a measured one.

Applicability is evaluated per head, and a result is supported, outside the head's domain, or refused. An unsupported estimate can be inspected by explicitly requesting extrapolation. A refusal is a distinct response type carrying no score and no posterior, and no display setting in the browser can reveal one, because there is nothing behind it to reveal. From v0.5.1 the folding heads refuse genome-length input outright and name the scanner instead. Each response also carries evidence semantics, the applicable warnings and model provenance. Readiness and inference resolve the same explicitly configured model and atlas paths and verify the selected assets against their manifests.

### Comparison workflows and structural display

Condition sweeps vary one input and hold the rest of the requested buffer fixed. Response summaries cover in-domain points only; with extrapolation disabled neither the line nor the point markers are drawn for unsupported points, and refused points are absent regardless of the display setting. Batch input accepts FASTA or one sequence per line. Empty records are kept in the accounting as named exclusions, repeated names remain separate entries, and internal row indices keep caller-supplied duplicate indices from collapsing the count. JSON and CSV exports retain input context and provenance.

AENNA-3D is built with React, TypeScript and Three.js and shows evidence beside structural archetypes. It is a model builder, not a structure predictor: it places nucleotides on the geometry each fold is known to adopt, at published helical parameters — 3.3 Å rise and 30° twist between G-tetrads at a 9.7 Å C1′ radius for parallel quadruplexes, and 3.1 Å per intercalated C:C⁺ step with 6.2 Å between successive pairs of one duplex at a 5.4 Å C1′ radius for i-motifs — so that a rendered picture and an exported coordinate file are the right shape and the right size. It does not refine, does not score and does not claim to reproduce any deposited structure, and every export is labelled as a model. The preceding release applied B-form duplex parameters to every fold, which misstates a quadruplex rise by roughly a factor of two and places the tetrad hydrogen-bond network where it cannot be; correcting that is why the geometry is described here at all. Selecting an evidence card selects the matching archetype. Figure 1 summarizes the workflow and the boundaries between evidence types.

## Results

### Retained internal model performance

Model results in this subsection come from the frozen v0.4.7 sidecars, retained unchanged; nothing was refitted for v0.5.1. Table 1 gives the four regression heads with measured targets. The G4 melting-temperature head reached an R² of 0.670, a root mean squared error (RMSE) of 7.76 °C and a mean absolute error of 5.58 °C; its 2,274 records fall into 400 sequence groups. The sequence-oriented i-motif transitional-pH head reached an R² of 0.592, an RMSE of 0.349 pH units and a mean absolute error of 0.255 over 160 records in 85 groups. Out-of-fold predictions recomputed for the benchmark agree with the frozen record to two decimal places (0.672 and 0.585).

The condition-dependent i-motif transitional-pH and melting-temperature heads reached R² values of 0.835 and 0.853. Both were evaluated with condition groups held out for C9 and Tel21C, and the two-construct scope is the central limitation on those two numbers. Mean held-out-group coverage of the empirical error widths ranged from 0.880 to 0.904 across the four regression heads, against a nominal 0.90.

The sequence-only ablation gives a target-specific comparison inside the same development framework (Figure 2). Removing condition features dropped the G4 melting-temperature R² from 0.672 to 0.256 and the median per-sequence buffer-response correlation from 0.90 to zero, since a model without buffer input cannot respond to buffer at all. For the condition-dependent i-motif heads the change was 0.835 to 0.625 for transitional pH and 0.853 to 0.063 for melting temperature. The sequence-oriented transitional-pH head barely moved, 0.592 to 0.590 — an unsurprising result for a head whose 160 training rows were measured in a single buffer, and a reminder that a near-zero ablation delta is a statement about the data rather than about the method.

The G4 topology head reached a balanced accuracy of 0.723 over 1,005 records. The G4 and i-motif binary fold classifiers reached areas under the receiver operating characteristic curve of 0.957 and 0.851 against dinucleotide-preserving shuffles of their own positives [@shuffle], which bounds how far those numbers can be read biophysically. Genomic proxy classification was less discriminative: 0.765 for G4 and 0.779 for i-motif peak association. The four-class locus model reached a balanced accuracy of 0.450, with low recall for the i-motif-only class. Supplementary Table S1 lists training size, grouping and target semantics for all 12 heads.

### Head-to-head comparison with published tools

Table 4 reports the benchmark. Three findings organise it.

First, on the tasks QuadCond was built for it leads the comparable field. For buffer-resolved G4 melting temperature, grouped out-of-fold R² was 0.672 (95% CI 0.606–0.729) with a mean absolute error of 5.56 °C and Spearman ρ 0.819. The G4STAB architecture retrained on the same grouped folds reached R² 0.478 (0.390–0.558), MAE 7.42 °C and ρ 0.702. The released G4STAB ensemble reached R² 0.835 on this collection, but it was trained on it, and we report that as an in-sample upper bound rather than as a competitive result. The protocol matters as much as the model: evaluated under random five-fold splitting, the protocol under which G4STAB was published, QuadCond's own head reaches R² 0.826 and ρ 0.905 — statistically indistinguishable from the released ensemble's in-sample figure and roughly 0.15 R² above its own grouped result. Any comparison that mixes the two protocols is measuring the protocol.

For i-motif folding against dinucleotide shuffles, QuadCond reached an AUROC of 0.851 (0.794–0.901) against 0.620 for iM-Seeker, 0.618 for a canonical C3+ regex count and 0.498 for GC content. For i-motif transitional pH, iM-Seeker's released folding strength correlates at ρ 0.966, but it was trained on these values; retrained on grouped folds over the 137 sequences for which it returns a motif call, it reached ρ 0.465 and R² 0.253, against ρ 0.730 and R² 0.484 for QuadCond on the same 137 sequences.

Second, condition response separates the field completely. Across 128 sequences measured in more than one buffer, QuadCond's median per-sequence Spearman correlation between predicted and measured Tm was 0.896, and it explained 75% of the variance in buffer-centred Tm — the component of the signal that depends on buffer alone with every sequence effect removed. The retrained G4STAB architecture reached 0.800 and 47%. Every buffer-blind tool reached exactly zero on both, not because it scored poorly but because a single score per sequence cannot vary with buffer. Against raw Tm those same tools reach modest rank correlations — G4mismatch 0.393, G4Hunter 0.340, pqsfinder 0.300, QGRS Mapper 0.197, DeepG4 0.106, G4Boost −0.006 — which reflects the sequence component of stability and says nothing about condition sensitivity.

Third, on single-substitution effects the gap is large. Over 1,020 pairs whose wild-type and mutant Tm were measured in the same buffer, QuadCond's predicted ΔTm correlated at ρ 0.537 (0.461–0.614) and recovered the correct sign in 80.9% of pairs differing by at least 2 °C. The variant scores in current use fared worse: Δpqsfinder 0.406, ΔQGRS 0.378, ΔG4mismatch 0.336, ΔDeepG4 0.315, and the ΔG4Hunter score that underlies G4SNVHunter's variant impact measure [@g4snvhunter] 0.036 with a sign accuracy of 0.483, which is chance. The retrained G4STAB architecture reached 0.331 and 0.696. This is the comparison most directly relevant to variant interpretation, and it is the one where a buffer-blind ΔG4Hunter is least usable.

QuadCond's G4 folding classifier reached an AUROC of 0.957 against dinucleotide shuffles, where the best sequence-based competitor reached 0.752 (QGRS Mapper) and G4Hunter reached 0.394. That last figure is not a failure of G4Hunter but a property of the task: a dinucleotide-preserving shuffle holds G-richness approximately constant, so a G-richness score cannot separate the classes, and a score below 0.5 indicates the residual composition signal runs the other way. We report the task because it is the one those heads were trained on, and we do not read it as evidence of genomic performance.

### Where the folding heads fail, and what replaces them

On genomic windows the folding heads perform badly, and this is the clearest negative result in the benchmark. Against G4-seq K⁺ positives with sequence-matched unobserved negatives (n = 6,002), the whole-window folding score reached an AUROC of 0.289, a motif-scan maximum reached 0.476, and a genomic variant of the head reached 0.391; G4Hunter's maximum 25-nt window reached 0.685 and G4mismatch reached 0.854. Against random genomic negatives (n = 5,793) the same heads reached 0.657, 0.687 and 0.607, where G4Hunter reached 0.966 and pqsfinder 0.945. A model trained on short synthetic oligonucleotides does not transfer to 124-nt genomic windows, and presenting one of these heads as a genome-wide scanner would have been an error of the kind this workbench exists to prevent.

Release v0.5.1 therefore refuses those heads on genome-length input and supplies a scanner trained for the task (Table 5). On chromosomes 2, 8 and 17, held out entirely from training, the scanner reached an AUROC of 0.942 against random negatives and 0.928 against sequence-matched unobserved negatives. The comparison that matters is the second: G4Hunter reached 0.672, pqsfinder 0.536, DeepG4 0.541 and a canonical G3 regex count 0.532, because removing the G-richness contrast removes most of what a motif score measures. Two convolutional models trained on G4-seq performed well but asymmetrically — G4mismatch reached 0.965 and 0.849 on the two negative sets, and G4detector's released models reached 0.970 and 0.818 for the random-negative variant but 0.441 and 0.928 for the PQ-negative variant, each strong on the split matching its own training negatives. Taking the lower of a tool's two scores as a summary of its robustness to the choice of negative, QuadCond's scanner is the highest at 0.928, ahead of G4mismatch at 0.849 and G4detector at 0.818.

The same ordering holds across species. On mouse G4-seq windows the scanner reached 0.928 and 0.775 on the two negative sets, against 0.930 and 0.554 for G4Hunter, 0.942 and 0.636 for G4mismatch, and 0.956 and 0.623 for G4detector's random-negative model. The scanner's lower score, 0.775, is again the highest of the group, and the cross-species drop from 0.928 to 0.775 is reported rather than smoothed: this is a human-trained model tested on mouse.

### Variant triage on promoter positive controls

Six promoter G4s with published structures and biology — MYC Pu27, KIT c-kit1 and c-kit2, VEGFA Pu22, BCL2 Pu39 and hTERT hT21 — were scanned in full and joined to Atlas variant impact scores, together with a length-matched chr22 window containing no quadruplex motif. Across the six loci 633 substitutions joined to a regulatory record, every one of which received a quadrant assignment: 65 high on both axes, 177 structural only, 94 regulatory only and 297 neither. Within the motifs themselves, 462 substitutions were evaluated, of which 135 abolished the motif; maximum absolute predicted ΔTm ranged from 5.58 °C at KIT c-kit1 to 9.38 °C at VEGFA Pu22. The chr22 control window produced 183 substitutions and no structural delta at all, which is the expected behaviour and the reason it was included.

The motif-loss ranking change introduced in v0.5.1 is visible here: 111 of the 633 joined variants abolish the motif and therefore have no ΔTm, and under the previous behaviour they carried no structural rank. Two of the six top-ranked variants across the panel are motif-destroying, so the effect of the change is not marginal.

### Genome-wide, structural disruption and predicted regulatory effect are independent

Across chromosomes 1–22 and X, 150 sampled motifs per chromosome yielded 109,125 Atlas records with no failed queries, joined to 2,707,107 scored substitutions. Table 6 summarises the classes.

The primary within-locus test is significant in the direction opposite to concordance. Motif-destroying substitutions carried a *lower* absolute impact score than substitutions in the flanking window of the same locus: for G4, a median paired difference of −0.0084 with motif-destroying substitutions higher at only 37.1% of 1,195 loci (p = 5.9 × 10⁻²⁰); for i-motifs, −0.0093 at 36.2% of 2,254 loci (p = 3.9 × 10⁻³⁴). The same comparison using every motif substitution rather than only the destroying ones gives 37.2% and 37.3%, which is indistinguishable. The deficit is therefore a property of the motif interval, not of motif destruction.

Sequence composition accounts for it. Motif-destroying positions are CpG-depleted relative to their own flanks — 3.3% against 5.9% for G4 and 3.1% against 5.2% for i-motifs — because a guanine tract contains no CpG dinucleotide, and CpG variants carry disproportionate predicted regulatory weight. In the adjusted model, with the flanking class as reference, the G4 motif-loss coefficient is −0.59 percentile points (95% CI −1.43 to +0.25, p = 0.17) and the i-motif coefficient −0.97 points (−1.54 to −0.39, p = 9.6 × 10⁻⁴), about one percentile point on 301,000 substitutions. Rank correlation between predicted destabilisation magnitude and impact score was 0.018 for G4 and 0.013 for i-motifs.

The distant GC-matched control behaves as the design predicted it would, and its failure is worth recording. It sits significantly *above* the flanking control (+4.42 percentile points, p = 3.6 × 10⁻⁷ for G4), at a mean GC of 0.703 against 0.769 and a CpG fraction of 11.4% against 5.9%, because a high-GC motif that cannot find a GC-matched partner within 2–20 kb contributes no control at all and the control set is drawn from the easier loci. Any analysis using the distant control alone measures that selection effect.

The readout itself is not flat. Measured on the flanking substitutions alone and on the same percentile scale, a substitution at a CpG dinucleotide shifts the absolute impact score by +9.67 percentile points for G4 loci and +10.11 for i-motif loci (Cliff's δ 0.19 and 0.20; median absolute impact 0.074 against 0.043), and substitution type separates strongly (joint p = 7 × 10⁻⁸⁶ and 1 × 10⁻¹⁸³). CpG context and substitution type together explain 0.69% and 0.65% of the rank; structural class explains 0.068% and 0.022%, and adds 0.056% and 0.026% over composition alone. The feature known to matter moves this readout about seven times further than the largest motif effect the data are compatible with, and about sixteen times the point estimate. The bound is therefore an equivalence statement rather than an absence of power: destroying a G4 shifts the predicted-impact percentile by at most 1.43 points in either direction, and an i-motif by at most 1.54.

### Software verification and reproducibility

The client type system separates refused responses from numerical responses, so a refused genomic score or locus posterior is a compile-time error to read rather than a runtime failure. A release verdict requires a non-synthetic service and the expected software version, model version, artifact hash and dataset fingerprint; a pending failure record is written before browser or server startup, so an early failure cannot leave an older verdict standing. Browser checks run against an explicitly synthetic service exercise the integration contract and are marked ineligible for release sign-off. The benchmark, the scanner training and evaluation, the positive-control run and the genome-wide analysis each ship as scripts with recorded inputs and outputs, and the manuscript figures and tables regenerate from frozen metadata without fitting a model.

## Discussion

QuadCond and AENNA-3D put questions that are normally spread across motif finders, stability predictors and structural viewers into one place, and keep the target and provenance of each answer attached to it. The benchmark reported here establishes what that buys and what it does not.

Where a shared target exists, the condition-aware heads lead: melting temperature under grouped validation, i-motif folding and transitional pH, and above all the response of predicted stability to buffer, which no sequence-only tool can express and which the retrained comparators express only partially. The single-substitution result is the one with the most immediate practical consequence. Variant-annotation pipelines for quadruplexes currently rest on differences in a buffer-blind motif score, and against 1,020 measured pairs that difference is at chance for sign and near zero for rank, where a condition-aware ΔTm recovers the sign in four pairs out of five.

Two features of the benchmark deserve emphasis beyond this tool. Released models evaluated on their own training collections reach figures that no honest out-of-sample protocol reproduces — 0.835 against 0.478 for one architecture on the same data — and the difference between grouped and random splitting moved QuadCond's own head by 0.15 R². Neither observation is a criticism of the tools concerned; both are reasons why cross-tool figures quoted from separate papers cannot be compared, and why we retrained rather than quoted.

The genomic result cuts the other way and is reported as such. The folding heads do not transfer from oligonucleotides to genomic windows, and against sequence-matched negatives they are worse than a motif score. The response was to refuse those heads on genomic input and train a scanner for the task, which does hold up on held-out chromosomes and, with a drop, across species. We note that the two convolutional G4-seq models each perform well on the negative set resembling their own training negatives and considerably worse on the other, which suggests that single-split genomic AUROCs are easy to over-read, and we would encourage reporting both.

### What the genome-wide null does and does not establish

Joining the structural axis to a regulatory one across 2.7 million substitutions returns no concordance, with a tight bound. The claim this licenses is about the regulatory model and not about biology: AlphaGenome variant impact carries no signal about G4 or i-motif structural disruption beyond what local sequence composition already explains. This design cannot separate the possibility that motif loss has no regulatory consequence from the possibility that a model trained on sequence-to-function data has no representation of a non-canonical secondary structure, and the second is at least as plausible, since nothing in that training objective requires one.

That distinction is not a weakness of the triage; it is the assumption the triage encodes. The workbench ranks a variant by the lower of its two percentile ranks precisely because the axes are treated as independent lines of evidence rather than as corroboration, and a near-zero incremental contribution of structural class to predicted impact is direct evidence for that independence. The corollary is a negative recommendation with some force: a predicted regulatory impact score should not be used as a proxy for structure-mediated regulatory effect, in either direction.

Testing the biology requires a readout that is not another sequence model. Fine-mapped expression quantitative trait loci, massively parallel reporter measurements over quadruplex variants, and allele-specific quadruplex sequencing would each support the same paired design — motif-destroying substitutions against same-locus flanking substitutions — with a measured rather than predicted outcome. We regard that as the natural next experiment.

### Remaining limitations

The condition-dependent i-motif heads remain the strongest numbers in Table 1 and the narrowest in scope, covering two constructs. Shuffled negatives and genomic proxies answer questions about their own constructed labels; a calibrated probability against dinucleotide shuffles is not a folding fraction, and an antibody peak association is not an equilibrium. Mutation differences are validated here against 1,020 measured pairs in rank and sign, which is considerably more than was previously available, but the paired-estimator spread reported beside each difference remains a statement about the ensemble rather than a calibrated interval for that difference. A further limitation is one of scope rather than accuracy: the conditions this model takes as input are bulk solution variables, and molecular binding partners are not among them. A recent report of peptide-directed folding of an RNA i-motif is the clearest illustration -- the peptide raised folding into the pH 5 range without materially changing the melting temperature, so a partner shifted the folding equilibrium along an axis that neither pH nor cation concentration describes [CITATION TO ADD: peptide-directed folding of the RNA i-motif, Chem. Sci., doi:10.1039/D6SC01203E]. No prediction here should be read as applying inside a protein-bound or ligand-bound complex. The ensemble free-energy calibration is weak — an R² of 0.074 with a residual standard deviation of 1.68 kcal/mol against melting temperatures converted under a two-state van 't Hoff model — and is presented only as a compressed, labelled scale. No prospective experimental validation of prioritised variants has been performed; measuring the melting behaviour of a small panel of top-ranked promoter substitutions is the most direct way to close that gap. The genome-wide analysis is a 150-motif-per-chromosome sample rather than the full set, which bounds precision but, given the width of the reported intervals, would not change the conclusion.

## Data and software availability

The source revision is QuadCond and AENNA-3D v0.5.1. The manuscript figures and tables regenerate from the supplied, unchanged model metadata without fitting a model. The principal model SHA-256 is 493640193ad70974d35582c55b4a721095bfcba50edc6827d9636d2f7b7be0c7, the sequence-only ablation model is 402e9db1f5086df6acb28d2b980dff2e2235e93213f2a2096bed0d59c5bfa100, the core atlas is 145711cd33002ecf6436c52318af55c5feb9e46814239bfb9e4b6f548f7302f6 and the genomic scanner is c5aebae3ddc3129143ce8018a3b4d3b6767e177f4cf3145890a5b822d02973e6. Every asset is verified against this manifest before it is opened, and a mismatch is an error rather than a warning. The benchmark scripts, result tables and figures are distributed with the source under `benchmarks/published_tools`, and the genome-wide analysis under `genomewide`. Two further directories support work that postdates this manuscript and carry no results reported here: `benchmarks/external_validation` holds the harness and provenance records for scoring collections published after the model artifact was frozen, and `benchmarks/ingest_audit` measures how much of a candidate dataset overlaps the atlas and the benchmark evaluation sets before any of it is used for training. Public source repository and archived release DOI: [VERIFY AND INSERT]. Public server: [VERIFY AND INSERT HTTPS URL]. QuadCond carries an MIT licence; AENNA-3D licensing and redistribution terms must be confirmed before public release. Original datasets remain subject to their source terms.

## Supplementary data

Supplementary Table S1 describes all prediction heads. Supplementary Table S2 records atlas source counts and evidence labels. Supplementary Table S3 gives the full benchmark result table including bootstrap intervals for every tool and task. The 1,020-pair single-substitution ΔTm set is released as a standalone benchmark file. The accompanying claim audit, reference metadata and reproducible figure script document the evidence behind this draft.

## Acknowledgements

[Name contributors and facilities only with their agreement.]

## Funding

[Provide funder names, grant numbers and funding for the open-access charge; state no external funding only if confirmed.]

## Author contributions

[Assign author-approved CRediT roles. A.D. and H.S. contributed equally; the shared-first-authorship statement above must match the journal's required wording. Do not infer contributions from software metadata.]

## Conflict of interest statement

[Provide the authors' confirmed declaration.]

## Disclosure of AI-assisted preparation

An AI coding and writing assistant was used to develop and revise the software, to run the benchmark and genome-wide analyses under the authors' direction, and to prepare this draft from the resulting outputs and retrieved literature. [Authors must verify the code, citations and scientific interpretation and adapt this disclosure to the journal's policy and their actual use.]

## Tables

### Table 1. Internal regression performance in the frozen model record

{{REGRESSION_TABLE}}

RMSE, root mean squared error; MAE, mean absolute error. Counts refer to training records, not independent measurements. Coverage is the mean from 20 held-out-group residual partitions at nominal 90%. The two condition-dependent i-motif heads cover only C9 and Tel21C. No confidence intervals for these performance estimates are supplied in the frozen sidecars.

### Table 2. User workflows and limits of interpretation

| Workflow | Output to inspect | Limit of interpretation |
| --- | --- | --- |
| Sequence evidence | Strand, training target, domain status and empirical error scale | No duplex competition or molecular occupancy inference |
| 3D selection | Fold archetype at published helical parameters, linked to an evidence card | Idealised model geometry; no predicted atomic structure or sequence-specific loop conformation |
| Condition sweep | Supported response segment and requested buffer grid | No interpolation claim across unsupported points; two-construct scope for condition-dependent i-motif heads |
| Mutation scan | Estimated ΔTm and paired estimator spread; motif-loss variants ranked by basis rather than magnitude | Rank and sign validated on 1,020 measured pairs; the estimator spread is not a calibrated interval for a difference |
| Genome scan | Called regions with window scores and canonical-motif rescue basis | Trained on G4-seq K⁺; an observed-propensity call, not a folding probability |
| Variant triage | Structural and regulatory percentile ranks and the quadrant | Two model outputs on different scales with no calibrated weighting; a triage label, not a mechanism |
| Batch | Per-record results, exclusions and provenance | No silent removal of invalid entries; record counts are not independent evidence counts |
| Genomic context | Association with assay peak labels | No physical folding probability or same-molecule co-occupancy |

### Table 3. Scope of related resources and of this workbench

| Resource | Question it answers | Buffer as input | Strand-resolved G4 and iM | Applicability or refusal reported |
| --- | --- | --- | --- | --- |
| QGRS Mapper [@qgrs] | Where are putative G4 motifs, and how are they scored | No | G4 only | No |
| G4Hunter [@g4hunter] | G-richness and skewness over a window | No | G4 only | No |
| pqsfinder [@pqsfinder] | Imperfection-tolerant search for potential quadruplex sequences | No | G4 only | No |
| G4-iM Grinder [@grinder] | G4 and iM sequence patterns, including higher-order arrangements | No | Both, by pattern | No |
| QuaDB [@quadb] | Putative quadruplex and i-motif sequences for a given identifier | No | Both, by pattern | No |
| DeepG4 [@deepg4] | Which G4 regions are active in a given cell type | No | G4 only | No |
| G4Boost [@g4boost] | G4 folding state and thermodynamic stability from sequence features | No | G4 only | No |
| G4ShapePredictor [@g4shape] | Which folding topology a G4 adopts | No | G4 only | No |
| G4STAB [@g4stab] | G4 melting temperature from sequence and salt | Salt concentration | G4 only | No |
| iM-Seeker [@im_web] | i-motif folding status and folding strength | No | iM only | No |
| G4Killer [@g4killer] | Which mutations abolish G4 sequence potential | No | G4 only | No |
| ONQUADRO [@onquadro] | Which quadruplex structures have been determined experimentally | Recorded per structure | Both, as deposited | Not applicable |
| G4LDB [@g4ldb] | Which ligands bind G4s and i-motifs | Recorded per assay | Both, as reported | Not applicable |
| QuadCond and AENNA-3D | Which strand-resolved estimate the evidence supports under a stated buffer, and where it stops | K⁺, Na⁺, Li⁺/NH₄⁺, Mg²⁺, pH, temperature, crowder, strand concentration | Both, routed per head | Yes, per head, with refusal distinct from extrapolation |

Entries describe the question each resource is designed to answer, taken from its cited description. Table 4 gives the quantitative comparison on the subset of tasks where a shared target exists.

### Table 4. Head-to-head performance on shared tasks

{{BENCHMARK_TABLE}}

All QuadCond figures are grouped out-of-fold. Intervals are 1,000-replicate bootstraps resampled over sequence groups. Rows marked in-sample are released models evaluated on data they were trained on and are upper bounds, not comparable results; for each such architecture the same estimator retrained on QuadCond's grouped folds is reported immediately below. Buffer-blind tools score identically across buffers by construction, which is recorded as zero condition response rather than omitted. The full table with every tool and both metrics is Supplementary Table S3.

### Table 5. Genomic window scanning against G4-seq K⁺, by negative set

{{SCANNER_TABLE}}

Held-out chromosomes are 2, 8 and 17, excluded entirely from scanner training; the mouse set is a cross-species test of a human-trained model. Random negatives are length-matched genomic windows; PQS-matched negatives are unobserved windows matched on predicted quadruplex-sequence content, which removes the G-richness contrast and is the harder control. The final column takes the lower of a tool's two AUROCs as a summary of how much its performance depends on which negative set is chosen.

### Table 6. Genome-wide substitution classes and predicted regulatory effect

{{GENOMEWIDE_TABLE}}

Percentile shifts are coefficients from the adjusted rank model with the flanking class as reference, controlling for substitution type, CpG context, GC content and chromosome, with standard errors clustered by motif. A positive shift means a higher absolute Atlas variant impact score than same-locus flanking substitutions. The distant control row is reported to document its imbalance, not as a comparison class.

## Figure legends

### Figure 1. Evidence-aware workflow in QuadCond and AENNA-3D

The input sequence and requested conditions are routed to strand-specific prediction heads. Applicability determines whether a numerical estimate is displayed, withheld as extrapolation or refused; genome-length input is refused by the folding heads and routed to the scanner. Biophysical, genomic proxy and auxiliary evidence retain separate interpretations, and comparison exports preserve their context. The 3D panel builds an idealised model of the fold at published helical parameters and does not predict atomic coordinates. This is an original schematic of the implementation, not an experimental result.

{{FIGURE_1}}

### Figure 2. Recorded contribution of condition features to four regression tasks

R² values are read directly from the frozen main-model and sequence-only sidecars. G4 melting temperature and sequence-oriented i-motif transitional pH use sequence-grouped evaluation; the two condition-dependent i-motif tasks use condition-grouped evaluation for C9 and Tel21C only. Bars show point estimates and not uncertainty intervals. The comparison is an internal ablation, not an external benchmark.

{{FIGURE_2}}

### Figure 3. G4 folding discrimination against dinucleotide-preserving shuffles

Receiver operating characteristic curves for QuadCond's grouped out-of-fold folding score and for every comparator run on the identical 5,214-sequence set. Negatives preserve the dinucleotide composition of their own positives, so scores that measure G-richness alone cannot separate the classes; this is a property of the task rather than of those tools, and the panel is included because it is the task the folding head was trained on.

{{FIGURE_3}}

### Figure 4. Topology, melting temperature and i-motif performance

Balanced accuracy for three-class G4 topology, R² for buffer-resolved G4 melting temperature and AUROC for i-motif folding, with bootstrap intervals. Released models evaluated on their own training collections are drawn in a separate style and labelled in-sample; the retrained-on-grouped-folds figure for the same architecture is shown beside each.

{{FIGURE_4}}

### Figure 5. Single-substitution effects on melting temperature

Predicted against measured ΔTm for 1,020 pairs whose wild-type and mutant melting temperatures were measured in the same buffer, for QuadCond and for each variant score in current use. The inset reports sign accuracy restricted to pairs differing by at least 2 °C. The ΔG4Hunter panel is the score underlying current quadruplex variant-annotation practice.

{{FIGURE_5}}

### Figure 6. Condition response and i-motif transitional pH

Left, per-sequence Spearman correlation between predicted and measured melting temperature across buffers, for the 128 sequences measured in more than one condition; buffer-blind tools are constant by construction and sit at zero. Right, predicted against measured i-motif transitional pH for QuadCond and for a retrained iM-Seeker architecture on the 137 sequences for which iM-Seeker returns a motif call.

{{FIGURE_6}}

## References

{{REFERENCES}}
