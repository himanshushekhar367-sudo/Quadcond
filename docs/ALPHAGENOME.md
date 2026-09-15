# The AlphaGenome join

AlphaGenome Atlas answers "does this variant change regulation". QuadCond
answers "does this variant change the local G-quadruplex or i-motif, under this
buffer". Neither is a mechanism on its own. Together they are a testable
hypothesis with a named structure in it, and the useful output is a *ranking
within one element* rather than another pathogenicity score.

## Why the join is exact

`scans.mutation_scan` enumerates every single-base substitution in a window.
The Atlas publishes a record for every single-base substitution in the genome.
Given the window's coordinates, those two enumerations are the same set, so each
structural row pairs with at most one regulatory record on
`(chromosome, position, reference, alternate)`.

There is no overlap heuristic here, no nearest-feature rule and no window
arithmetic to be trusted — which matters, because the failure mode of an
approximate join is a table where every number still looks reasonable. The
coordinate is checked rather than assumed: an off-by-one is not a degraded
result, it is a different variant, and the tests pin exactly that case.

## Getting the regulatory axis

Three sources, tried in this order.

**1. The published AVI Tabix bundle — prefer this.** Google DeepMind publishes
AVI SNV scores as an 88.5 GB Tabix bundle **licensed for commercial and
non-commercial use**, alongside non-commercial splicing (20.6 GB) and SHAP
feature-importance (283.9 GB) bundles, from
<https://deepmind.google.com/science/alphagenome/downloads>. Tabix is why the
size does not matter: the index turns "every SNV in this 30-nt element" into a
seek and a few kilobytes read.

```bash
export ALPHAGENOME_AVI_TABIX=/data/alphagenome/avi_snv_scores.tsv.gz   # plus its .tbi
quadcond variant GGGTTAGGGTTAGGGTTAGGG --chromosome chr8 --start 128748301
```

Column positions are read from the file's own header, never hardcoded — a fixed
column order is a silent mis-join the first time the publisher adds a field, and
a mis-joined score is worse than none. Needs either `pysam` or the `tabix`
binary.

**2. A table you exported.** Any CSV/TSV with chromosome, position, ref and alt
columns, in wide form (one column per scorer) or long form (`scorer`/`score`).
For working from a subset, or from a network that cannot reach Google at all.

```bash
quadcond variant … --atlas-table my_region_scores.tsv
```

**3. The live Atlas API.** Needs a key (`ALPHAGENOME_API_KEY`, issued by Google
DeepMind, non-commercial use) and reaches `gdmscience.googleapis.com` over gRPC.
Queried by interval rather than per variant, because one element is ~90
substitutions and 90 round trips to answer a question the service answers in one
is the mistake the mutation scan already had removed.

```bash
pip install alphagenome
quadcond variant … --live
```

With none of the three configured the structural half still runs and the
response says the regulatory axis is absent. **A missing score is not a low
score**, and a variant with no record is reported unclassified rather than
ranked last.

## The 2×2, and what it is not

| | QuadCond low | QuadCond high |
|---|---|---|
| **Atlas low** | neither | structural only |
| **Atlas high** | regulatory only | **both** |

The bottom-right cell is the point: the regulatory model expects an effect *and*
this head expects the local structure to change, which is a candidate
structure-mediated regulatory mechanism and a reason to put the variant on an
instrument.

It is a triage device, not a classifier. Neither threshold was fitted — by
default each axis is split at its 75th percentile *within the request*, and
absolute thresholds can be supplied instead and are recorded in the run record.
Nothing here has been benchmarked against measured regulatory variants.

**There is deliberately no combined score.** The two axes are model output on
different scales and no weighting between them has been calibrated, so summing
or multiplying them would manufacture the kind of number the rest of this
project spends its effort refusing to emit. What is reported instead is
`combined_rank`: the *smaller* of the two within-request percentile ranks, so a
variant ranks high only when it ranks high on both. A percentile is a statement
about the other rows in the same table and is not comparable across requests.

## The validation this needs before it is a claim

Assemble regulatory variants with experimentally characterised G4s or i-motifs
at the same locus, freeze the protocol, then compare three predictors on the
same set: Atlas alone, QuadCond alone, and the join. Report whether the join
recovers experimentally supported variants at a higher rate than either axis,
with an overlap audit against QuadCond's training records and against
AlphaGenome's. Until that exists, `combined_rank` is a prioritisation index and
the quadrant is a triage label.

The structural axis also carries every limitation it carries elsewhere:
mutation-effect prediction has not been validated here as its own task, and the
paired estimator spread is ensemble disagreement rather than a validated
interval on a difference.

## Terms

Atlas predictions are model output, not measurements, and AlphaGenome has not
been validated for and is not approved for clinical use. The AVI static download
is permissively licensed; the splicing, feature-attribution and raw-feature
tiers are non-commercial. QuadCond records the path and header it read and makes
no claim about which release produced a file it was pointed at — the file is the
provenance.

## Operational notes

The service resolves one source at startup and caches it, so two rows of one
study cannot come from two different sources without saying so. The client
cannot name a path or post a key: a prediction endpoint that reads any path a
caller sends is a file-disclosure endpoint with a scientific interface.

`gdmscience.googleapis.com` does not resolve on every network. If the live
backend is unreachable, that is reported as a source failure and never as a
locus with no regulatory effect.
