# Genome-wide: do structure-changing SNVs carry larger regulatory effects?

This is the biological test the AlphaGenome join was built for, at genome scale
and without scoring all 9 billion variants. Every canonical G4 and i-motif on
chr1-22 and chrX is enumerated; every single-base substitution inside those
motifs is scored by QuadCond; a matched control window is drawn for every motif;
AlphaGenome AVI scores are read for both sets from the published Tabix bundle;
and the comparison is made with the locus held fixed.

Everything below is one bash session (WSL2 Ubuntu on this machine is fine).

## 0. Environment

```bash
sudo apt-get update && sudo apt-get install -y tabix          # htslib, for indexing
conda create -y -n qcgw python=3.11 && conda activate qcgw

cd /mnt/c/Users/pc/Downloads/Quadcond/quadcond-repo
pip install -e . -r requirements-model.txt
pip install pysam pandas pyarrow scipy statsmodels matplotlib

quadcond assets status        # model + atlas_core must say "ok"
```

## 1. Reference genome

```bash
export FASTA=/mnt/c/Users/pc/GRCh38.p13.genome.fa
samtools faidx "$FASTA"       # or: python -c "import pysam;pysam.faidx('$FASTA')"
```

The FASTA must be GRCh38/hg38 and use `chr` names (the GENCODE file above does).
AlphaGenome Atlas is on the same assembly; a mismatch here is a wrong answer,
not a warning.

Observed on this machine: chr22 gives 10,737 G4 and 22,410 i-motif motifs with
21,775 matched controls, and step 2 scored 3.0 M SNVs for it; chr1 gives 99,116
motifs and 8.9 M SNVs. Step 1 and step 2 need no network and no AVI file.

## 2. The AlphaGenome AVI bundle

Download the AVI SNV score bundle (and its `.tbi`) from
<https://deepmind.google.com/science/alphagenome/downloads>, onto a disk with
~90 GB free — a drive mounted under WSL, not the container:

```bash
export AVI=/mnt/d/alphagenome/avi_snv_scores.tsv.gz    # adjust to the real name
ls -l "$AVI" "$AVI".tbi || tabix -s 1 -b 2 -e 2 -S 1 "$AVI"   # only if no index shipped
```

Nothing here assumes the column order: `quadcond.alphagenome` reads the file's
own header and fails loudly if it cannot find chromosome/position/ref/alt.

### No bundle yet? Run the pilot through the live API

`run_all.sh` now stops before step 1 if the bundle is not on disk, because steps
1-2 are the long ones. Without it, use the sampled live-API path — it needs only
`ALPHAGENOME_API_KEY` and produces the same file format, so step 4 is unchanged:

```bash
python genomewide/gw_01_motifs.py --fasta "$FASTA" --chroms chr22
python genomewide/gw_02_structural.py --chroms chr22 --workers $(( $(nproc) - 1 ))
python genomewide/gw_03_avi_live.py --chroms chr22 --max-motifs 500   # 1,000 queries
python genomewide/gw_04_stats.py --fasta "$FASTA" --chroms chr22
```

It asks the Atlas for `AVI_SCORE` alone. The full default scorer set returns
every track for every variant — about 20 MB for a 21-nt window, which overruns
the gRPC client's 4 MB receive cap and fails with `RESOURCE_EXHAUSTED`. The cap
is raised here as well (`--max-message-mb`, default 256).

`gw_03_avi_live.py` samples motifs that actually have a motif-destroying SNV, at
random with a fixed seed, so the pilot is not spent on windows where nothing
happens — and it is a **sample**: report `--max-motifs` and `--strategy` in the
methods. One interval query per motif and per control is ~66,000 calls for chr22
alone and ~2.2 M genome-wide, which is why the bundle is the way to do this
properly.

## 3. Pilot on one chromosome (do this first)

```bash
export GW_OUT=$PWD/gw_out
CHROMS=chr22 bash run_all.sh
cat "$GW_OUT/stats/summary.json"
```

chr22 finishes in well under an hour and produces the same tables as the full
run. Check `stats/summary.json` looks sane before committing to the genome.

## 4. Full run

```bash
bash run_all.sh                       # all chromosomes
# or, to halve the structural step at the cost of power:
SAMPLE=0.5 bash run_all.sh
```

Each step writes one file per chromosome and skips chromosomes already done, so
the run is resumable: interrupt it and start it again.

Rough scale on a desktop CPU: ~1.1 M canonical motifs genome-wide (G4 + i-motif,
both strands), ~40 M motif SNVs, ~800-1,500 SNV predictions per second per core.
That is a long overnight run on 8 cores for step 2; steps 1, 3 and 4 are each
one to three hours. Disk: ~10-20 GB of output.

## What comes out

```
$GW_OUT/motifs/<chrom>.tsv.gz       motifs + their matched control windows
$GW_OUT/structural/<chrom>.tsv.gz   per SNV: state (motif_lost / retained), wild-type value, delta
$GW_OUT/avi/<chrom>.tsv.gz          per SNV: AlphaGenome scores, motif set and control set
$GW_OUT/stats/summary.json          every test below
$GW_OUT/stats/{G4,iM}_avi_by_class.png
```

## The tests, and which one is the result

1. **|AVI| by structural class** (motif_lost / destabilising / moderate /
   neutral / stabilising / control), with Mann-Whitney and Cliff's delta.
2. **Enrichment in the top 1 % of |AVI|**, the threshold taken from the control
   SNVs, as a Fisher odds ratio per class.
3. **Adjusted logistic model**: class + substitution type + CpG context + GC +
   chromosome, standard errors clustered by motif. G>A at CpG is both the
   commonest G4-breaking substitution and a mutational hotspot, so the
   unadjusted contrast in 1-2 is confounded by design.
4. **Within-motif paired test** (`within_motif_lost_vs_retained`): for each
   motif, the mean |AVI| of its motif-destroying SNVs against the mean |AVI| of
   its motif-retaining SNVs, Wilcoxon signed-rank. **This is the headline.**
   Both sets sit in the same regulatory element, so promoter activity, chromatin
   state and conservation are held fixed and only the structural consequence
   differs.
5. **Dose-response**: Spearman between predicted destabilisation and |AVI|
   among motif-retaining SNVs.
6. **Leave-one-chromosome-out**: the odds ratio from test 2, per chromosome. A
   result that only appears on a few chromosomes is not a result.

## What the result can and cannot say

Both axes are model predictions. AlphaGenome's AVI is not a measurement, and
QuadCond's delta is a predicted melting-temperature change benchmarked on
oligonucleotides (`benchmarks/published_tools`: Spearman 0.54 on 1,020 measured
single-substitution pairs, direction correct 81 % at |dTm| >= 2 degC). A positive
result is evidence that two independent models agree about a class of variants,
which is a reason to test those variants experimentally. It is not evidence that
any particular variant changes expression through a G-quadruplex.

The honest negative result is worth as much: if motif-destroying SNVs carry no
more predicted regulatory effect than their motif-retaining neighbours, the 2x2
prioritisation in the workbench has no basis at genome scale, and the paper
should say so.

## Result of the sampled pilot (150 motifs per chromosome, chr1-22 + chrX)

109,125 AVI records from the live Atlas API, 0 failed queries, joined to
2,707,107 scored SNVs. `gw_out/stats/summary.json` and
`gw_out/stats/sensitivity.json` hold the numbers below.

**The primary within-locus test is significant in the direction opposite to the
hypothesis.** Motif-destroying SNVs score *lower* on |AVI| than SNVs in the
flanking window of the same locus: G4 median difference -0.0084, higher in
37.1 % of 1,195 loci, p = 5.9e-20; iM -0.0093, 36.2 % of 2,254 loci,
p = 3.9e-34.

**That is not about motif loss.** `within_locus_motif_vs_flank`, which pools
every motif SNV rather than only the destroying ones, gives 37.2 % and 37.3 %
-- indistinguishable. The whole motif interval sits below its own flanks, and
destroying the motif adds nothing on top of being inside it.

**Sequence composition explains it.** Motif-destroying positions are
CpG-depleted (G4 3.3 %, iM 3.1 %, against 5.9 % and 5.2 % in the flanks),
because a G-run contains no CG dinucleotide. Adjusting for CpG, substitution
type, GC and chromosome, and clustering by motif, the G4 contrast collapses to
-0.0059 percentile [-0.0143, +0.0025], p = 0.17; iM stays at -0.0097
[-0.0154, -0.0039], about one percentile point out of a hundred on 301,000
SNVs. Dose-response between predicted destabilisation and |AVI| is Spearman
0.018 (G4) and 0.013 (iM).

**The readout is not flat.** `gw_05_sensitivity.py`, using flank rows only, on
the same percentile scale: a SNV at a CpG shifts |AVI| by +0.097 percentile
(G4) and +0.101 (iM), Cliff's delta 0.19-0.20, median |AVI| 0.074 against
0.043. Substitution type separates too (joint p = 7e-86 and 1e-183). CpG and
substitution type explain 0.69 % and 0.65 % of the |AVI| rank; motif class
explains 0.068 % and 0.022 %, and adds 0.056 % and 0.026 % over composition
alone. The feature that is known to matter moves the readout roughly sixteen
times further than the largest effect motif destruction is compatible with.

**Equivalence bound.** Destroying a G4 shifts the |AVI| percentile by at most
1.43 points in either direction; for an i-motif, at most 1.54.

**What the distant control was doing.** The GC-matched control 2-20 kb away
sits *above* the flank (+0.044 percentile, p = 3.6e-7 for G4), at GC 0.703
against 0.769 and CpG 11.4 % against 5.9 %. Every number computed against it --
the whole of `gw_out/stats_noflank/` -- is reading that confound. Use the flank.

`per_chromosome_OR_motif_lost_vs_control` contains zeros, infinities and NaNs:
it is a tail statistic on a handful of events per chromosome. It is noise, not
heterogeneity, and should not be reported.

### What this licenses saying, and what it does not

The claim the data support is about the model, not about biology: **AlphaGenome's
variant-effect predictions carry no signal about G4 or i-motif structural
disruption beyond what local sequence composition already explains.** This
design cannot separate "G4 loss has no regulatory consequence" from
"AlphaGenome has no representation of G4 folding", and the second is at least
as likely, since nothing in its training objective asks it to model a
non-canonical secondary structure.

That is still a useful result for the workbench. The 2x2 triage assumes the
structural and regulatory axes are non-redundant; a near-zero incremental R^2
and a 1.4-percentile equivalence bound are direct evidence that they are. The
axes should be read as independent lines of evidence, and AVI should not be
used as a proxy for structure-mediated regulatory effect.

Testing the biology needs a readout that is not another sequence model:
fine-mapped eQTL credible sets, MPRA measurements over G4 variants, or
allele-specific G4 ChIP. Motif-destroying SNVs against same-locus flank SNVs,
same paired design, measured outcome.

### Scope

This is a sample: 150 motifs per chromosome, motif window plus 100 bp flanks,
one AVI scorer (AVI_SCORE). Record `--max-motifs`, `--flank` and `--strategy`
in any methods section. The full genome-wide set needs the AVI Tabix bundle
(`gw_03_avi.py`), not the live API.
