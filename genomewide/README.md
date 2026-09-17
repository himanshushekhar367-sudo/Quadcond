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
