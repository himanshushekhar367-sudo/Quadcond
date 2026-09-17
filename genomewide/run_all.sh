#!/usr/bin/env bash
# Genome-wide QuadCond x AlphaGenome analysis, end to end.
#   bash run_all.sh            # all chromosomes  (long: see README for timings)
#   CHROMS="chr22" bash run_all.sh   # one chromosome, for a pilot
set -euo pipefail

: "${FASTA:?set FASTA to an indexed GRCh38 FASTA}"
: "${AVI:?set AVI to the AlphaGenome AVI tabix file (or a pattern with {chrom})}"
export GW_OUT="${GW_OUT:-$PWD/gw_out}"
CHROMS="${CHROMS:-$(echo chr{1..22} chrX)}"
WORKERS="${WORKERS:-$(( $(nproc) - 1 ))}"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "== 1/4 motifs + matched controls"
python "$HERE/gw_01_motifs.py" --fasta "$FASTA" --chroms $CHROMS

echo "== 2/4 QuadCond structural effect of every SNV in every motif"
python "$HERE/gw_02_structural.py" --chroms $CHROMS --workers "$WORKERS" ${SAMPLE:+--sample $SAMPLE}

echo "== 3/4 AlphaGenome AVI scores for motif and control SNVs"
python "$HERE/gw_03_avi.py" --avi "$AVI" --chroms $CHROMS

echo "== 4/4 statistics"
python "$HERE/gw_04_stats.py" --fasta "$FASTA" --chroms $CHROMS

echo "done -> $GW_OUT/stats/summary.json"
