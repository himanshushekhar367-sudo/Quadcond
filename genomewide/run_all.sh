#!/usr/bin/env bash
# Genome-wide QuadCond x AlphaGenome analysis, end to end.
#   bash run_all.sh            # all chromosomes  (long: see README for timings)
#   CHROMS="chr22" bash run_all.sh   # one chromosome, for a pilot
set -euo pipefail

: "${FASTA:?set FASTA to an indexed GRCh38 FASTA}"
: "${AVI:?set AVI to the AlphaGenome AVI tabix file (or a pattern with {chrom}); for a sampled pilot without the bundle use gw_03_avi_live.py, see README}"
export GW_OUT="${GW_OUT:-$PWD/gw_out}"
CHROMS="${CHROMS:-$(echo chr{1..22} chrX)}"
WORKERS="${WORKERS:-$(( $(nproc) - 1 ))}"
HERE="$(cd "$(dirname "$0")" && pwd)"

# Fail before step 1 rather than after step 2: steps 1-2 are the long ones and
# there is no point spending them if the regulatory axis is not on disk.
first_avi="$(printf %s "$AVI" | sed "s/{chrom}/$(printf %s "$CHROMS" | awk '{print $1}')/")"
if [ ! -f "$first_avi" ]; then
  echo "AVI bundle not found: $first_avi" >&2
  echo "Download it from https://deepmind.google.com/science/alphagenome/downloads," >&2
  echo "or run a sampled pilot through the live API instead:" >&2
  echo "  python $HERE/gw_01_motifs.py --fasta \"\$FASTA\" --chroms \$CHROMS" >&2
  echo "  python $HERE/gw_02_structural.py --chroms \$CHROMS --workers $WORKERS" >&2
  echo "  python $HERE/gw_03_avi_live.py --chroms \$CHROMS --max-motifs 500" >&2
  echo "  python $HERE/gw_04_stats.py --fasta \"\$FASTA\" --chroms \$CHROMS" >&2
  exit 1
fi

echo "== 1/4 motifs + matched controls"
python "$HERE/gw_01_motifs.py" --fasta "$FASTA" --chroms $CHROMS

echo "== 2/4 QuadCond structural effect of every SNV in every motif"
python "$HERE/gw_02_structural.py" --chroms $CHROMS --workers "$WORKERS" ${SAMPLE:+--sample $SAMPLE}

echo "== 3/4 AlphaGenome AVI scores for motif and control SNVs"
python "$HERE/gw_03_avi.py" --avi "$AVI" --chroms $CHROMS

echo "== 4/4 statistics"
python "$HERE/gw_04_stats.py" --fasta "$FASTA" --chroms $CHROMS

echo "done -> $GW_OUT/stats/summary.json"
