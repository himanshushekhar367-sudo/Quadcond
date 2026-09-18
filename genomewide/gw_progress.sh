#!/usr/bin/env bash
# Progress of a running gw_03_avi_live.py pilot.
#
#   bash genomewide/gw_progress.sh [logfile]
#   watch -n 30 bash genomewide/gw_progress.sh avi_flank.log
#
# Finished chromosomes come from the output directory; the current chromosome
# from the last progress line in the log; the rate from the mtimes of the
# finished files, which is the only timing the run actually records.
set -u
LOG=${1:-avi_flank.log}
OUT="${GW_OUT:-$PWD/gw_out}/avi"
TOTAL_CHR=${TOTAL_CHR:-23}
PER_CHR=${PER_CHR:-150}

mapfile -t files < <(ls -1t "$OUT"/*.tsv.gz 2>/dev/null)
done_chr=${#files[@]}
cur_line=$(grep -oE 'chr[0-9XYxy]+: [0-9]+/[0-9]+ motifs' "$LOG" 2>/dev/null | tail -1)
cur_done=$(sed -E 's#.*: ([0-9]+)/[0-9]+ motifs#\1#' <<<"${cur_line:-0/0 motifs}")
[[ "$cur_done" =~ ^[0-9]+$ ]] || cur_done=0

motifs_done=$(( done_chr * PER_CHR + cur_done ))
motifs_total=$(( TOTAL_CHR * PER_CHR ))
pct=$(awk -v a="$motifs_done" -v b="$motifs_total" 'BEGIN{printf "%.1f", 100*a/b}')

running="no (finished or stopped)"
pgrep -f gw_03_avi_live >/dev/null && running="yes (pid $(pgrep -f gw_03_avi_live | tr '\n' ' '))"

echo "run active     : $running"
echo "chromosomes    : $done_chr/$TOTAL_CHR done"
echo "current        : ${cur_line:-waiting for the first progress line}"
echo "motifs queried : $motifs_done/$motifs_total  (${pct}%)"

if (( done_chr >= 2 )); then
  newest=$(stat -c %Y "${files[0]}")
  oldest=$(stat -c %Y "${files[$((done_chr-1))]}")
  span=$(( newest - oldest ))
  if (( span > 0 )); then
    rate=$(( span / (done_chr - 1) ))
    left=$(( TOTAL_CHR - done_chr ))
    eta=$(( rate * left - (cur_done * rate / PER_CHR) ))
    (( eta < 0 )) && eta=0
    printf "pace           : %dm %ds per chromosome\n" $((rate/60)) $((rate%60))
    printf "estimated left : %dh %dm  (done about %s)\n" $((eta/3600)) $(((eta%3600)/60)) "$(date -d "+${eta} seconds" +%H:%M)"
  fi
fi

if (( done_chr > 0 )); then
  last="${files[0]}"
  echo -n "last file      : $(basename "$last")  "
  zcat "$last" 2>/dev/null | awk -F'\t' 'NR>1{c[$1]++} END{printf "%d motif, %d flank, %d control rows\n", c["motif"], c["flank"], c["control"]}'
fi
echo "recent log     :"; tail -3 "$LOG" 2>/dev/null | sed 's/^/  /'
