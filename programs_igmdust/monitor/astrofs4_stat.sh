#!/bin/bash
# astrofs4_stat.sh — resident stat collector for the igmprogress host cards.
# Runs ON astrofs4 under nohup; writes CPU load / memory / SLURM-allocated
# cores into the shared home every 2 min, where astro4's igm_updater.sh
# merges it into status.json.  No network, no credentials.
#   nohup bash astrofs4_stat.sh > ~/astrofs4_stat.log 2>&1 &
set -u
OUT=~/.igmstats
mkdir -p "$OUT"
for i in $(seq 1 1000); do                      # ~33 h, relaunchable
  { cat /proc/loadavg; nproc; free -b | awk 'NR==2{print $2, $7}'; \
    squeue -h -o %C -t R 2>/dev/null | awk '{s+=$1} END {print "ALLOC", s+0}'; \
  } > "$OUT/astrofs4.stat.tmp" 2>/dev/null
  mv "$OUT/astrofs4.stat.tmp" "$OUT/astrofs4.stat"
  sleep 120
done
echo "ASTROFS4_STAT_ENDED $(date)"
