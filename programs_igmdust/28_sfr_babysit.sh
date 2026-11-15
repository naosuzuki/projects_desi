#!/bin/bash
# 28_sfr_babysit.sh — star-formation split wave on astro4: BGS sf
# (blue cloud) / passive (red sequence) halves from the z-conditional
# weighted-median g-r split (02_paircount --fg-csplit), x merged QSOs,
# full lw config, 3 shards each, memory- and load-gated, one shard at
# a time (v5 pattern).  Both halves share the parent randoms: the
# per-z median split preserves each half's n(z) by construction.
# Launches IMMEDIATELY (no hold): it may coexist with the z-bin wave,
# so the RAM gate is raised to 50G and --chunk 1000 shrinks the
# transient pair-index gathers (the bright-s0 OOM was a 2.7 GB gather
# at chunk 2000 with six jobs peaking together).
#   on astro4:  nohup bash 28_sfr_babysit.sh > ~/sfr_babysit.log 2>&1 &
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
L=~/slurmlogs
P=~/data/igmdust/paircount
mkdir -p "$L"

one () {  # csplit shard
  local cs=$1 sh=$2 ns=3
  local npz=$P/paircount_bgs${cs}_mqcsfdlw_s${ns}.${sh}_qso.npz
  local log=$L/igmdust_sfr_${cs}_s${sh}.log
  [ -s "$npz" ] && { echo "SKIP $cs s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg bgs .*--fg-csplit $cs .*--shard $sh( |\$)" \
      >/dev/null && { echo "SKIP $cs s$sh: running"; return; }
  local k=0
  while :; do                                   # RAM + CPU gates
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 50 ] && [ "$LD" -lt 48 ] && break
    k=$((k+1)); [ $k -gt 360 ] && { echo "GIVEUP $cs s$sh: RAM/load"; return; }
    sleep 60
  done
  echo "LAUNCH $cs s$sh $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 02_paircount.py --fg bgs --fg-csplit $cs \
    --qso-cat $QSO --mwcorr --qweight --chunk 1000 \
    --rp-min 0.01 --rp-max 45 --nshard $ns --shard $sh > "$log" 2>&1 &
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -q "shard $sh/$ns: fg" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK $cs s$sh $(date '+%T')" \
                           || echo "DIED $cs s$sh (pass will retry)"
}

for pass in 1 2 3 4; do
  echo "== pass $pass $(date '+%F %T')"
  for cs in passive sf; do
    for s in 0 1 2; do one $cs $s; done
  done
done
echo "SFR_BABYSIT_DONE $(date)"
# afterwards: 02b_merge_shards.py --name bgssf_mqcsfdlw --nshard 3
#             02b_merge_shards.py --name bgspassive_mqcsfdlw --nshard 3
#             10_ml_legacy.py per half; compare A_V profiles --
#             Menard+10 found reddening around both blue and red
#             galaxies; the passive half is also the cross-check for
#             the 125 kpc light-contamination bin (red galaxies are
#             brighter at fixed z).
