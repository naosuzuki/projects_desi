#!/bin/bash
# 25_massbin_babysit.sh — mass-evolution wave, ASTRO4 FALLBACK: BGS
# foregrounds in 3 equal weighted-count absolute-W1-magnitude
# (luminosity ~ stellar mass) bins over the full 0.05 < z < 0.45
# range (02_paircount --nmassbin 3 --massbin i), lw config,
# memory/load-gated, one job at a time.  BGS only: the LRG clustering
# catalog carries no fluxes, and LRGs span < 0.5 dex in M* anyway.
# PREFERRED PATH: 26_massbin_fs4.sbatch on astrofs4, which runs in
# parallel with astro4's waves (shared home: both write the same npz,
# and the SKIP-if-npz-exists guard keeps them from colliding only if
# the two are never launched together — run ONE of 25/26, not both).
# This fallback WAITS for the z-bin wave (24_) by requiring all nine
# z-bin npz on disk, so the astro4 waves never contend.
#   on astro4:  nohup bash 25_massbin_babysit.sh > ~/slurmlogs/massbin_babysit.log 2>&1 &
# N.B. 3 bins is the DR1 S/N sweet spot; with DR2 statistics rerun
# with NB=5 below.
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
L=~/slurmlogs
P=~/data/igmdust/paircount
NB=3
mkdir -p "$L"

# ---- hold until the z-bin wave is done (all 9 z-bin npz on disk) ----
# never gate on a launcher's done marker: it appears when the
# launcher's passes end, while its jobs can still be running
i=0
while :; do
  n=$(ls $P/paircount_bgsz{0,1,2}_mqcsfdlw_s2.{0,1}_qso.npz \
         $P/paircount_lrgz{0,1,2}_mqcsfdlw_qso.npz 2>/dev/null | wc -l)
  [ "$n" -ge 9 ] && break
  i=$((i+1)); [ $i -gt 4320 ] && { echo "GIVEUP: z-bin wave never finished"; exit 1; }
  sleep 300
done
echo "ZBIN WAVE DONE -> mass-bin wave $(date '+%F %T')"

one () {  # massbin shard
  local mb=$1 sh=$2 ns=3
  local npz=$P/paircount_bgsm${mb}_mqcsfdlw_s${ns}.${sh}_qso.npz
  local log=$L/igmdust_massbin_m${mb}_s${sh}.log
  [ -s "$npz" ] && { echo "SKIP m$mb s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg bgs .*--nmassbin $NB --massbin $mb .*--shard $sh( |\$)" \
      >/dev/null && { echo "SKIP m$mb s$sh: running"; return; }
  local k=0
  while :; do                                   # RAM + CPU gates (v5)
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 45 ] && [ "$LD" -lt 48 ] && break
    k=$((k+1)); [ $k -gt 180 ] && { echo "GIVEUP m$mb s$sh: RAM/load"; return; }
    sleep 60
  done
  echo "LAUNCH m$mb s$sh $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 02_paircount.py --fg bgs --qso-cat $QSO \
    --mwcorr --qweight --chunk 2000 --rp-min 0.01 --rp-max 45 \
    --nmassbin $NB --massbin $mb --nshard $ns --shard $sh > "$log" 2>&1 &
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -q "shard $sh/$ns: fg" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK m$mb s$sh $(date '+%T')" \
                           || echo "DIED m$mb s$sh (pass will retry)"
}

for pass in 1 2 3 4; do
  echo "== pass $pass $(date '+%F %T')"
  for mb in 0 1 2; do
    for s in 0 1 2; do one $mb $s; done
  done
done
echo "MASSBIN_BABYSIT_DONE $(date)"
# afterwards, per bin: 02b_merge_shards.py --name bgsm<i>_mqcsfdlw --nshard 3
# then 10_ml_legacy.py per bin, then fit the Menard-template amplitude
# A(M) per bin; bin 0 = most massive.  The three bins share the parent
# n(z) only approximately (randoms follow each bin), so quote each
# bin's median z alongside A(M).
