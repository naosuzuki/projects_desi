#!/bin/bash
# 24_zbin_babysit.sh — redshift-evolution wave: BGS and LRG foregrounds
# in 3 equal weighted-count z bins each (02_paircount --nzbin 3 --zbin i),
# lw config, memory/load-gated, one job at a time.  Six A(z) points over
# 0.05 < z_fg < 1.1.  WAITS for the brightness-split wave (23_) to
# finish before launching anything, so the two waves never contend.
#   on astro4:  nohup bash 24_zbin_babysit.sh > ~/zbin_babysit.log 2>&1 &
# N.B. 3 bins is the DR1 S/N sweet spot; with DR2 statistics rerun with
# NB=5 below.
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
L=~/slurmlogs
P=~/data/igmdust/paircount
NB=3
mkdir -p "$L"

# ---- hold until the split wave is done (marker, or all 6 split npz) ----
i=0
while :; do
  grep -q SPLIT_BABYSIT_DONE ~/split_babysit.log 2>/dev/null && break
  n=$(ls $P/paircount_bgs{bright,faint}_mqcsfdlw_s3.*_qso.npz 2>/dev/null | wc -l)
  [ "$n" -ge 6 ] && break
  i=$((i+1)); [ $i -gt 2880 ] && { echo "GIVEUP: split wave never finished"; exit 1; }
  sleep 300
done
echo "SPLIT WAVE DONE -> z-bin wave $(date '+%F %T')"

one () {  # fg zbin nshard shard
  local fg=$1 zb=$2 ns=$3 sh=$4
  local tag; [ "$ns" = 1 ] && tag="" || tag="_s${ns}.${sh}"
  local npz=$P/paircount_${fg}z${zb}_mqcsfdlw${tag}_qso.npz
  local log=$L/igmdust_zbin_${fg}z${zb}_s${sh}.log
  [ -s "$npz" ] && { echo "SKIP $fg z$zb s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg $fg .*--nzbin $NB --zbin $zb .*--shard $sh( |\$)" \
      >/dev/null && { echo "SKIP $fg z$zb s$sh: running"; return; }
  local k=0
  while :; do                                   # RAM + CPU gates (v5)
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 45 ] && [ "$LD" -lt 48 ] && break
    k=$((k+1)); [ $k -gt 180 ] && { echo "GIVEUP $fg z$zb s$sh: RAM/load"; return; }
    sleep 60
  done
  echo "LAUNCH $fg z$zb s$sh $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 02_paircount.py --fg $fg --qso-cat $QSO \
    --mwcorr --qweight --chunk 2000 --rp-min 0.01 --rp-max 45 \
    --nzbin $NB --zbin $zb --nshard $ns --shard $sh > "$log" 2>&1 &
  local pid=$! j=0
  local peak="shard $sh/$ns: fg"            # printed only when sharded
  [ "$ns" = 1 ] && peak="jackknife regions:"
  while kill -0 $pid 2>/dev/null && ! grep -q "$peak" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK $fg z$zb s$sh $(date '+%T')" \
                           || echo "DIED $fg z$zb s$sh (pass will retry)"
}

for pass in 1 2 3 4; do
  echo "== pass $pass $(date '+%F %T')"
  for zb in 0 1 2; do
    one bgs $zb 2 0
    one bgs $zb 2 1
  done
  for zb in 0 1 2; do one lrg $zb 1 0; done
done
echo "ZBIN_BABYSIT_DONE $(date)"
# afterwards, per bin: 02b_merge_shards.py --name bgsz<i>_mqcsfdlw --nshard 2
# (lrg bins are unsharded), then 10_ml_legacy.py --npz ... per bin, then
# fit the Menard-template amplitude A(z) per bin; convert each bin with
# the extinction-curve factor at lambda_obs/(1+z_fg) before interpreting
# evolution.
