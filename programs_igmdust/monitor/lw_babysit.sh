#!/bin/bash
# lw_babysit.sh — relaunch missing lw-recount jobs one at a time,
# memory-gated (2026-11-11). A job is skipped if its npz exists or a
# matching process runs; otherwise wait for >=45 GB available, launch,
# then wait until its log shows the post-peak line (sharding done /
# regions built) or the process dies. 3 outer passes catch re-kills.
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
L=~/slurmlogs
P=~/data/igmdust/paircount

one () {  # fg nshard shard log npz pastpeak_regex
  local fg=$1 ns=$2 sh=$3 log=$4 npz=$5 peak=$6
  # BGS MUST run chunked: unchunked, its low-z 45 Mpc apertures build
  # ~80M-pair temp arrays per fg batch (118 GB RSS seen 2026-11-12).
  local X=""; [ "$fg" = bgs ] && X="--chunk 2000"
  [ -f "$npz" ] && { echo "SKIP $fg s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg $fg .*--nshard $ns --shard $sh\$" \
      >/dev/null && { echo "SKIP $fg s$sh: running"; return; }
  if [ "$ns" = 1 ]; then
    pgrep -f "02_paircount.py --fg $fg --qso-cat.*--rp-max 45\$" \
        >/dev/null && { echo "SKIP $fg: running"; return; }
  fi
  local i=0
  while :; do
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 45 ] && [ "$LD" -lt 48 ] && break   # RAM + CPU gates
    i=$((i+1)); [ $i -gt 120 ] && { echo "GIVEUP $fg s$sh: RAM/load"; return; }
    sleep 60
  done
  echo "LAUNCH $fg s$sh $(date '+%T') avail=${A}G load=${LD}"
  if [ "$ns" = 1 ]; then
    nohup nice -n 19 python3 02_paircount.py --fg $fg --qso-cat $QSO \
      --mwcorr --qweight --rp-min 0.01 --rp-max 45 > $log 2>&1 &
  else
    nohup nice -n 19 python3 02_paircount.py --fg $fg --qso-cat $QSO \
      --mwcorr --qweight --rp-min 0.01 --rp-max 45 $X \
      --nshard $ns --shard $sh > $log 2>&1 &
  fi
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -qE "$peak" $log; do
    j=$((j+1)); [ $j -gt 120 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK $fg s$sh $(date '+%T')" \
                           || echo "DIED $fg s$sh (pass will retry)"
}

for pass in 1 2 3; do
  echo "== pass $pass $(date '+%F %T')"
  for s in 0 1 2 3; do
    one photgal 4 $s $L/igmdust_pgw_lw_s$s.log \
        $P/paircount_photgal_mqcsfdlw_s4.${s}_qso.npz "shard $s/4: fg"
  done
  # BGS is back here (2026-11-12): astrofs4's 20G cgroup cap OOM-killed
  # the load phase twice (arrays 90 at 12G, 97 at 20G) — the pre-shard
  # full-catalog phase needs more than 20G no matter the shard count.
  for s in 0 1 2 3 4 5; do
    one bgs 6 $s $L/igmdust_bqw_lw_s$s.log \
        $P/paircount_bgs_mqcsfdlw_s6.${s}_qso.npz "shard $s/6: fg"
  done
  one lrg 1 0 $L/igmdust_lqw_lw.log \
      $P/paircount_lrg_mqcsfdlw_qso.npz "jackknife regions"
  sleep 300
done
echo "BABYSIT_DONE $(date '+%F %T')"
