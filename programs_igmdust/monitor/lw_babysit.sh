#!/bin/bash
# lw_babysit.sh — relaunch missing lw-recount jobs one at a time,
# memory-gated (2026-11-11, v5 2026-11-12). A job is skipped if its npz
# exists or a matching process runs; otherwise wait for the RAM gate,
# launch, then wait until its log shows the post-peak line (sharding
# done / regions built) or the process dies. 6 outer passes catch
# re-kills.
# v5: PER-FOREGROUND RAM GATES. A single 45G gate killed the wave twice
# on 2026-11-12 (12:17, 13:31 OOM): chunked BGS shards grow to ~20G
# steady (the 9G reading was early in the run), so three of them hold
# 60G, and a replica (photgal) shard — 40G steady, 30-55G load peak —
# launched into a 56-71G gap dies once everything grows. photgal now
# requires 80G free (peak 55 + margin); it simply waits until the BGS
# queue drains. BGS runs before photgal for the same reason.
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
  local MIN=45; [ "$fg" = photgal ] && MIN=80    # per-fg RAM gate (v5)
  local i=0
  while :; do
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge "$MIN" ] && [ "$LD" -lt 48 ] && break   # RAM + CPU gates
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

for pass in 1 2 3 4 5 6; do
  echo "== pass $pass $(date '+%F %T')"
  # BGS first (v5): 6 x ~20G drain fastest; photgal's 80G gate then
  # opens as they finish. BGS stays off astrofs4: its 20G cgroup cap
  # OOM-killed the load phase twice (arrays 90 at 12G, 97 at 20G).
  for s in 0 1 2 3 4 5; do
    one bgs 6 $s $L/igmdust_bqw_lw_s$s.log \
        $P/paircount_bgs_mqcsfdlw_s6.${s}_qso.npz "shard $s/6: fg"
  done
  for s in 0 1 2 3; do
    one photgal 4 $s $L/igmdust_pgw_lw_s$s.log \
        $P/paircount_photgal_mqcsfdlw_s4.${s}_qso.npz "shard $s/4: fg"
  done
  one lrg 1 0 $L/igmdust_lqw_lw.log \
      $P/paircount_lrg_mqcsfdlw_qso.npz "jackknife regions"
  sleep 300
done
echo "BABYSIT_DONE $(date '+%F %T')"
