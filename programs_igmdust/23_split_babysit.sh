#!/bin/bash
# 23_split_babysit.sh — BGS brightness-split wave (the 125 kpc
# contamination discriminator, see 22_dip125_check.py and ms.tex
# sec:systematics): bgsbright / bgsfaint halves x merged QSOs with the
# full lw config, 3 shards each, memory- and load-gated, one shard at a
# time (lw_babysit.sh v5 pattern).  If the dip is bright-foreground
# light contamination it should track the bright half.
#   on astro4:  nohup bash 23_split_babysit.sh > ~/split_babysit.log 2>&1 &
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
L=~/slurmlogs
P=~/data/igmdust/paircount
mkdir -p "$L"

one () {  # split shard
  local sp=$1 sh=$2 ns=3
  local npz=$P/paircount_bgs${sp}_mqcsfdlw_s${ns}.${sh}_qso.npz
  local log=$L/igmdust_bsplit_${sp}_s${sh}.log
  [ -s "$npz" ] && { echo "SKIP $sp s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg bgs --fg-split $sp .*--shard $sh( |\$)" \
      >/dev/null && { echo "SKIP $sp s$sh: running"; return; }
  local i=0
  while :; do                                   # RAM + CPU gates (v5)
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 45 ] && [ "$LD" -lt 48 ] && break
    i=$((i+1)); [ $i -gt 180 ] && { echo "GIVEUP $sp s$sh: RAM/load"; return; }
    sleep 60
  done
  echo "LAUNCH $sp s$sh $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 02_paircount.py --fg bgs --fg-split $sp \
    --qso-cat $QSO --mwcorr --qweight --chunk 2000 \
    --rp-min 0.01 --rp-max 45 --nshard $ns --shard $sh > "$log" 2>&1 &
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -q "shard $sh/$ns: fg" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK $sp s$sh $(date '+%T')" \
                           || echo "DIED $sp s$sh (pass will retry)"
}

for pass in 1 2 3 4; do
  echo "== pass $pass $(date '+%F %T')"
  for sp in bright faint; do
    for s in 0 1 2; do one $sp $s; done
  done
done
echo "SPLIT_BABYSIT_DONE $(date)"
# afterwards: 02b_merge_shards.py --name bgsbright_mqcsfdlw --nshard 3
#             02b_merge_shards.py --name bgsfaint_mqcsfdlw  --nshard 3
#             10_ml_legacy.py --npz .../paircount_bgs{bright,faint}_mqcsfdlw_qso.npz
