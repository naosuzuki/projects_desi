#!/bin/bash
# 33_pathab_babysit.sh — McCleary with DESI Legacy data, both paths,
# on astro4 (no redMaGiC run needed in either):
#
#  PATH A  her exact samples, Legacy photometry.  redMaGiC highdens
#          foregrounds (0.15<z<0.45) x highlum crayons (0.65<z<0.90),
#          the crayon colors re-measured from DR11 (Tractor) sweeps:
#          32_rmg_dr11match.py (sweep match, I/O bound, ~2-5 h), then
#          des/01_crayon_paircount.py --dr11-colors with the flags of
#          the completed MOF reproduction.  Third photometry pipeline
#          on the same galaxies (SOF vs MOF moved the answer x2).
#  PATH B  her foreground SELECTION, DESI spectroscopy.  BGS passive
#          half & 0.15<z<0.45 & M_r<-20.2 (02_paircount --fg-mcc) x
#          merged QSOs, lw config, 3 shards.
#
# The match starts at once in the background; Path B's shards run
# meanwhile (gated, one at a time); Path A's count waits for the match.
#   on astro4:  nohup bash 33_pathab_babysit.sh > ~/slurmlogs/pathab_babysit.log 2>&1 &
cd ~/github/projects_desi/programs_igmdust || exit 1
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
DES=/home/suzuki/data/DES/catalog/redmagic
L=~/slurmlogs
P=~/data/igmdust/paircount
MATCH=$DES/rmg_crayons_dr11.fits
mkdir -p "$L"

gate () {  # label -> returns 1 on give-up
  local k=0
  while :; do
    A=$(free -g | awk '/^Mem:/{print $7}')
    LD=$(awk '{printf "%d", $1}' /proc/loadavg)
    [ "$A" -ge 50 ] && [ "$LD" -lt 48 ] && return 0
    k=$((k+1)); [ $k -gt 360 ] && { echo "GIVEUP $1: RAM/load"; return 1; }
    sleep 60
  done
}

# ---- Path A step 1: sweep match, in the background -------------------
if [ ! -s "$MATCH" ] && ! pgrep -f "32_rmg_dr11match[.]py" >/dev/null; then
  echo "LAUNCH A-match $(date '+%F %T')"
  nohup nice -n 19 python3 32_rmg_dr11match.py \
    > "$L/igmdust_rdm_match.log" 2>&1 &
fi

pathb () {  # shard
  local sh=$1 ns=3
  local npz=$P/paircount_bgsmcc_mqcsfdlw_s${ns}.${sh}_qso.npz
  local log=$L/igmdust_mccb_s${sh}.log
  [ -s "$npz" ] && { echo "SKIP B s$sh: npz exists"; return; }
  pgrep -f "02_paircount.py --fg bgs --fg-mcc .*--shard $sh( |\$)" \
      >/dev/null && { echo "SKIP B s$sh: running"; return; }
  gate "B s$sh" || return
  echo "LAUNCH B s$sh $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 02_paircount.py --fg bgs --fg-mcc \
    --qso-cat $QSO --mwcorr --qweight --chunk 1000 \
    --rp-min 0.01 --rp-max 45 --nshard $ns --shard $sh > "$log" 2>&1 &
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -q "shard $sh/$ns: fg" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK B s$sh $(date '+%T')" \
                           || echo "DIED B s$sh (pass will retry)"
}

patha () {
  local npz=$P/paircount_rmgdr11csfd_qso.npz
  local log=$L/igmdust_rdm_count.log
  [ -s "$npz" ] && { echo "SKIP A-count: npz exists"; return; }
  pgrep -f "01_crayon_paircount.py .*--dr11-colors" >/dev/null \
      && { echo "SKIP A-count: running"; return; }
  # hold for the match: output file on disk AND the matcher gone
  local i=0
  while [ ! -s "$MATCH" ] || pgrep -f "32_rmg_dr11match[.]py" >/dev/null; do
    i=$((i+1)); [ $i -gt 720 ] && { echo "GIVEUP A-count: no match file"; return; }
    sleep 60
  done
  gate "A-count" || return
  echo "LAUNCH A-count $(date '+%T') avail=${A}G load=${LD}"
  nohup nice -n 19 python3 des/01_crayon_paircount.py --fg redmagic \
    --dr11-colors $MATCH --crayon-zmin 0.65 --rp-min 0.01 --rp-max 45 \
    --mwcorr > "$log" 2>&1 &
  local pid=$! j=0
  while kill -0 $pid 2>/dev/null && ! grep -q "jackknife regions:" "$log"; do
    j=$((j+1)); [ $j -gt 150 ] && break
    sleep 60
  done
  kill -0 $pid 2>/dev/null && echo "PASTPEAK A-count $(date '+%T')" \
                           || echo "DIED A-count (pass will retry)"
}

for pass in 1 2 3; do
  echo "== pass $pass $(date '+%F %T')"
  for s in 0 1 2; do pathb $s; done
  patha
done
echo "PATHAB_BABYSIT_PASSES_DONE $(date)  (jobs may still be running)"
# afterwards:
#   B: 02b_merge_shards.py --name bgsmcc_mqcsfdlw --nshard 3, then
#      10_ml_legacy.py --npz .../paircount_bgsmcc_mqcsfdlw_qso.npz
#   A: des/08_ml_estimator.py on paircount_rmgdr11csfd_qso.npz, and
#      compare against paircount_rmgmofcsfdml (the MOF reproduction)
#      and the McCleary law (07_assess_mccleary.py).
