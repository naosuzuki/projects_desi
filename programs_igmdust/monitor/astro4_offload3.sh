#!/bin/bash
# astro4_offload3.sh — McCleary iteration 3: MOF photometry.
# Runs ON astro4 (network + CPU):
#   1. fetch Y3 Gold MOF colors for the 850k crayons (00b --mof, NOIRLab TAP)
#   2. crayon pair counts with MOF colors + CSFD re-correction (= her config)
#   3. automated assessment -> verdict on the igmprogress page
#   nohup bash astro4_offload3.sh > ~/astro4_offload3.log 2>&1 &
set -u
PY=/home/suzuki/miniconda3/bin/python3
DES=/home/suzuki/github/projects_desi/programs_igmdust/des
PCDIR=/home/suzuki/data/igmdust/paircount
mkdir -p ~/slurmlogs
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
export OMP_NUM_THREADS=20 OPENBLAS_NUM_THREADS=20

echo "== iteration 3 (MOF) start $(date)"
cd "$DES" && $PY 00b_fetch_y3gold_colors.py --mof \
  > ~/des_y3gold_mof.log 2>&1
RC=$?
echo "MOF fetch rc=$RC $(date)"
grep -c FAILED ~/des_y3gold_mof.log || true
[ $RC -ne 0 ] && { echo "fetch failed; rerun me (parts are cached)"; exit 1; }

cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
$PY "$DES/01_crayon_paircount.py" \
    --fg redmagic --crayon-zmin 0.65 --rp-min 0.01 --rp-max 45 \
    --mof --mwcorr \
    > ~/slurmlogs/igmdust_rmc_mof.log 2>&1
RC=$?
echo "MOF+CSFD crayon run rc=$RC $(date)"
[ $RC -ne 0 ] && exit 1

$PY "$DES/07_assess_mccleary.py" \
    --npz "$PCDIR/paircount_rmgmofcsfd_qso.npz" \
    > ~/slurmlogs/igmdust_ass_mof.log 2>&1
echo "assessment rc=$? $(date)"
grep -h MCCLEARY_REPRODUCTION ~/slurmlogs/igmdust_ass_mof.log || true
echo "== iteration 3 done $(date)"
