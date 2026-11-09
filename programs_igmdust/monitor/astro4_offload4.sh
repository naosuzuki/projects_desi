#!/bin/bash
# astro4_offload4.sh — McCleary iteration 4: covariance-weighted ML
# estimator (the final ladder rung).  Reruns the MOF+CSFD pair count with
# per-region accumulations saved, re-solves E with the full color-color
# jackknife covariance (08), and re-assesses (07) -> verdict on the page.
#   nohup bash astro4_offload4.sh > ~/astro4_offload4.log 2>&1 &
set -u
PY=/home/suzuki/miniconda3/bin/python3
DES=/home/suzuki/github/projects_des/programs_igmdust
PCDIR=/home/suzuki/data/igmdust/paircount
mkdir -p ~/slurmlogs
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
export OMP_NUM_THREADS=20 OPENBLAS_NUM_THREADS=20

echo "== iteration 4 (ML estimator) start $(date)"
$PY "$DES/01_crayon_paircount.py" \
    --fg redmagic --crayon-zmin 0.65 --rp-min 0.01 --rp-max 45 \
    --mof --mwcorr \
    > ~/slurmlogs/igmdust_rmc_mof.log 2>&1
RC=$?
echo "MOF+CSFD rerun (jk arrays) rc=$RC $(date)"
[ $RC -ne 0 ] && exit 1

$PY "$DES/08_ml_estimator.py" \
    --npz "$PCDIR/paircount_rmgmofcsfd_qso.npz" \
    > ~/slurmlogs/igmdust_ml.log 2>&1
RC=$?
echo "ML solve rc=$RC $(date)"
[ $RC -ne 0 ] && { tail -3 ~/slurmlogs/igmdust_ml.log; exit 1; }

$PY "$DES/07_assess_mccleary.py" \
    --npz "$PCDIR/paircount_rmgmofcsfdml_qso.npz" \
    > ~/slurmlogs/igmdust_ass_ml.log 2>&1
echo "assessment rc=$? $(date)"
grep -h MCCLEARY_REPRODUCTION ~/slurmlogs/igmdust_ass_ml.log || true
echo "== iteration 4 done $(date)"
