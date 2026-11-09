#!/bin/bash
# astro4_offload5.sh — Legacy fix wave: propagate the three validated DES
# ingredients to the Legacy program (ms.tex sec. Legacy diagnostics).
#   A) LRG wide (10 kpc - 45 Mpc) CSFD rerun -> 04 (anchor 15 Mpc) -> cov-ML
#   B) BGS wide CSFD rerun in 6 parallel shards -> merge -> 04 -> cov-ML,
#      plus a star-null-subtracted ML solve if the null profile has landed.
# Run ON astro4:  nohup bash astro4_offload5.sh > ~/astro4_offload5.log 2>&1 &
set -u
PY=/home/suzuki/miniconda3/bin/python3
EXT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext
QCAT=$EXT/qso_merged_dr11.fits
PC=/home/suzuki/data/igmdust/paircount
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4   # 7 procs + snl/lcr share 64
mkdir -p ~/slurmlogs
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1

echo "== offload5 start $(date)"

# ---- A: LRG wide CSFD (cheap: small theta_max at z > 0.4) ----------------
(
  $PY 02_paircount.py --fg lrg --qso-cat "$QCAT" --mwcorr \
      --rp-min 0.01 --rp-max 45 --chunk 10000 \
      > ~/slurmlogs/igmdust_lqw_astro4.log 2>&1
  echo "lrg wide rc=$? $(date)"
  $PY 04_dustmag.py --fg lrg_mqcsfd --anchor 15 \
      > ~/slurmlogs/igmdust_lqwfit.log 2>&1
  echo "lrg 04 rc=$? $(date)"
  $PY 10_ml_legacy.py --npz "$PC/paircount_lrg_mqcsfd_qso.npz" \
      > ~/slurmlogs/igmdust_lqwml.log 2>&1
  echo "lrg ML rc=$? $(date)"
) &

# ---- B: BGS wide CSFD, 6 interleaved shards ------------------------------
for i in 0 1 2 3 4 5; do
  $PY 02_paircount.py --fg bgs --qso-cat "$QCAT" --mwcorr \
      --rp-min 0.01 --rp-max 45 --chunk 1000 --nshard 6 --shard $i \
      > ~/slurmlogs/igmdust_bqw_s$i.log 2>&1 &
done
wait   # shards and branch A
$PY 02b_merge_shards.py --name bgs_mqcsfd --nshard 6 \
    > ~/slurmlogs/igmdust_bqwmerge.log 2>&1
echo "bgs merge rc=$? $(date)"
$PY 04_dustmag.py --fg bgs_mqcsfd --anchor 15 \
    > ~/slurmlogs/igmdust_bqwfit.log 2>&1
echo "bgs 04 rc=$? $(date)"
$PY 10_ml_legacy.py --npz "$PC/paircount_bgs_mqcsfd_qso.npz" \
    > ~/slurmlogs/igmdust_bqwml.log 2>&1
echo "bgs ML rc=$? $(date)"

# star-null-subtracted ML, once the designed offset arbiter has landed
NULLNPZ=$PC/paircount_bgsstars_qso.npz
if [ -s "$NULLNPZ" ]; then
  $PY 10_ml_legacy.py --npz "$PC/paircount_bgs_mqcsfd_qso.npz" \
      --null "$NULLNPZ" > ~/slurmlogs/igmdust_bqwmlnull.log 2>&1
  echo "bgs ML(null) rc=$? $(date)"
else
  echo "star null npz not present yet; ML(null) deferred"
fi
echo "== offload5 done $(date)"
