#!/bin/bash
# 03_run_astro4.sh — igmdust/HSC: first measurement wave on astro4.
#   A) star null: 8.1M i<22.5 point sources x BGS (systematics arbiter)
#   B) quiescent crayons x BGS in 2 shards -> merge
# Run ON astro4:
#   ( nohup bash ~/github/projects_desi/programs_igmdust/hsc/03_run_astro4.sh \
#       > ~/hsc_wave1.log 2>&1 & )
# Cov-ML re-solve with HSC R coefficients follows once 10_ml grows a
# --rcoef option (jk accumulations are saved by 02 regardless).
set -u
PY=/home/suzuki/miniconda3/bin/python3
PROG=/home/suzuki/github/projects_desi/programs_igmdust/hsc
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
mkdir -p ~/slurmlogs
cd "$PROG" || exit 1
echo "== hsc wave1 start $(date)"

$PY 02_hsc_paircount.py --mode stars \
    > ~/slurmlogs/igmdust_hsn_astro4.log 2>&1 &
for i in 0 1; do
  $PY 02_hsc_paircount.py --mode crayons --nshard 2 --shard $i \
      > ~/slurmlogs/igmdust_hcr_s$i.log 2>&1 &
done
wait
echo "pair counts done $(date)"
$PY ~/github/projects_desi/programs_igmdust/02b_merge_shards.py \
    --name hscbgscsfd --nshard 2 > ~/slurmlogs/igmdust_hcrmerge.log 2>&1
echo "merge rc=$? $(date)"
echo "== hsc wave1 done $(date)"
