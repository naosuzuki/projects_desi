#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_pg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=/home/suzuki/slurmlogs/%x_%j.log
#
# 13_slurm_menard.sbatch — igmdust: the Menard-replica photometric-
# foreground run on astrofs4 (power-budget fix, 2026-11-08: BGS is 17x
# sparser than Menard+10's photometric galaxies; this run uses every
# r<21 DR11 galaxy as foreground -> ~17x the pair counts inside 500 kpc).
#
# Submit ONCE with 12_submit_menard_astrofs4.sh (typed in your astrofs4
# shell; never via ssh remote command).  Jobs: scan -> {count0, count1}
# -> finish, chained with SLURM dependencies.  Each job keeps to the
# 10-of-24-cores / 20-of-61-GB safety rule; two run side by side, so the
# count phase uses 20 cores while logins stay possible.
set -u
PY=/home/suzuki/miniconda3/bin/python3
EXT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext
RAND=/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11/randoms
QCAT=$EXT/qso_merged_dr11.fits
PC=/home/suzuki/data/igmdust/paircount
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
echo "== menard-replica job: arg=$1 host=$(hostname) $(date)"

wait_randoms() {    # 00d on astro4 downloads + cuts; marker lands last
  until [ -f "$RAND/.downloaded" ] && [ -s "$EXT/photran_dr11.fits" ]; do
    echo "  waiting for photran_dr11.fits (00d on astro4)... $(date)"
    sleep 300
  done
}

count_shard() {
  $PY 02_paircount.py --fg photgal --qso-cat "$QCAT" --mwcorr \
      --rp-min 0.01 --rp-max 45 --chunk 2000 --nshard 4 --shard "$1" \
      > ~/slurmlogs/igmdust_pgw_s$1.log 2>&1
  echo "shard $1 rc=$? $(date)"
}

case "$1" in
  scan)
    exec $PY 11_sweepgal.py --mode galaxies --nproc 10 ;;
  count0)
    wait_randoms; count_shard 0; count_shard 1 ;;
  count1)
    wait_randoms; count_shard 2; count_shard 3 ;;
  finish)
    $PY 02b_merge_shards.py --name photgal_mqcsfd --nshard 4 \
        > ~/slurmlogs/igmdust_pgmerge.log 2>&1
    echo "merge rc=$? $(date)"
    $PY 04_dustmag.py --fg photgal_mqcsfd --anchor 15 \
        > ~/slurmlogs/igmdust_pgfit.log 2>&1
    echo "04 rc=$? $(date)"
    $PY 10_ml_legacy.py --npz "$PC/paircount_photgal_mqcsfd_qso.npz" \
        > ~/slurmlogs/igmdust_pgml.log 2>&1
    echo "ML rc=$? $(date)"
    NULLNPZ=$PC/paircount_bgsstars_qso.npz
    if [ -s "$NULLNPZ" ]; then
      $PY 10_ml_legacy.py --npz "$PC/paircount_photgal_mqcsfd_qso.npz" \
          --null "$NULLNPZ" > ~/slurmlogs/igmdust_pgmlnull.log 2>&1
      echo "ML(null) rc=$? $(date)"
    else
      echo "star null npz not present yet; ML(null) deferred"
    fi ;;
  *) echo "unknown arg: $1 (scan|count0|count1|finish)"; exit 1 ;;
esac
