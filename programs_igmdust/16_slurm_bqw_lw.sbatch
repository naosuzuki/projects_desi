#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_bqwlw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=12G
#SBATCH --array=0-5
#SBATCH --output=/home/suzuki/slurmlogs/%x_%A_%a.log
#
# 16_slurm_bqw_lw.sbatch — BGS × QSOs wide lw recount (locus residuals
# + per-QSO ivar weights) on astrofs4: 6 shards as one array job, while
# astro4 runs the replica shards + LRG of the same wave (2026-11-11).
# 12G × 3 cores per task: SLURM packs ~3 alongside the HSC job (job 87)
# and queues the rest — the scheduler is the memory/CPU babysitter.
# --chunk 2000 is the proven 20G-cap setting (13_slurm_menard).
#
# Submit ONCE (typed in your astrofs4 shell; never via ssh remote cmd):
#     sbatch ~/github/projects_desi/programs_igmdust/16_slurm_bqw_lw.sbatch
set -u
PY=/home/suzuki/miniconda3/bin/python3
QCAT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-3}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-3}
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
S=${SLURM_ARRAY_TASK_ID}
echo "== bqw lw shard $S host=$(hostname) $(date)"
$PY 02_paircount.py --fg bgs --qso-cat "$QCAT" --mwcorr --qweight \
    --rp-min 0.01 --rp-max 45 --chunk 2000 --nshard 6 --shard "$S" \
    > ~/slurmlogs/igmdust_bqw_lw_s$S.log 2>&1
echo "shard $S rc=$? $(date)"
