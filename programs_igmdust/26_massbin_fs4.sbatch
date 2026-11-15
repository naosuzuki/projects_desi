#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_massbin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=45G
#SBATCH --array=0-8
#SBATCH --output=/home/suzuki/slurmlogs/%x_%A_%a.log
#
# 26_massbin_fs4.sbatch — mass-evolution wave ON ASTROFS4: BGS in 3
# equal weighted-count absolute-W1-magnitude (luminosity ~ stellar
# mass) bins over the full 0.05 < z < 0.45 range, lw config, 3 shards
# per bin.  Runs in PARALLEL with astro4's split / z-bin waves: the
# mass wave needs no z-bin outputs, only the shared catalogs, so the
# serialization on astro4 was purely a RAM constraint there.
#
# Sizing: 45G means SLURM runs ONE array task at a time inside
# astrofs4's ~59.4G allocatable (22G jobs pended at 2026-11-08's
# watermark test); 20 cores leaves headroom for the srun login
# wrapper (the 2026-11-06 lockout rule).  --chunk 1000 (astro4 runs
# use 2000) halves the transient pair-index arrays as the OOM margin
# under the hard 45G cap.
#
# Submit by TYPING in your astrofs4 shell (never via ssh command):
#   sbatch ~/github/projects_desi/programs_igmdust/26_massbin_fs4.sbatch
#
# afterwards (anywhere on the shared home):
#   02b_merge_shards.py --name bgsm<i>_mqcsfdlw --nshard 3  per bin,
#   then 10_ml_legacy.py per bin, then fit the Menard-template
#   amplitude A(M) per bin; bin 0 = most massive.  Quote each bin's
#   median z alongside A(M) (the bins share the parent n(z) only
#   approximately; randoms follow each bin).
set -u
PY=/home/suzuki/miniconda3/bin/python3
PROG=/home/suzuki/github/projects_desi/programs_igmdust
QSO=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
P=/home/suzuki/data/igmdust/paircount
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mb=$((SLURM_ARRAY_TASK_ID / 3))
sh=$((SLURM_ARRAY_TASK_ID % 3))
npz=$P/paircount_bgsm${mb}_mqcsfdlw_s3.${sh}_qso.npz
[ -s "$npz" ] && { echo "SKIP m$mb s$sh: npz exists"; exit 0; }
[ -r "$QSO" ] || { echo "FATAL: QSO catalog not readable from astrofs4"; exit 1; }
cd "$PROG" || exit 1
echo "mass bin $mb shard $sh -> $npz"
exec $PY 02_paircount.py --fg bgs --qso-cat $QSO --mwcorr --qweight \
  --chunk 1000 --rp-min 0.01 --rp-max 45 \
  --nmassbin 3 --massbin $mb --nshard 3 --shard $sh
