#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_desbin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=18G
#SBATCH --array=0-3
#SBATCH --output=/home/suzuki/slurmlogs/%x_%A_%a.log
#
# 29_desbin_fs4.sbatch — DES binned-foreground studies ON ASTROFS4:
# the DESI BGS x DES Y3 crayon measurement (01_crayon_paircount, MOF
# + CSFD, 45 Mpc) repeated in 2 equal weighted-count z bins and 2
# absolute-W1 luminosity bins (DES-overlap statistics are ~2-3 sigma
# full-sample, so 2 bins, not 3 — Nao 2026-11-15).  No SFR split for
# DES (dropped; redMaGiC-side samples are passive by construction).
#
# Sizing: 18G x 2 cores BACKFILLS into the 4 cores / ~23G that the
# two running igmdust_massbin tasks (10c/18G each) leave free, so
# these run alongside the mass wave instead of queuing behind it,
# while still keeping 2 cores for the login wrapper (10+10+2 = 22 of
# 24; 18*3 = 54G of ~59.4G allocatable — the <20G astrofs4 rule).
#
# Submit by TYPING in your astrofs4 shell (never via ssh command):
#   sbatch ~/github/projects_desi/programs_igmdust/29_desbin_fs4.sbatch
#
# afterwards: 08_ml_estimator-style cov-ML per npz (des/08 or
# 10_ml_legacy with the DECam responses), then A(z) / A(M) amplitudes
# against the Legacy bins at each bin's median z.
set -u
PY=/home/suzuki/miniconda3/bin/python3
PROG=/home/suzuki/github/projects_desi/programs_igmdust/des
P=/home/suzuki/data/igmdust/paircount
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
case $SLURM_ARRAY_TASK_ID in
  0) ARGS="--nzbin 2 --zbin 0";    NPZ=descrayonz0mofcsfd ;;
  1) ARGS="--nzbin 2 --zbin 1";    NPZ=descrayonz1mofcsfd ;;
  2) ARGS="--nmassbin 2 --massbin 0"; NPZ=descrayonm0mofcsfd ;;
  3) ARGS="--nmassbin 2 --massbin 1"; NPZ=descrayonm1mofcsfd ;;
esac
npz=$P/paircount_${NPZ}_qso.npz
[ -s "$npz" ] && { echo "SKIP task $SLURM_ARRAY_TASK_ID: npz exists"; exit 0; }
cd "$PROG" || exit 1
echo "DES bin task $SLURM_ARRAY_TASK_ID ($ARGS) -> $npz"
exec $PY 01_crayon_paircount.py --fg bgs --mof --mwcorr \
  --rp-min 0.01 --rp-max 45 $ARGS
