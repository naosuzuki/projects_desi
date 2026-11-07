#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=/home/suzuki/slurmlogs/%x_%j.log
#
# 03_slurm_paircount.sbatch — igmdust: run 02_paircount.py as a SLURM batch
# job on the cluster (submitted FROM astrofs4, the submit host).
#
# Usage (type this yourself in an interactive astrofs4 shell):
#   sbatch ~/github/projects_desi/programs_igmdust/03_slurm_paircount.sbatch selftest
#   sbatch ~/github/projects_desi/programs_igmdust/03_slurm_paircount.sbatch lrg
#   sbatch ~/github/projects_desi/programs_igmdust/03_slurm_paircount.sbatch bgs
#   squeue -u suzuki          # watch;   logs in ~/slurmlogs/
#
# SAFETY (learned 2026-11-06 the hard way): on astrofs4 every login is itself
# a SLURM job, so a batch job must NEVER request the whole node.  This script
# asks for 10 of 24 cores and 20 of 61 GB, leaving plenty for logins and a
# second job.  Two of these can run side by side.
#   - Do NOT raise --cpus-per-task above ~10 or --mem above ~24G.
#   - Submit by TYPING sbatch in your astrofs4 shell; never via ssh remote
#     command (the login wrapper duplicates those).

set -u
PY=/home/suzuki/miniconda3/bin/python3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}

ARG=${1:-selftest}
echo "== igmdust job: arg=$ARG host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?} $(date)"
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1

# environment / mount sanity (printed to the job log)
$PY -c "import numpy, scipy, astropy; print('env OK:', numpy.__version__, scipy.__version__, astropy.__version__)" || exit 1
ls -d /data/astrofs2_1/suzuki/data/DESI/catalog/dr1 >/dev/null 2>&1 \
  && echo "mount OK: /data/astrofs2_1" \
  || echo "WARN: /data/astrofs2_1 not visible on $(hostname) (selftest still fine)"

case "$ARG" in
  selftest)  exec $PY 02_paircount.py --selftest ;;
  bgs|lrg)   exec $PY 02_paircount.py --fg "$ARG" ;;
  descrayon) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/01_crayon_paircount.py ;;
  speczcalib) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/02_speczcalib.py ;;
  qsodes) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/03_qso_desoverlap.py ;;
  lrgcrayon) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/01_crayon_paircount.py \
    --fg lrg ;;
  rmgcrayon) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/01_crayon_paircount.py \
    --fg redmagic --crayon-zmin 0.65 --rp-min 0.01 --rp-max 45 ;;
  menarddes) exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/05_qso_menard_des.py ;;
  sweepmatch) exec $PY 06_sweepmatch.py ;;
  redmapper) shift; exec $PY \
    /home/suzuki/github/projects_des/programs_igmdust/04_redmapper_crayons.py \
    "$@" ;;
  bgsmerged)
    QCAT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
    $PY 02_paircount.py --fg bgs --qso-cat "$QCAT" --chunk 10000 || exit 1
    exec $PY 04_dustmag.py --fg bgs_mq ;;
  *) echo "unknown arg: $ARG (selftest|bgs|lrg|descrayon|speczcalib|qsodes)"
     exit 1 ;;
esac
