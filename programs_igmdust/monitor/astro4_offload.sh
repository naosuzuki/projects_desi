#!/bin/bash
# astro4_offload.sh — run CPU-bound igmdust jobs on astro4 (which is idle)
# instead of waiting for astrofs4 cores.  Run ON astro4 under nohup:
#   nohup bash astro4_offload.sh > ~/astro4_offload.log 2>&1 &
#
# Stage 1 (now):  CSFD iteration-2 crayon run + McCleary assessment.
#                 No dependencies; inputs (DES crayons, CSFD maps) already
#                 on disk.  Logs go to ~/slurmlogs with the igmdust_rmc /
#                 igmdust_ass names so igm_status.py picks them up and the
#                 verdict appears on the igmprogress web page automatically.
#                 SLURM jobs 75/76 stay queued on astrofs4 as a cross-check
#                 (do NOT scancel 75 — job 77's after:75 release depends
#                 on it starting).
# Stage 2 (auto): as soon as job 60 writes qso_merged_dr11.fits, run the
#                 DES Menard-recipe quasar job here (then job 69 on
#                 astrofs4 is redundant; its rerun is harmless).
set -u
PY=/home/suzuki/miniconda3/bin/python3
DES=/home/suzuki/github/projects_desi/programs_igmdust/des
NC=$(( $(nproc) - 2 )); [ $NC -lt 4 ] && NC=4   # leave cores for updater
export OMP_NUM_THREADS=$NC OPENBLAS_NUM_THREADS=$NC
mkdir -p ~/slurmlogs
# outputs are cwd-relative, same convention as 03_slurm_paircount.sbatch
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1

echo "== astro4 offload start $(date), $NC threads"

# ---- stage 1: CSFD crayons + verdict --------------------------------
$PY "$DES/01_crayon_paircount.py" \
    --fg redmagic --crayon-zmin 0.65 --rp-min 0.01 --rp-max 45 --mwcorr \
    > ~/slurmlogs/igmdust_rmc_astro4.log 2>&1
RC=$?
echo "CSFD run finished rc=$RC $(date)"
if [ $RC -eq 0 ]; then
  $PY "$DES/07_assess_mccleary.py" \
      > ~/slurmlogs/igmdust_ass_astro4.log 2>&1
  echo "assessment finished rc=$? $(date)"
  grep -h MCCLEARY_REPRODUCTION ~/slurmlogs/igmdust_ass_astro4.log || true
fi

# ---- stage 2: DES quasars once job 60's output exists ----------------
QCAT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/qso_merged_dr11.fits
for i in $(seq 1 120); do          # wait up to 10 h
  if [ -s "$QCAT" ]; then
    # job 60 writes the file at the very end; settle 3 min to be safe
    sleep 180
    echo "QCAT present, starting menarddes $(date)"
    $PY "$DES/05_qso_menard_des.py" \
        > ~/slurmlogs/igmdust_mqd_astro4.log 2>&1
    echo "menarddes finished rc=$? $(date)"
    break
  fi
  sleep 300
done
echo "== astro4 offload done $(date)"
