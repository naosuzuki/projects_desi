#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_hsc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=18G
#SBATCH --output=/home/suzuki/slurmlogs/%x_%j.log
#
# 04_slurm_hsc_astrofs4.sbatch — igmdust/HSC wave 1 on astrofs4, with
# deliberate I/O separation (Nao, 2026-11-08): the inputs are STAGED to
# a disk served by a different machine than astrofs2 (which the sweep
# scan is hammering), preferring astrofs4's own local disks so the HSC
# reads never touch the network at all.
#
# Submit by TYPING in your astrofs4 shell (never via ssh command):
#   sbatch ~/github/projects_desi/programs_igmdust/hsc/04_slurm_hsc_astrofs4.sbatch
#
# 10 cores / 18 GB: allocatable memory is < RealMemory (24G and even
# 22G pended with "Resources" at 59.4G total); 18G puts the 3-job
# total at 55,296 MB, the watermark proven to fit on 2026-11-08, and
# coexists with replica job 79 (16+20+18G) once the scan ends.  With
# the Menard-replica chain also queued, SLURM
# runs two 10-core jobs at a time, keeping 20 of 24 cores busy and
# logins always possible (the 2026-11-06 lockout rule).
# 2026-11-08 OOM rerun: raw-memmap loader + --chunk 64 (job 82 was
# OOM-killed: astropy allocated 43 GB slicing the 18 GB crayon file,
# and 2000-fg query chunks held ~250M pair indices).
set -u
PY=/home/suzuki/miniconda3/bin/python3
PROG=/home/suzuki/github/projects_desi/programs_igmdust/hsc
HSC_SRC=/home/suzuki/data/HSC/catalog
DR1_SRC=/data/astrofs2_1/suzuki/data/DESI/catalog/dr1
echo "== hsc wave1: host=$(hostname) $(date)"

# ---- pick a stage volume NOT served by astrofs2 ---------------------
DEST=""
for cand in /data/astrofs4_2/suzuki /data/astrofs4_1/suzuki \
            /data/astrofs5_1/suzuki; do
  mkdir -p "$cand/igmdust_stage" 2>/dev/null || continue
  avail=$(df -Pm "$cand" 2>/dev/null | awk 'NR==2{print $4}')
  [ "${avail:-0}" -gt 60000 ] || continue
  DEST="$cand/igmdust_stage"; break
done
[ -n "$DEST" ] || { echo "no stage volume with >60G free"; exit 1; }
echo "== staging to $DEST ($(df -Ph $DEST | awk 'NR==2{print $4}') free)"

mkdir -p "$DEST/HSC" "$DEST/dr1"
rsync -a --info=progress2 "$HSC_SRC/" "$DEST/HSC/" || exit 1
rsync -a "$DR1_SRC"/BGS_ANY_*_clustering.dat.fits \
         "$DR1_SRC"/BGS_ANY_*_clustering.ran.fits "$DEST/dr1/" || exit 1
echo "== staged $(du -sh $DEST | cut -f1) $(date)"

cd "$PROG" || exit 1
# ---- stars null first, ALONE (its pair buffers need the room) ------
OMP_NUM_THREADS=10 OPENBLAS_NUM_THREADS=10 \
  $PY 02_hsc_paircount.py --mode stars --chunk 64 --hsc-dir "$DEST/HSC" \
  --base "$DEST/dr1" > ~/slurmlogs/igmdust_hsn_astrofs4.log 2>&1
echo "stars rc=$? $(date)"
# ---- crayon shards 0+1 in parallel (5+5 of 10 cores) ----------------
OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5 \
  $PY 02_hsc_paircount.py --mode crayons --nshard 2 --shard 0 --chunk 64 \
  --hsc-dir "$DEST/HSC" --base "$DEST/dr1" \
  > ~/slurmlogs/igmdust_hcr_s0.log 2>&1 &
OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5 \
  $PY 02_hsc_paircount.py --mode crayons --nshard 2 --shard 1 --chunk 64 \
  --hsc-dir "$DEST/HSC" --base "$DEST/dr1" \
  > ~/slurmlogs/igmdust_hcr_s1.log 2>&1 &
wait
echo "shards done $(date)"
$PY /home/suzuki/github/projects_desi/programs_igmdust/02b_merge_shards.py \
    --name hscbgscsfd --nshard 2 > ~/slurmlogs/igmdust_hcrmerge.log 2>&1
echo "merge rc=$? $(date)"
echo "== hsc wave1 done $(date)"
