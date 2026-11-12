#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_hsclrg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/home/suzuki/slurmlogs/%x_%j.log
#
# 05_slurm_hsc_lrg_astrofs4.sbatch — igmdust/HSC wave 1b: LRG foregrounds
# (z 0.40-1.10) x HSC quiescent crayons + star null.  Sized to fill the
# headroom astrofs4 has beside two 10-core jobs (2 cores / 16G keeps
# 2 cores free for the srun login wrapper — the 2026-11-06 rule).
# 2026-11-08 OOM rerun: raw-memmap loader + --chunk 64 (job 83 was
# OOM-killed at 14G: 2000-fg query chunks against 8.1M stars / dense
# crayons held hundreds of millions of pair indices).
#
# I/O: reuses the HSC catalogs that job igmdust_hsc stages to an
# astrofs4-local disk (waits for them, or stages itself after 4 h),
# and rsyncs only the small LRG clustering files.  Nothing here reads
# astrofs2 during pair counting.
#
# Submit by TYPING in your astrofs4 shell (never via ssh command):
#   sbatch ~/github/projects_desi/programs_igmdust/hsc/05_slurm_hsc_lrg_astrofs4.sbatch
set -u
PY=/home/suzuki/miniconda3/bin/python3
PROG=/home/suzuki/github/projects_desi/programs_igmdust/hsc
HSC_SRC=/home/suzuki/data/HSC/catalog
DR1_SRC=/data/astrofs2_1/suzuki/data/DESI/catalog/dr1
echo "== hsc wave1b (LRG fg): host=$(hostname) $(date)"

# ---- same stage-volume pick as 04 (finds the one 04 already chose) --
DEST=""
for cand in /data/astrofs4_2/suzuki /data/astrofs4_1/suzuki \
            /data/astrofs5_1/suzuki; do
  mkdir -p "$cand/igmdust_stage" 2>/dev/null || continue
  avail=$(df -Pm "$cand" 2>/dev/null | awk 'NR==2{print $4}')
  [ "${avail:-0}" -gt 60000 ] || continue
  DEST="$cand/igmdust_stage"; break
done
[ -n "$DEST" ] || { echo "no stage volume with >60G free"; exit 1; }
echo "== using stage $DEST"
mkdir -p "$DEST/HSC" "$DEST/dr1"

# ---- wait for job igmdust_hsc's HSC staging (sizes match source) ----
staged() {
  local n s d
  for n in hsc_dr4_igmdust_669956 hsc_dr4_igmdust_669957 \
           hsc_dr4_igmdust_669937; do
    s=$(stat -c%s "$HSC_SRC/$n" 2>/dev/null) || return 1
    d=$(stat -c%s "$DEST/HSC/$n" 2>/dev/null) || return 1
    [ "$s" = "$d" ] || return 1
  done
}
t0=$(date +%s)
until staged; do
  if [ $(( $(date +%s) - t0 )) -ge 14400 ]; then
    echo "stage wait timed out after 4 h — staging myself $(date)"
    rsync -a "$HSC_SRC/" "$DEST/HSC/" || exit 1
    break
  fi
  echo "waiting for HSC stage... $(date)"; sleep 300
done
echo "== HSC catalogs staged $(date)"

# ---- LRG clustering files (small, one-time NFS read) ----------------
rsync -a "$DR1_SRC"/LRG_NGC_clustering.dat.fits \
         "$DR1_SRC"/LRG_SGC_clustering.dat.fits \
         "$DR1_SRC"/LRG_*_[01]_clustering.ran.fits "$DEST/dr1/" || exit 1
echo "== LRG files staged $(date)"

cd "$PROG" || exit 1
# ---- star null first (cheaper), then crayons, both on 2 cores -------
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
  $PY 02_hsc_paircount.py --fg lrg --mode stars --chunk 64 \
  --hsc-dir "$DEST/HSC" --base "$DEST/dr1" \
  > ~/slurmlogs/igmdust_hls_astrofs4.log 2>&1
echo "lrg stars rc=$? $(date)"
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
  $PY 02_hsc_paircount.py --fg lrg --mode crayons --chunk 64 \
  --hsc-dir "$DEST/HSC" --base "$DEST/dr1" \
  > ~/slurmlogs/igmdust_hlr_astrofs4.log 2>&1
echo "lrg crayons rc=$? $(date)"
echo "== hsc wave1b done $(date)"
