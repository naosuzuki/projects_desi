#!/bin/bash
# astro4_offload2.sh — second wave of astro4 delegation so everything
# finishes before morning.  Division of labor (no job runs on both
# machines — the SLURM duplicates 71/74 are scancel'ed):
#   astro4   : CSFD verdict (offload 1), menarddes (offload 1, gated),
#              lrgmerged   (here, gated on job 60's QSO catalog),
#              legacycrayon(here, gated on job 72's lrg_dr11.fits)
#   astrofs4 : job 60 + 72 scans (NFS-bound, must stay), then SLURM 61
#              (bgsmerged, afterok:60) and 73 (star null, afterok:72)
# Run ON astro4:  nohup bash astro4_offload2.sh > ~/astro4_offload2.log 2>&1 &
set -u
PY=/home/suzuki/miniconda3/bin/python3
EXT=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext
QCAT=$EXT/qso_merged_dr11.fits
LRGCAT=$EXT/lrg_dr11.fits
export OMP_NUM_THREADS=20 OPENBLAS_NUM_THREADS=20   # 2 jobs + offload1 share 64
mkdir -p ~/slurmlogs
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1

wait_for() {  # wait_for <file>; file is written atomically at job end,
  for i in $(seq 1 150); do        # up to 12.5 h
    [ -s "$1" ] && { sleep 180; return 0; }   # settle time after appearing
    sleep 300
  done
  echo "TIMEOUT waiting for $1"; return 1
}

( wait_for "$QCAT" && {
    echo "lrgmerged start $(date)"
    $PY 02_paircount.py --fg lrg --qso-cat "$QCAT" --chunk 10000 \
        > ~/slurmlogs/igmdust_lmq_astro4.log 2>&1 \
      && $PY 04_dustmag.py --fg lrg_mq \
        >> ~/slurmlogs/igmdust_lmq_astro4.log 2>&1
    echo "lrgmerged done rc=$? $(date)"
  } ) &

( wait_for "$LRGCAT" && {
    echo "legacycrayon start $(date)"
    $PY 07_crayon_legacy.py > ~/slurmlogs/igmdust_lcr_astro4.log 2>&1
    echo "legacycrayon done rc=$? $(date)"
  } ) &

wait
echo "== astro4 offload2 done $(date)"
