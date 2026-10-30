#!/bin/zsh
# Runs detached after run_walk_v2.sh (2026-10-29): when the shared v2 walk has stored every coadd's native spectra
# (all 4 walk logs end in 'done:'), DESI v2 stacks from the native store (3 shards), then merge (SDSS v2 + DESI v2)
# -> quasar/v2/qso_stacks.npz and EM-PCA v2 (global + one basis per normalization window). Log: after_walk_v2.log.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar
log() { echo "$(date '+%F %T') $*"; }
empca_all() {   # $1 stacks file, $2 tag
  EMPCA_SERIES_LAST=0 OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --tag $2
  EMPCA_SERIES_LAST=0 OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 0 0.8 --tag ${2}_z0008
  EMPCA_SERIES_LAST=0 OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 0.8 2.0 --tag ${2}_z0820
  OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 2.0 9.9 --tag ${2}_z2099
}
until [ $(grep -l '^done:' $S/walk_v2_?.log 2>/dev/null | wc -l) -ge 4 ] && ! pgrep -f '02_passQ_desi.py --all --with-native' >/dev/null; do sleep 120; done
log "v2 walk done: DESI v2 stacks from the native store"
for i in 0 1 2; do
  OPENBLAS_NUM_THREADS=1 nice -n 10 python3 05_stack.py --v2 --survey desi --shard $i/3 >> $S/qso_stack_v2_desi$i.log 2>&1 &
done
wait
until [ $(grep -l '^done:' $S/qso_stack_v2_sdss?.log 2>/dev/null | wc -l) -ge 3 ]; do sleep 60; done
log "merge v2 + EM-PCA v2"
python3 05_stack.py --v2 --merge && empca_all $P/v2/qso_stacks.npz v2
log "v2 finished"
