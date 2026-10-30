#!/bin/zsh
# Runs detached after run_stack.sh (2026-10-28 01:00): when every DESI shard has finished its dark coadds (marker
# acc_desi_shard*of3.npz.dark_done), merge -> qso_stacks_dark.npz and EM-PCA v1_dark (global + one basis per
# normalization window); when all shards are done, merge -> qso_stacks.npz and EM-PCA v1. Log: after_stack.log.
cd ~/github/projects_desi/programs_quasars
A=~/data/sdssdesi_products/quasar/stack_acc; S=~/data/sdssdesi_status
log() { echo "$(date '+%F %T') $*"; }
empca_all() {   # $1 stacks file, $2 tag
  OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --tag $2
  OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 0 0.8 --tag ${2}_z0008
  OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 0.8 2.0 --tag ${2}_z0820
  OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $1 --zrange 2.0 9.9 --tag ${2}_z2099
}
until [ $(ls $A/acc_desi_shard*of3.npz.dark_done 2>/dev/null | wc -l) -ge 3 ]; do sleep 60; done
log "dark coadds done on all shards: provisional merge + EM-PCA v1_dark"
python3 05_stack.py --merge --out qso_stacks_dark && empca_all ~/data/sdssdesi_products/quasar/qso_stacks_dark.npz v1_dark
log "v1_dark finished"
until [ $(grep -l '^done:' $S/qso_stack_desi?.log 2>/dev/null | wc -l) -ge 3 ] && ! pgrep -f '05_stack.py --survey desi' >/dev/null; do sleep 60; done
log "all DESI coadds done: final merge + EM-PCA v1"
python3 05_stack.py --merge && empca_all ~/data/sdssdesi_products/quasar/qso_stacks.npz v1
log "v1 finished"
