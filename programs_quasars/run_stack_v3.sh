#!/bin/zsh
# v3 quasar stacks (Nao 2026-10-31): no redshift axis, every member chained to 3000-3050 A with the DESI-fitted ratios
# (04_cells.py --chain -> qso_cells_v3_all; 05_stack.py --v3). DESI from the native store on exdisk3 (4 shards), SDSS from
# the spPlates (2 shards), accumulators on exdisk3, then merge -> quasar/v3/qso_stacks.npz. Resumable. Log: qso_stack_v3.log.
# QSO_TEMPLATE=<stacks npz> (2026-11-01): template-amplitude normalization (05_stack.py --template) instead of the window mean.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; ACC=/Volumes/exdisk3/sdssdesi_qso_v3_acc
TPL=${QSO_TEMPLATE:+--template $QSO_TEMPLATE}
# QSO_TELLURIC=1 (run 5, 2026-11-02): the Stellar thread's telluric factor divided out of every member (05_stack.py --telluric, needs the template)
TPL="$TPL ${QSO_TELLURIC:+--telluric}"
log() { echo "$(date '+%F %T') $*"; }
mkdir -p $ACC
log "v3 stacks start"
for i in 0 1 2 3; do
  OPENBLAS_NUM_THREADS=1 nice -n 10 python3 05_stack.py --v3 --survey desi --shard $i/4 --accdir $ACC ${=TPL} >> $S/qso_stack_v3_desi$i.log 2>&1 &
done
for i in 0 1; do
  OPENBLAS_NUM_THREADS=1 nice -n 10 python3 05_stack.py --v3 --survey sdss --shard $i/2 --accdir $ACC ${=TPL} >> $S/qso_stack_v3_sdss$i.log 2>&1 &
done
wait
n=$(grep -l '^done:' $S/qso_stack_v3_desi?.log $S/qso_stack_v3_sdss?.log 2>/dev/null | wc -l)
if [ $n -lt 6 ]; then log "only $n of 6 shards done; not merging"; exit 1; fi
log "merge"
python3 05_stack.py --v3 --merge --accdir $ACC
log "v3 stacks finished"
