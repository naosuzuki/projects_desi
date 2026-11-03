#!/bin/zsh
# v3 run 5 (Nao 2026-11-02, 'apply the telluric correction the Stellar thread is finding and redo your stack and PCA'):
# the run-4 cells and bins, every member divided by the Stellar thread's telluric factor (programs_telluric/tellapply.py,
# model = the member's run-4 stack template at its amplitude), re-stack, luminosity plots, EM-PCA tagged v5.
#  1. run 4 must have finished (after_stack_v3.log) and tellapply.py must exist
#  2. run 4 kept aside as *_run4 (stacks, cells csv, template table, exdisk3 accumulators, logs); cells parquet unchanged
#  3. QSO_TEMPLATE=v3/qso_stacks_run4.npz QSO_TELLURIC=1 run_stack_v3.sh; QSO_EMTAG=v5 after_stack_v3.sh
# Log: ~/data/sdssdesi_status/qso_v3_telluric.log
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar
log() { echo "$(date '+%F %T') $*"; }
[ -f ../programs_telluric/tellapply.py ] || { log "no programs_telluric/tellapply.py yet"; exit 1; }
grep -q 'after_stack_v3 finished' $S/after_stack_v3.log 2>/dev/null || { log "run 4 not finished"; exit 1; }
if [ ! -f $P/v3/qso_stacks_run4.npz ]; then
  log "keeping run 4 aside as *_run4"
  mv $P/v3/qso_stacks.npz $P/v3/qso_stacks_run4.npz; mv $P/v3/qso_stacks_cells.csv $P/v3/qso_stacks_cells_run4.csv
  mv $P/v3/qso_norm_table_template.parquet $P/v3/qso_norm_table_template_run4.parquet
  mv /Volumes/exdisk3/sdssdesi_qso_v3_acc /Volumes/exdisk3/sdssdesi_qso_v3_acc_run4
  mv $S/qso_stack_v3.log $S/qso_stack_v3_run4.log; mv $S/after_stack_v3.log $S/after_stack_v3_run4.log
  for f in $S/qso_stack_v3_desi?.log $S/qso_stack_v3_sdss?.log; do mv $f ${f%.log}_run4.log; done
fi
log "stacks (telluric-corrected members, run-4 templates)"; QSO_TEMPLATE=$P/v3/qso_stacks_run4.npz QSO_TELLURIC=1 zsh run_stack_v3.sh >> $S/qso_stack_v3.log 2>&1
log "after_stack (v5)"; QSO_EMTAG=v5 zsh after_stack_v3.sh >> $S/after_stack_v3.log 2>&1
log "run_v3_telluric finished"
