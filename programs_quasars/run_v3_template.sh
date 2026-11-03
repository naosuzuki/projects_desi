#!/bin/zsh
# v3 re-stack with the template-amplitude normalization (Nao 2026-11-01 'Yes, template norm'): each member normalized by the
# least-squares amplitude of its stack's mean spectrum (the window-normalized v3 stacks, kept as *_winnorm) over all its
# pixels; then after_stack_v3.sh (plots with the sparse edges dropped, paper figures, pivoted EM-PCA with --min-total).
# Log: ~/data/sdssdesi_status/qso_v3_template.log. Resumable.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar
log() { echo "$(date '+%F %T') $*"; }
if [ ! -f $P/v3/qso_stacks_winnorm.npz ]; then
  log "keeping the window-normalized v3 aside as *_winnorm"
  mv $P/v3/qso_stacks.npz $P/v3/qso_stacks_winnorm.npz; mv $P/v3/qso_stacks_cells.csv $P/v3/qso_stacks_cells_winnorm.csv
  mv $P/qso_empca_v3.fits $P/qso_empca_v3_winnorm.fits; mv $P/qso_empca_v3.png $P/qso_empca_v3_winnorm.png
  for k in simple ew slope; do mv $P/checks/qso_lum_stacks_v3_$k.png $P/checks/qso_lum_stacks_v3_${k}_winnorm.png; mv $P/checks/qso_lum_stacks_v3_$k.csv $P/checks/qso_lum_stacks_v3_${k}_winnorm.csv; done
  mv /Volumes/exdisk3/sdssdesi_qso_v3_acc /Volumes/exdisk3/sdssdesi_qso_v3_acc_winnorm
  mv $S/qso_stack_v3.log $S/qso_stack_v3_run3_winnorm.log; mv $S/after_stack_v3.log $S/after_stack_v3_run3_winnorm.log
  for f in $S/qso_stack_v3_desi?.log $S/qso_stack_v3_sdss?.log; do mv $f ${f%.log}_run3_winnorm.log; done
fi
log "stacks (template normalization)"; QSO_TEMPLATE=$P/v3/qso_stacks_winnorm.npz zsh run_stack_v3.sh >> $S/qso_stack_v3.log 2>&1
log "after_stack"; zsh after_stack_v3.sh >> $S/after_stack_v3.log 2>&1
log "run_v3_template finished"
