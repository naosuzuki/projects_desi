#!/bin/zsh
# v3 run 4 (2026-11-01, after the template-normalized run 3): luminosity bins from the corrected template amplitudes
# (14_template_rebin.py), members normalized by the amplitude of the mixed templates (fit-chained where they exist).
#  1. keep run 3 aside as *_tplwin (stacks, cells, EM-PCA, plots, exdisk3 accumulators, logs)
#  2. 04_cells.py --survey all --chain --template-bins --template-table v3/qso_norm_table_template_corr.parquet
#  3. run_stack_v3.sh with QSO_TEMPLATE=v3/qso_stacks_tplmix.npz, then after_stack_v3.sh
# Log: ~/data/sdssdesi_status/qso_v3_rebin.log. 14_template_rebin.py must have been run first (it reads the run-3 products).
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar
log() { echo "$(date '+%F %T') $*"; }
[ -f $P/v3/qso_norm_table_template_corr.parquet ] && [ -f $P/v3/qso_stacks_tplmix.npz ] || { log "run 14_template_rebin.py first"; exit 1; }
if [ ! -f $P/v3/qso_stacks_tplwin.npz ]; then
  log "keeping run 3 aside as *_tplwin"
  mv $P/v3/qso_stacks.npz $P/v3/qso_stacks_tplwin.npz; mv $P/v3/qso_stacks_cells.csv $P/v3/qso_stacks_cells_tplwin.csv
  mv $P/v3/qso_norm_table_template.parquet $P/v3/qso_norm_table_template_tplwin.parquet
  mv $P/qso_cells_v3_all.parquet $P/qso_cells_v3_all_tplwin.parquet; mv $P/qso_cells_v3_all.txt $P/qso_cells_v3_all_tplwin.txt 2>/dev/null
  mv $P/qso_empca_v3.fits $P/qso_empca_v3_tplwin.fits; mv $P/qso_empca_v3.png $P/qso_empca_v3_tplwin.png
  for k in simple ew slope; do mv $P/checks/qso_lum_stacks_v3_$k.png $P/checks/qso_lum_stacks_v3_${k}_tplwin.png; mv $P/checks/qso_lum_stacks_v3_$k.csv $P/checks/qso_lum_stacks_v3_${k}_tplwin.csv; done
  mv /Volumes/exdisk3/sdssdesi_qso_v3_acc /Volumes/exdisk3/sdssdesi_qso_v3_acc_tplwin
  mv $S/qso_stack_v3.log $S/qso_stack_v3_run3_tplwin.log; mv $S/after_stack_v3.log $S/after_stack_v3_run3_tplwin.log
  for f in $S/qso_stack_v3_desi?.log $S/qso_stack_v3_sdss?.log; do mv $f ${f%.log}_run3_tplwin.log; done
fi
if [ ! -f $P/qso_cells_v3_all.parquet ]; then
  log "cells (template-amplitude bins)"; python3 04_cells.py --survey all --chain --template-bins --template-table $P/v3/qso_norm_table_template_corr.parquet > $S/qso_cells_v3_rebin.log 2>&1 || exit 1
  tail -4 $S/qso_cells_v3_rebin.log
fi
log "stacks (mixed templates)"; QSO_TEMPLATE=$P/v3/qso_stacks_tplmix.npz zsh run_stack_v3.sh >> $S/qso_stack_v3.log 2>&1
log "after_stack"; zsh after_stack_v3.sh >> $S/after_stack_v3.log 2>&1
log "run_v3_rebin finished"
