#!/bin/zsh
# v3 re-stack with the measured 3000-3050 A normalization (Nao 2026-10-31 'Yes, re-stack'):
#  1. measurement pass (05_stack.py --v3 --measure, 4 DESI + 2 SDSS shards, exdisk3) -> v3/qso_norm_table.parquet
#  2. the fit-chained v3 products are kept aside as *_fitchain (stacks, cells, EM-PCA, plots, accumulators)
#  3. 04_cells.py --survey all --chain (measured ratio wherever covered, fitted R(M) otherwise)
#  4. run_stack_v3.sh + after_stack_v3.sh (stacks, luminosity plots, paper figures, EM-PCA)
# Log: ~/data/sdssdesi_status/qso_v3_measured.log. Resumable: rerun; done shards / steps are skipped.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar; M=/Volumes/exdisk3/sdssdesi_qso_v3_meas
log() { echo "$(date '+%F %T') $*"; }
mkdir -p $M
if [ ! -f $P/v3/qso_norm_table.parquet ]; then
  log "measurement pass start"
  for i in 0 1 2 3; do
    OPENBLAS_NUM_THREADS=1 nice -n 10 python3 05_stack.py --v3 --measure --survey desi --shard $i/4 --accdir $M >> $S/qso_meas_v3_desi$i.log 2>&1 &
  done
  for i in 0 1; do
    OPENBLAS_NUM_THREADS=1 nice -n 10 python3 05_stack.py --v3 --measure --survey sdss --shard $i/2 --accdir $M >> $S/qso_meas_v3_sdss$i.log 2>&1 &
  done
  wait
  n=$(grep -l '^done:' $S/qso_meas_v3_desi?.log $S/qso_meas_v3_sdss?.log 2>/dev/null | wc -l)
  if [ $n -lt 6 ]; then log "only $n of 6 measurement shards done; stopping"; exit 1; fi
  log "measurement merge"; python3 05_stack.py --v3 --measure --merge --accdir $M || exit 1
fi
if [ ! -f $P/v3/qso_stacks_fitchain.npz ]; then
  log "keeping the fit-chained v3 aside as *_fitchain"
  mv $P/v3/qso_stacks.npz $P/v3/qso_stacks_fitchain.npz; mv $P/v3/qso_stacks_cells.csv $P/v3/qso_stacks_cells_fitchain.csv
  mv $P/qso_cells_v3_all.parquet $P/qso_cells_v3_all_fitchain.parquet; mv $P/qso_cells_v3_all.txt $P/qso_cells_v3_all_fitchain.txt 2>/dev/null
  mv $P/qso_empca_v3.fits $P/qso_empca_v3_fitchain.fits; mv $P/qso_empca_v3.png $P/qso_empca_v3_fitchain.png
  for k in simple ew slope; do mv $P/checks/qso_lum_stacks_v3_$k.png $P/checks/qso_lum_stacks_v3_${k}_fitchain.png; mv $P/checks/qso_lum_stacks_v3_$k.csv $P/checks/qso_lum_stacks_v3_${k}_fitchain.csv; done
  mv /Volumes/exdisk3/sdssdesi_qso_v3_acc /Volumes/exdisk3/sdssdesi_qso_v3_acc_fitchain
  mv $S/qso_stack_v3.log $S/qso_stack_v3_run2_fitchain.log; mv $S/after_stack_v3.log $S/after_stack_v3_run2_fitchain.log
  for f in $S/qso_stack_v3_desi?.log $S/qso_stack_v3_sdss?.log; do mv $f ${f%.log}_run2_fitchain.log; done
fi
if [ ! -f $P/qso_cells_v3_all.parquet ]; then
  log "cells (measured chaining)"; python3 04_cells.py --survey all --chain > $S/qso_cells_v3_measured.log 2>&1 || exit 1
  tail -4 $S/qso_cells_v3_measured.log
fi
log "stacks"; zsh run_stack_v3.sh >> $S/qso_stack_v3.log 2>&1
log "after_stack"; zsh after_stack_v3.sh >> $S/after_stack_v3.log 2>&1
log "run_v3_measured finished"
