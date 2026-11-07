#!/bin/zsh
# Stellar version 8 (hr_all.py): wait for the store, then lines, stack, pca and the figures. Log: ~/data/sdssdesi_status/hr_all_v8.log
cd "$(dirname "$0")"
LOG=~/data/sdssdesi_status/hr_all_v8.log
while pgrep -f 'hr_all.py store' >/dev/null; do sleep 60; done
HR_NPROC=3 python3 hr_all.py store >> $LOG 2>&1 || { echo "FAILED store" >> $LOG; exit 1; }      # picks up any coadd left over
for s in lines stack pca; do
  echo "=== $s $(date)" >> $LOG
  python3 hr_all.py $s >> $LOG 2>&1 || { echo "FAILED $s" >> $LOG; exit 1; }
done
echo "=== plots $(date)" >> $LOG
STAR_HR=all python3 hr_plots.py >> $LOG 2>&1 || { echo "FAILED plots" >> $LOG; exit 1; }
echo "=== ALL DONE $(date)" >> $LOG
