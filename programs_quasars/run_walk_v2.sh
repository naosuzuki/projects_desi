#!/bin/zsh
# v2 shared DESI walk (Nao, 2026-10-29): ONE read of every coadd for galaxies (galaxy_hook v2: exact sub-pixel rest frame
# -> /Volumes/exdisk3/sdssdesi_restframe_v2, ~/data/sdssdesi_passC_v2) and quasars + stars (native-pixel store
# -> /Volumes/exdisk3/sdssdesi_native/desi). Pass-Q outputs exist, so no quasar re-measurement. Resumable.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status
for i in 0 1 2 3; do
  if pgrep -f "02_passQ_desi.py --all --with-native.*--shard $i/4" >/dev/null; then echo "shard $i running"; continue; fi
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --with-native --with-galaxies --gal-store /Volumes/exdisk3/sdssdesi_restframe_v2 --shard $i/4 >> $S/walk_v2_$i.log 2>&1 &
  nohup caffeinate -i -w $! >/dev/null 2>&1 &
  sleep 3
done
