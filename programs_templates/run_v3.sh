#!/bin/zsh
# Galaxy release v3 re-read (Nao 2026-10-31): DESI coadds (exdisk2, read-only) -> both grids (MABS_G, M_AB(5450)), 4 shards.
# Usage: ./run_v3.sh 0 1      start the given shards (resumable; keep N=4)
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; mkdir -p ~/data/sdssdesi_passC_v3
for i in "$@"; do
  if pgrep -f "08_passC_stack.py --all --v3 --shard $i/4" >/dev/null; then echo "shard $i running"; continue; fi
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 08_passC_stack.py --all --v3 --shard $i/4 >> $S/gal_v3_desi_$i.log 2>&1 &
  nohup caffeinate -i -w $! >/dev/null 2>&1 &
  echo "$(date '+%F %T') v3 DESI shard $i/4 started (pid $!)" | tee -a $S/gal_v3.log
  sleep 2
done
