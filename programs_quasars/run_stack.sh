#!/bin/zsh
# Quasar stacks (05_stack.py): 3 SDSS shards (internal disk) + 3 DESI shards (exdisk2, read-only). Resumable: rerun
# after a stop/reboot (checkpoints every 100 files; SIGTERM finishes the current file and saves). Then: 05_stack.py --merge
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status
for sv in sdss desi; do for i in 0 1 2; do
  if pgrep -f "05_stack.py --survey $sv --shard $i/3" >/dev/null; then echo "$sv $i/3 already running"; continue; fi
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 05_stack.py --survey $sv --shard $i/3 >> $S/qso_stack_${sv}$i.log 2>&1 &
  sleep 2
done; done
for p in $(pgrep -f '05_stack.py --survey'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
echo "stack workers: $(pgrep -f '05_stack.py --survey' | wc -l | tr -d ' ')"
