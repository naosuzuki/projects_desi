#!/bin/zsh
# Launch (or resume) pass C with N workers, after pass Q DESI has finished with the coadd disk.
# Usage: ./run_passC.sh [N] [STORE_DIR]       e.g.  ./run_passC.sh 4 /Volumes/exdisk3/sdssdesi_restframe
# Resumable: each shard skips the files already in its accumulator checkpoint (acc/acc_shard<i>of<N>.npz);
# keep N the same across restarts (the shard split depends on N).
set -e
N=${1:-4}; STORE=${2:-}
cd ~/github/projects_desi/programs_templates
mount | grep -q ' /Volumes/exdisk2 ' || { echo "exdisk2 not mounted"; exit 1; }
sp=$(system_profiler SPUSBDataType 2>/dev/null | awk 'BEGIN{RS="\n\n"} /25EE|25ee/ {print}' | grep -m1 Speed | sed 's/.*Up to //')
case "$sp" in *Gb*) ;; *) [ -z "$FORCE" ] && { echo "exdisk2 link is '$sp' (not USB 3) - aborting; FORCE=1 to override"; exit 1; };; esac
if pgrep -f 'python3 02_passQ_desi.py --all' >/dev/null; then echo "pass Q still running - pass C waits (one walk at a time)"; exit 1; fi
if pgrep -f 'python3 08_passC_stack.py' >/dev/null; then echo "pass C already running"; exit 1; fi
mkdir -p ~/data/sdssdesi_status
STOREARG=(); [ -n "$STORE" ] && { mkdir -p "$STORE"; STOREARG=(--store "$STORE"); }
echo "$(date '+%F %T') pass C start: $N workers, store=${STORE:-none}" | tee -a ~/data/sdssdesi_status/passC.log
for ((i=0; i<N; i++)); do
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 08_passC_stack.py --all --shard $i/$N $STOREARG >> ~/data/sdssdesi_status/passC_w${i}of${N}.log 2>&1 &
  sleep 3
done
sleep 2
for p in $(pgrep -f 'python3 08_passC_stack.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
echo "workers: $(pgrep -f 'python3 08_passC_stack.py --all' | wc -l | tr -d ' ')  caffeinate: $(pgrep -x caffeinate | wc -l | tr -d ' ')"
