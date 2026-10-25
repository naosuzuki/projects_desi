#!/bin/zsh
# Resume pass B with N workers (default 4), keep-awake holds, and the monitor if it is not running.
# Usage: ./resume_passB.sh [N]     (requires /Volumes/exdisk2 mounted; existing outputs are skipped)
set -e
N=${1:-4}
cd ~/github/projects_desi/programs_templates
if ! mount | grep -q ' /Volumes/exdisk2 '; then echo "exdisk2 is not mounted - aborting"; exit 1; fi
if ! ls /Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main >/dev/null 2>&1; then echo "exdisk2 mounted but coadd tree unreadable - aborting"; exit 1; fi
if pgrep -f '02_passB_measure.py --all' >/dev/null; then echo "workers already running:"; pgrep -fl '02_passB_measure.py --all'; exit 1; fi
mkdir -p ~/data/sdssdesi_status
echo "$(date '+%F %T') resume with $N workers; done so far: $(ls ~/data/sdssdesi_passB/*/ 2>/dev/null | grep -c parquet)" | tee -a ~/data/sdssdesi_status/resume.log
for ((i=0; i<N; i++)); do
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passB_measure.py --all --shard $i/$N >> ~/data/sdssdesi_status/passB_w${i}of${N}.log 2>&1 &
  sleep 3
done
sleep 2
for p in $(pgrep -f '02_passB_measure.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
if ! pgrep -f 04_monitor.py >/dev/null; then nohup python3 04_monitor.py --port 8765 >> ~/data/sdssdesi_status/monitor.log 2>&1 & fi
sleep 1
echo "workers: $(pgrep -f '02_passB_measure.py --all' | wc -l | tr -d ' ')  caffeinate: $(pgrep -x caffeinate | wc -l | tr -d ' ')  monitor: http://127.0.0.1:8765/monitor.html"
