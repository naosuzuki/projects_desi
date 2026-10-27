#!/bin/zsh
# Resume the SDSS-DESI pipeline after the planned reboot (2026-10-26). Run once, from any terminal:
#   ~/data/sdssdesi_status/resume_after_reboot.sh
# 1. checks the external disks (exdisk2 = coadds, exdisk1 = overflow), 2. starts the monitor + guard on 8765,
# 3. relaunches the combined DESI walk (quasars + stars + exposure map; 4 shards; existing outputs skipped),
# 4. keep-awake holds. Galaxy pass C and anything else are started by their own threads afterwards.
LOG=~/data/sdssdesi_status/resume.log
if ! mount | grep -q ' /Volumes/exdisk2 '; then echo "exdisk2 not mounted - plug it in (USB 3 port!) and rerun"; exit 1; fi
if ! mount | grep -q ' /Volumes/exdisk1 ' && ! mount | grep -q ' /Volumes/My Book '; then echo "WARNING: overflow disk (exdisk1 / My Book) not mounted - 969 coadds will be reported MISSING"; fi
sp=$(system_profiler SPUSBDataType 2>/dev/null | awk 'BEGIN{RS="\n\n"} /My Book 25EE/ {print}' | grep -m1 Speed | sed 's/.*Up to //'); echo "exdisk2 link: ${sp:-?}"
case "$sp" in *Gb*) ;; *) echo "ABORT: exdisk2 is not on a USB 3 link (${sp:-unknown}); re-plug it on a USB 3 port and rerun (FORCE=1 to override; a USB 2 relaunch costs 3x on the remaining walk)"; [ "${FORCE:-0}" = "1" ] || exit 1;; esac
if ! ls /Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main >/dev/null 2>&1; then echo "ABORT: coadd tree unreadable on exdisk2"; exit 1; fi
if ! pgrep -f 'monitor.py --port 8765' >/dev/null; then
  (cd ~/github/projects_desi/programs_monitor && nohup python3 monitor.py --port 8765 >> ~/data/sdssdesi_status/monitor.log 2>&1 &)
  echo "monitor started: http://127.0.0.1:8765/monitor.html"
fi
if pgrep -f 'python3 02_passQ_desi.py --all' >/dev/null; then echo "walk already running"; exit 0; fi
cd ~/github/projects_desi/programs_quasars
echo "$(date '+%F %T') resume walk with 4 shards; done so far Q=$(ls ~/data/sdssdesi_passQ/*/ 2>/dev/null | grep -c parquet)" | tee -a $LOG
for i in 0 1 2 3; do
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --with-stars --expfibermap --with-galaxies --shard $i/4 >> ~/data/sdssdesi_status/passQ_desi$i.log 2>&1 &
  sleep 3
done
for p in $(pgrep -f 'python3 02_passQ_desi.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
echo "workers: $(pgrep -f 'python3 02_passQ_desi.py --all' | wc -l | tr -d ' ')"
