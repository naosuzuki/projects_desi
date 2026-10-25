#!/bin/zsh
# Unattended autoscaler for pass B. Every 2 min decide how many workers should run from
#   - the exdisk2 USB link speed (USB 2 = 480 Mb/s: the disk is the bottleneck, 2 workers max;
#     USB 3 = 5 Gb/s: 4 workers, or 3 while the Prospector job is still using the performance cores)
#   - whether the Prospector job (projects_hsc 21_prosp_fit.py) is still running
# and restart the workers through resume_passB.sh when the target changes. Exits when pass B is complete.
# Log: ~/data/sdssdesi_status/autoscale.log.   Start: nohup ./autoscale_passB.sh > /dev/null 2>&1 &
cd ~/github/projects_desi/programs_templates
LOG=~/data/sdssdesi_status/autoscale.log
TOTAL=34508
log() { echo "$(date '+%F %T') $*" | tee -a $LOG; }
current_n() { pgrep -f '02_passB_measure.py --all' | wc -l | tr -d ' '; }
link() { system_profiler SPUSBDataType 2>/dev/null | awk 'BEGIN{RS="\n\n"} /BSD Name: disk5/ {print}' | grep -m1 Speed | sed 's/.*Up to //'; }
log "autoscaler started (current workers: $(current_n))"
last_target=0
while true; do
  done_n=$(ls ~/data/sdssdesi_passB/*/ 2>/dev/null | grep -c parquet)
  if [ "$done_n" -ge "$TOTAL" ]; then log "pass B complete ($done_n files) - autoscaler exiting"; exit 0; fi
  if ! mount | grep -q ' /Volumes/exdisk2 '; then log "exdisk2 not mounted - waiting"; sleep 120; continue; fi
  sp=$(link); prosp=$(pgrep -f 21_prosp_fit | wc -l | tr -d ' ')
  case "$sp" in *Gb*) usb3=1;; *) usb3=0;; esac
  # USB 3: 4 workers (3 while Prospector holds the performance cores). USB 2: the link is the limit, but a
  # 10-min test gave 18/min with 4 workers vs 13-16/min with 2, so use 4 once the CPU is free, else 2.
  if [ $usb3 -eq 1 ]; then if [ "$prosp" -gt 0 ]; then target=3; else target=4; fi
  else if [ "$prosp" -gt 0 ]; then target=2; else target=4; fi; fi
  n=$(current_n)
  if [ "$target" -ne "$last_target" ] || [ "$n" -eq 0 ]; then
    log "link=${sp:-?} prospector=$prosp workers=$n -> target=$target (done $done_n/$TOTAL)"
    if [ "$n" -ne "$target" ]; then
      for p in $(pgrep -f '02_passB_measure.py --all'); do kill -TERM $p; done; sleep 5; pkill -x caffeinate
      ./resume_passB.sh $target >> $LOG 2>&1 || log "resume failed"
    fi
    last_target=$target
  fi
  sleep 120
done
