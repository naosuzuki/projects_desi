#!/bin/zsh
# Pre-shutdown: stop all pipeline processes cleanly and unmount the external disks, so the Mac can be powered off.
# Usage:  ./preshutdown.sh            run now
#         ./preshutdown.sh 04:15     wait until that clock time (today/tomorrow), then run (detached with nohup)
LOG=~/data/sdssdesi_status/shutdown.log
log() { echo "$(date '+%F %T') $*" | tee -a $LOG; }
if [ -n "$1" ]; then
  target=$(date -j -f '%H:%M' "$1" '+%s' 2>/dev/null); now=$(date +%s)
  [ "$target" -le "$now" ] && target=$((target + 86400))
  log "pre-shutdown armed for $(date -r $target '+%F %T') (in $(( (target-now)/60 )) min)"
  sleep $((target - now))
fi
log "=== pre-shutdown start ==="
log "state: passB outputs $(ls ~/data/sdssdesi_passB/*/ 2>/dev/null | grep -c parquet); passQ outputs $(ls ~/data/sdssdesi_passQ/*/ 2>/dev/null | grep -c parquet); expfibermap $(ls ~/data/sdssdesi_expfibermap/*/ 2>/dev/null | grep -c parquet)"
# 1. stop pipeline processes (SIGTERM; outputs are per-file atomic, in-flight files are redone on resume)
for pat in 'python3 02_passQ_desi.py --all' 'python3 02_passB_measure.py --all' 'python3 08_passC_stack.py' 'python3 11_passG_sdss.py' 'autoscale_passB.sh' 'after_passB.sh' 'monitor.py --port' '04_monitor.py'; do
  n=$(pgrep -f "$pat" | wc -l | tr -d ' ')
  if [ "$n" -gt 0 ]; then pkill -TERM -f "$pat"; log "SIGTERM -> $n process(es) matching '$pat'"; fi
done
# graceful stop: the walk (with the galaxy hook) and pass C finish the current coadd and checkpoint -> wait up to 3 min
for i in $(seq 1 36); do pgrep -f 'python3 (02_passQ_desi|02_passB_measure|08_passC_stack|11_passG_sdss)' >/dev/null || break; sleep 5; done
for pat in 'python3 02_passQ_desi.py --all' 'python3 02_passB_measure.py --all' 'python3 08_passC_stack.py' 'python3 11_passG_sdss.py'; do
  n=$(pgrep -f "$pat" | wc -l | tr -d ' '); [ "$n" -gt 0 ] && { pkill -KILL -f "$pat"; log "SIGKILL -> $n stubborn process(es) '$pat'"; }
done
pkill -x caffeinate && log "caffeinate holds released"
# 2. flush and unmount the external disks (they are read-only for us, but a clean unmount avoids a journal replay)
sync; sleep 2
for v in exdisk2 exdisk1 exdisk3; do
  if mount | grep -q " /Volumes/$v "; then
    if diskutil unmount "/Volumes/$v" >> $LOG 2>&1; then log "unmounted /Volumes/$v"; else sleep 15; diskutil unmount force "/Volumes/$v" >> $LOG 2>&1 && log "force-unmounted /Volumes/$v" || log "FAILED to unmount /Volumes/$v - check lsof"; fi
  else log "/Volumes/$v was not mounted"; fi
done
log "remaining pipeline processes: $(pgrep -fl 'passQ|passB|monitor.py|autoscale' | wc -l | tr -d ' ')"
log "=== ready for shutdown. After reboot: see docs/RESTART_after_shutdown.md ==="
