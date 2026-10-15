#!/bin/zsh
# 03_monitor_transfer.sh — watch the DESI DR2 healpix Globus transfer to exdisk2.
#
#   zsh 03_monitor_transfer.sh [minutes=60] [task_id]
#
# Every 10 min appends one line to logs/transfer_rate.log:
#   time, files on disk, GB used / free on exdisk2, Globus bytes, MB/s (disk growth, last 10 min), task status
#   (API = Globus API unreachable), faults, hours until only 300 GB are left
# Exits after <minutes>, or early (with PROBLEM:) if the task stops being ACTIVE,
# new faults appear, or nothing arrives for 30 min.
setopt NULL_GLOB
mins=${1:-60}
TASK=${2:-07cf18f4-c48e-11f1-b8f0-0affd5e180af}
G=~/.local/bin/globus
DEST=/Volumes/exdisk2/desiredux
LOG=${0:A:h}/logs/transfer_rate.log
mkdir -p ${LOG:h}

tmo() { local s=$1; shift; "$@" & local p=$!; ( sleep $s; kill $p ) >/dev/null 2>&1 & local w=$!
  wait $p; local rc=$?; kill $w 2>/dev/null; return $rc }      # macOS has no timeout(1)
# 2026-10-13: transfer.api.globus.org stopped answering while the data kept flowing -> on an API
# failure fall back to disk growth (status API, bytes -1) instead of exiting with garbage.
read_task() { local j; j=$(tmo 60 $G task show $TASK -F json 2>/dev/null)
  st=$(print -r -- $j | python3 -c 'import json,sys;d=json.load(sys.stdin);print(d["status"],d["bytes_transferred"],d["faults"],d["request_time"])' 2>/dev/null)
  [ -n "$st" ] || st="API -1 -1 -" }
dused() { df -k /Volumes/exdisk2 | tail -1 | awk '{print $3*1024}' }

read_task; set -- ${=st}; f0=$3; t0=$(date +%s); d0=$(dused)
[ $f0 -lt 0 ] && f0=0
idle=0
for i in $(seq 1 $((mins / 10))); do
  sleep 600
  read_task; set -- ${=st}; tstat=$1; bytes=$2; faults=$3
  now=$(date +%s); d1=$(dused)
  r10=$(( (d1 - d0) / (now - t0) / 1000000 ))            # rate from disk growth (works without the API)
  nf=$(find $DEST -name '*.fits' | wc -l | tr -d ' ')
  used=$(df -g /Volumes/exdisk2 | tail -1 | awk '{print $3}')
  avail=$(df -k /Volumes/exdisk2 | tail -1 | awk '{print $4*1024}')
  full=$(python3 -c "r=max($r10,1)*1e6;h=($avail-3e11)/r/3600;print(f'{h:.0f} h')")
  gb=$([ $bytes -ge 0 ] && echo $((bytes/1000000000)) || echo n/a)
  line="$(date '+%F %T')  files $nf  disk ${used} GB  free $((avail/1000000000)) GB  globus ${gb} GB  ${r10} MB/s (10 min)  task $tstat  faults $faults  300GB-left in $full"
  print -r -- "$line" >> $LOG
  [ "$tstat" = ACTIVE ] || [ "$tstat" = API ] || { echo "PROBLEM: task status $tstat"; echo "$line"; exit 0; }
  [ $faults -gt $f0 ] && { echo "PROBLEM: faults $f0 -> $faults"; echo "$line"; tmo 60 $G task event-list $TASK | head -5; exit 0; }
  if [ $d1 -le $d0 ]; then idle=$((idle + 1)); else idle=0; fi
  [ $idle -ge 3 ] && { echo "PROBLEM: no bytes for 30 min"; echo "$line"; exit 0; }
  d0=$d1; t0=$now
done
tail -n $((mins / 10)) $LOG
