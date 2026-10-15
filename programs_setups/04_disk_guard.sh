#!/bin/zsh
# 04_disk_guard.sh — stop the DESI DR2 Globus transfer before exdisk2 is full (Nao, 2026-10-14: "set the stop at 50 GB").
#
#   zsh 04_disk_guard.sh [min_free_GB=50] [task_id]
#
# Every 60 s: if free space on /Volumes/exdisk2 < min_free_GB, cancel the Globus task (logged, then exit).
# Exits by itself when the task is no longer ACTIVE (SUCCEEDED / FAILED / cancelled elsewhere).
# A cancelled task is resumed later by re-running 02_transfer_healpix.sh (checksum sync skips files on disk).
# Log: logs/disk_guard.log
MIN_GB=${1:-50}
TASK=${2:-07cf18f4-c48e-11f1-b8f0-0affd5e180af}
G=~/.local/bin/globus
LOG=${0:A:h}/logs/disk_guard.log
mkdir -p ${LOG:h}
say() { print -r -- "$(date '+%F %T') $*" >> $LOG }
tmo() { local s=$1; shift; "$@" & local p=$!; ( sleep $s; kill $p ) >/dev/null 2>&1 & local w=$!
  wait $p; local rc=$?; kill $w 2>/dev/null; return $rc }

say "start: cancel task $TASK when exdisk2 free < $MIN_GB GB"
n=0
while true; do
  free_gb=$(df -k /Volumes/exdisk2 2>/dev/null | tail -1 | awk '{printf "%d", $4*1024/1e9}')
  if [ -z "$free_gb" ]; then say "exdisk2 not mounted?"; sleep 60; continue; fi
  if [ $free_gb -lt $MIN_GB ]; then
    say "free $free_gb GB < $MIN_GB GB: cancelling task"
    if tmo 120 $G task cancel $TASK >> $LOG 2>&1; then
      say "cancelled. files on disk: $(find /Volumes/exdisk2/desiredux -name '*.fits' | wc -l | tr -d ' ')"
      exit 0
    fi
    say "cancel failed (API?); retrying in 60 s"
  fi
  if (( n % 10 == 0 )); then                     # task status every 10 min
    st=$(tmo 60 $G task show $TASK -F json 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["status"])' 2>/dev/null)
    say "free $free_gb GB, task ${st:-API?}"
    if [ -n "$st" ] && [ "$st" != ACTIVE ]; then say "task $st: guard exits"; exit 0; fi
  fi
  n=$((n + 1))
  sleep 60
done
