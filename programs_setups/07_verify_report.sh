#!/bin/zsh
# 07_verify_report.sh <task_id> <disk> — wait for one verify task to finish, then log its counts and save the list of
# files it re-downloaded (corrupted copies replaced). Used when 06_verify_redownload.sh's waiter was killed
# (2026-10-17 power outage). Survives nothing itself: rerun after a reboot.
HERE=${0:A:h}; G=~/.local/bin/globus; T=$1; DISK=${2:-exdisk2}; LOG=$HERE/logs/verify.log
say() { print -r -- "$(date '+%F %T') $*" >> $LOG }
until st=$($G task show $T -F json 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["status"])' 2>/dev/null); [ -n "$st" ] && [ "$st" != ACTIVE ] && [ "$st" != INACTIVE ]; do sleep 900; done
info=$($G task show $T -F json | python3 -c 'import json,sys;d=json.load(sys.stdin);print(d["files"],d["files_transferred"],d["files_skipped"],d["faults"])')
$G api transfer GET /task/$T/successful_transfers -Q limit=1000 2>/dev/null | python3 -c '
import json,sys
print("\n".join(x["destination_path"] for x in json.load(sys.stdin).get("DATA",[])))' > $HERE/logs/verify_redownloaded_$DISK.txt
say "verify $DISK task $T: $st  (files, re-downloaded, skipped-identical, faults) = $info"
say "VERIFY DONE"
