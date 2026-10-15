#!/bin/zsh
# 06_verify_redownload.sh — after the download (exdisk2 + overflow on exdisk1), check every DESI DR2 coadd against
# NERSC and re-fetch the corrupted ones (Nao, 2026-10-15: "check if every file is intact ... If they are corrupted,
# erase them and redownload them").
#
#   zsh 06_verify_redownload.sh        # waits for the overflow task (logs/overflow_exdisk1.log) to finish
#
# Method: one Globus task per disk with --sync-level checksum over the files that disk holds. Globus computes the
# checksum at NERSC and locally; a file whose checksum differs (corrupted) is transferred again over the bad copy,
# i.e. erased and re-downloaded; files that match are left alone. --verify-checksum checks each re-download.
# Batches only list files present on that disk, so nothing new lands on the full exdisk2.
# Result: logs/verify.log (+ logs/verify_redownloaded_<disk>.txt = the files that were re-fetched).
HERE=${0:A:h}
G=~/.local/bin/globus
SRC=9d6d994a-6d04-11e5-ba46-22000b92c6ec
DST=$(cat ~/.globusonline/lta/client-id.txt)
LOG=$HERE/logs/verify.log
say() { print -r -- "$(date '+%F %T') $*" >> $LOG }
status() { $G task show $1 -F json 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["status"])' 2>/dev/null }

# 2026-10-15: rerunnable after a reboot. Older verify tasks are cancelled first (Globus' default 24 h deadline is too
# short for ~24 TB of checksums), and the new ones get a 10-day deadline.
for st in ACTIVE INACTIVE; do
  $G task list --filter-status $st -F json 2>/dev/null | python3 -c '
import json,sys
for t in json.load(sys.stdin)["DATA"]:
    if (t.get("label") or "").startswith("DESI DR2 verify"): print(t["task_id"])' | while read t; do
    $G task cancel $t >/dev/null 2>&1 && say "cancelled old verify task $t"
  done
done
DEADLINE=$(date -u -v+10d '+%Y-%m-%dT%H:%M:%S')
say "waiting for the overflow transfer"
until grep -qE "submitted overflow task|nothing to send" $HERE/logs/overflow_exdisk1.log 2>/dev/null; do sleep 300; done
T2=$(grep -o "submitted overflow task [0-9a-f-]*" $HERE/logs/overflow_exdisk1.log | awk '{print $4}')
if [ -n "$T2" ]; then
  until st=$(status $T2); [ -n "$st" ] && [ "$st" != ACTIVE ] && [ "$st" != INACTIVE ]; do sleep 300; done
  say "overflow task $T2: $st"
fi

python3 - "$HERE" <<'EOF' >> $LOG
import os, sys
here = sys.argv[1]
rows = [l.split() for l in open(f'{here}/manifests/healpix_main_darkbright.batch')]
ov = f'{here}/manifests/overflow_exdisk1.batch'
import shlex
over = [tuple(shlex.split(l)) for l in open(ov) if l.strip()] if os.path.exists(ov) else []
oset = {s for s, d in over}
a = [(s, d) for s, d in rows if s not in oset and os.path.exists(d)]
b = [(s, d) for s, d in over if os.path.exists(d)]
open(f'{here}/manifests/verify_exdisk2.batch', 'w').write(''.join(f'{s} {d}\n' for s, d in a))
open(f'{here}/manifests/verify_exdisk1.batch', 'w').write(''.join(f'"{s}" "{d}"\n' for s, d in b))
miss = [s for s, d in rows if s not in oset and not os.path.exists(d)] + [s for s, d in over if not os.path.exists(d)]
print(f'verify: exdisk2 {len(a)} files, exdisk1 {len(b)} files; missing everywhere {len(miss)}')
open(f'{here}/logs/verify_missing_everywhere.txt', 'w').write('\n'.join(miss) + '\n')
EOF

typeset -A TASKS
for disk in exdisk2 exdisk1; do
  B=$HERE/manifests/verify_$disk.batch
  [ -s $B ] || continue
  t=$($G transfer --batch $B --label "DESI DR2 verify $disk (checksum sync)" --sync-level checksum --deadline $DEADLINE \
       --preserve-timestamp --verify-checksum --jmespath 'task_id' --format unix $SRC $DST)
  TASKS[$disk]=$t
  say "submitted verify task $t for $disk ($(wc -l < $B | tr -d ' ') files)"
done

for disk t in ${(kv)TASKS}; do
  until st=$(status $t); [ -n "$st" ] && [ "$st" != ACTIVE ] && [ "$st" != INACTIVE ]; do sleep 600; done
  info=$($G task show $t -F json | python3 -c 'import json,sys;d=json.load(sys.stdin);print(d["files"],d["files_transferred"],d["files_skipped"],d["faults"])')
  # files actually (re)transferred = checksum mismatch -> corrupted copies that were replaced
  $G api transfer GET /task/$t/successful_transfers -Q limit=1000 2>/dev/null | python3 -c '
import json,sys
d=json.load(sys.stdin); print("\n".join(x["destination_path"] for x in d.get("DATA",[])))' > $HERE/logs/verify_redownloaded_$disk.txt
  say "verify $disk task $t: $st  (files, re-downloaded, skipped-identical, faults) = $info"
done
say "VERIFY DONE"
