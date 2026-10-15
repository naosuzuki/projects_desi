#!/bin/zsh
# 05_overflow_to_exdisk1.sh — send the DESI DR2 coadds that did not fit on exdisk2 to exdisk1 (Nao, 2026-10-15:
# "I created a space on exdisk1"; exdisk2 = 24.0 TB, DR2 main dark+bright ~24.03 TB).
#
#   zsh 05_overflow_to_exdisk1.sh            # waits until 04_disk_guard.sh has cancelled the exdisk2 task
#
# 1. waits for "cancelled" in logs/disk_guard.log (or for the exdisk2 task to stop being ACTIVE)
# 2. every manifest file that is absent on exdisk2, or present with a size different from NERSC (a file cut by the
#    cancel), goes into manifests/overflow_exdisk1.batch, destination /Volumes/My Book/desiredux/... (same tree);
#    cut files on exdisk2 are removed so only one, complete copy exists
# 3. submits one Globus task (checksum sync + verify) and writes the list to
#    desiredux/DR2_overflow_on_exdisk1.txt on both disks
HERE=${0:A:h}
G=~/.local/bin/globus
TASK=07cf18f4-c48e-11f1-b8f0-0affd5e180af
SRC=9d6d994a-6d04-11e5-ba46-22000b92c6ec
DST=$(cat ~/.globusonline/lta/client-id.txt)
LOG=$HERE/logs/overflow_exdisk1.log
say() { print -r -- "$(date '+%F %T') $*" | tee -a $LOG }

say "waiting for the exdisk2 task to stop"
until grep -q "cancelled" $HERE/logs/disk_guard.log 2>/dev/null; do
  st=$($G task show $TASK -F json 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["status"])' 2>/dev/null)
  [ -n "$st" ] && [ "$st" != ACTIVE ] && break
  sleep 120
done
sleep 60                                         # let the endpoint close the last files
say "exdisk2 task stopped; building the overflow list"

python3 - "$HERE" <<'EOF' | tee -a $LOG
import os, sys, json, subprocess
here = sys.argv[1]
rows = [l.split() for l in open(f'{here}/manifests/healpix_main_darkbright.batch')]
known = {}
for f in ('logs/missing_now_sizes.txt',):
    for l in open(f'{here}/{f}'):
        a = l.split(); known[a[1]] = int(a[0])
G = os.path.expanduser('~/.local/bin/globus'); SRC = '9d6d994a-6d04-11e5-ba46-22000b92c6ec'
def nersc_size(src):
    d, b = os.path.split(src)
    out = subprocess.run([G, 'ls', '-l', '-F', 'json', f'{SRC}:{d}/'], capture_output=True, text=True).stdout
    return [x['size'] for x in json.loads(out)['DATA'] if x['name'] == b][0]
over, cut = [], []
for src, dst in rows:
    if not os.path.exists(dst):
        over.append((src, dst)); continue
    if src in known:                              # was missing at 08:14: check it arrived complete
        if os.path.getsize(dst) != known[src]:
            cut.append(dst); over.append((src, dst))
for dst in cut:
    os.remove(dst)
new = '/Volumes/My Book/desiredux'
with open(f'{here}/manifests/overflow_exdisk1.batch', 'w') as fb:
    for src, dst in over:
        fb.write(f'"{src}" "{dst.replace("/Volumes/exdisk2/desiredux", new)}"\n')   # quoted: "My Book" has a space
lst = '\n'.join(d.replace('/Volumes/exdisk2/desiredux/', '') for s, d in over) + '\n'
for root in ('/Volumes/exdisk2/desiredux', new):
    open(f'{root}/DR2_overflow_on_exdisk1.txt', 'w').write(
        f'# DESI DR2 loa healpix coadds stored on exdisk1 ({new}) because exdisk2 was full, 2026-10-15\n' + lst)
tot = sum(known.get(s, 0) for s, d in over)
print(f'overflow: {len(over)} files (~{tot/1e9:.0f} GB known sizes); removed {len(cut)} cut files on exdisk2')
EOF

n=$(wc -l < $HERE/manifests/overflow_exdisk1.batch)
[ $n -gt 0 ] || { say "nothing to send"; exit 0; }
T2=$($G transfer --batch $HERE/manifests/overflow_exdisk1.batch --label "DESI DR2 loa overflow to exdisk1" \
     --sync-level checksum --preserve-timestamp --verify-checksum --jmespath 'task_id' --format unix $SRC $DST)
say "submitted overflow task $T2 ($n files) -> /Volumes/My Book/desiredux"
