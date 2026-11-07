#!/bin/bash
# igm_updater.sh — self-contained monitor updater. Runs on astro4 under
# nohup, so the igmprogress page keeps updating even when every laptop
# connection is gone. Uploads via curl sftp with ~/.netrc credentials
# (verified host key in ~/.ssh/known_hosts; no password appears here).
#   nohup bash igm_updater.sh > ~/igm_updater.log 2>&1 &
set -u
M=~/github/projects_desi/programs_igmdust/monitor
SC=/tmp/igmstatus
DEST=sftp://astrophysics.physics.fsu.edu/home/suzuki/public_html/igmprogress
mkdir -p "$SC"
for i in $(seq 1 500); do                       # ~25 h, relaunchable
  python3 "$M/igm_status.py" > "$SC/base.json" 2>/dev/null
  { cat /proc/loadavg; nproc; free -b | awk 'NR==2{print $2, $7}'; } \
    > "$SC/astro4.stat" 2>/dev/null
  python3 "$M/merge_status.py" "$SC/base.json" astro4 "$SC/astro4.stat" \
    > /dev/null 2>&1
  [ -s "$SC/status.json" ] && curl -sn -m 30 -T "$SC/status.json" \
    "$DEST/status.json" || echo "upload failed $(date)"
  sleep 180
done
echo "UPDATER_ENDED $(date)"
