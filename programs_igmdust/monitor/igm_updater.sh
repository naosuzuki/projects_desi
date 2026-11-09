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
  # astrofs4's own collector (astrofs4_stat.sh) drops its stats in the
  # shared home; include them only while fresh (< 10 min old)
  HOSTARGS=(astro4 "$SC/astro4.stat")
  FS4=~/.igmstats/astrofs4.stat
  if [ -f "$FS4" ] && [ -n "$(find "$FS4" -mmin -10 2>/dev/null)" ]; then
    HOSTARGS+=(astrofs4 "$FS4")
  fi
  python3 "$M/merge_status.py" "$SC/base.json" "${HOSTARGS[@]}" \
    > /dev/null 2>&1
  [ -s "$SC/status.json" ] && curl -sn -m 30 -T "$SC/status.json" \
    "$DEST/status.json" || echo "upload failed $(date)"
  sleep 180
done
echo "UPDATER_ENDED $(date)"
