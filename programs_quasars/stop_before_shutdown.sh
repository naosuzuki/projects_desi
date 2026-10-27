#!/bin/zsh
# Graceful stop of every SDSS-DESI batch job before a planned shutdown (Nao: shutdown 2026-10-26 ~04:30).
# Scheduled: this script sleeps until the target time, then stops the workers (SIGTERM), the keep-awake holds and the
# monitor, and writes shutdown.log. Every pass writes each output file as .tmp + rename, so a stopped pass loses at most
# the files in progress; resume_after_reboot.sh relaunches everything (existing outputs are skipped).
TARGET="${1:-2026-10-26 04:15:00}"
LOG=~/data/sdssdesi_status/shutdown.log
now=$(date +%s); tgt=$(date -j -f '%Y-%m-%d %H:%M:%S' "$TARGET" +%s)
echo "$(date '+%F %T') scheduled stop at $TARGET ($(( (tgt-now)/60 )) min from now)" >> $LOG
[ $tgt -gt $now ] && sleep $((tgt-now))
echo "$(date '+%F %T') stopping: walk files done Q=$(ls ~/data/sdssdesi_passQ/*/ 2>/dev/null | grep -c parquet) S=$(ls ~/data/sdssdesi_passS/*/ 2>/dev/null | grep -c parquet) E=$(ls ~/data/sdssdesi_expfibermap/*/ 2>/dev/null | grep -c parquet)" >> $LOG
pkill -f 'after_passB.sh' 2>/dev/null
for pat in 'python3 02_passQ_desi.py --all' 'python3 02_passB_measure.py --all' 'python3 step2_3_pass_desi.py' 'python3 08_passC_stack.py --all' 'python3 06_passB_summary'; do
  for p in $(pgrep -f "$pat"); do kill -TERM $p 2>/dev/null && echo "$(date '+%F %T') TERM $p ($pat)" >> $LOG; done
done
for i in $(seq 1 24); do pgrep -f 'python3 02_passQ_desi.py --all' >/dev/null || break; sleep 5; done   # up to 2 min: the walk finishes its coadd and writes the galaxy checkpoint
for pat in 'python3 02_passQ_desi.py --all' 'python3 02_passB_measure.py --all'; do for p in $(pgrep -f "$pat"); do kill -KILL $p 2>/dev/null && echo "$(date '+%F %T') KILL $p" >> $LOG; done; done
pkill -x caffeinate 2>/dev/null
pkill -f 'monitor.py --port' 2>/dev/null
rm -f ~/data/sdssdesi_passQ/*/*.tmp ~/data/sdssdesi_passS/*/*.tmp ~/data/sdssdesi_expfibermap/*/*.tmp 2>/dev/null
echo "$(date '+%F %T') all stopped; safe to shut down. After reboot: ~/data/sdssdesi_status/resume_after_reboot.sh" >> $LOG
