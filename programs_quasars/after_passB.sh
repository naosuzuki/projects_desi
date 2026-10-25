#!/bin/zsh
# After pass B completes (34508 files): approved cleanup, then pass Q DESI with 4 shards (USB 3 link). Detached.
LOG=~/data/sdssdesi_status/after_passB.log
TOTAL=34508
while true; do
  n=$(ls ~/data/sdssdesi_passB/*/ 2>/dev/null | grep -c parquet)
  if [ "$n" -ge "$TOTAL" ] && ! pgrep -f 'python3 02_passB_measure.py --all' >/dev/null; then break; fi
  sleep 120
done
echo "$(date '+%F %T') pass B complete ($n files)" >> $LOG
rm -rf ~/data/sdssdesi_passB_test ~/data/sdssdesi_passB_test2 ~/data/sdssdesi_products/zlookup_main.parquet
echo "$(date '+%F %T') cleanup done (passB_test, passB_test2, zlookup_main.parquet removed)" >> $LOG
cd ~/github/projects_desi/programs_quasars
for i in 0 1 2 3; do
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --with-stars --expfibermap --shard $i/4 >> ~/data/sdssdesi_status/passQ_desi$i.log 2>&1 &
  sleep 3
done
for p in $(pgrep -f 'python3 02_passQ_desi.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
echo "$(date '+%F %T') pass Q DESI started with 4 shards" >> $LOG
