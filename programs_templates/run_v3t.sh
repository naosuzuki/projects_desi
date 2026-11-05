#!/bin/zsh
# Galaxy release v3t (Nao 2026-11-02): v3 re-read WITH the telluric correction (tellhook.py -> programs_telluric/tellapply.py, per spectrum,
# before stacking). DESI coadds on exdisk2 (read-only) -> ~/data/sdssdesi_passC_v3t (symlink to exdisk3), 4 shards; SDSS spPlates -> ~/data/sdssdesi_passG_sdss_v3t.
# Usage: ./run_v3t.sh desi 0 1 2 3      start DESI shards (resumable; keep N=4)
#        ./run_v3t.sh sdss dr8 dr17     start the SDSS surveys (one shard each)
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status
for d in ~/data/sdssdesi_passC_v3t ~/data/sdssdesi_passG_sdss_v3t; do
  if [ ! -e $d ]; then mkdir -p /Volumes/exdisk3/$(basename $d) && ln -s /Volumes/exdisk3/$(basename $d) $d; fi
done
what=$1; shift
for i in "$@"; do
  if [ "$what" = desi ]; then
    pat="08_passC_stack.py --all --v3 --telluric --shard $i/4"; log=$S/gal_v3t_desi_$i.log
  else
    pat="11_passG_sdss.py --v3 --telluric --survey $i"; log=$S/gal_v3t_sdss_$i.log
  fi
  if pgrep -f "$pat" >/dev/null; then echo "$pat running"; continue; fi
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 ${=pat} >> $log 2>&1 &
  nohup caffeinate -i -w $! >/dev/null 2>&1 &
  echo "$(date '+%F %T') v3t $what $i started (pid $!)" | tee -a $S/gal_v3t.log
  sleep 2
done
