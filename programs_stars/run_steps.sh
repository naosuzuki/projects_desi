#!/bin/zsh
# Stars thread: the steps in order (Nao, 2026-10-24). Each step is resumable; passes run detached and niced.
# Usage: ./run_steps.sh <step>    steps: 1 | 2-3sdss | 2b | 3b | 2-3desi | 5 | 6 | chain | monitor
set -e
cd ~/github/projects_desi/programs_stars
S=~/data/sdssdesi_status; mkdir -p $S
case "$1" in
  1)        python3 step1_desi_catalog.py; python3 step1_desi_hr.py; python3 step1_sdss_catalog.py ;;
  2-3sdss)  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 step2_3_pass_sdss.py --survey all --all >> $S/passS_sdss.log 2>&1 &
            echo "SDSS pass started: tail -f $S/passS_sdss.log" ;;
  2b)       python3 step2b_spcal.py --survey sdss | tee $S/step2b_sdss.txt; [ -d ~/data/sdssdesi_passS ] && python3 step2b_spcal.py --survey desi | tee $S/step2b_desi.txt ;;
  chain)    nohup zsh -c 'while pgrep -f "[s]tep2_3_pass_sdss.py" >/dev/null; do sleep 30; done; python3 step2b_spcal.py --survey sdss > '$S'/step2b_sdss.txt 2>&1; python3 step5_hr_cells_check.py --survey sdss --nmin 100 --snmin 5 > '$S'/step5_sdss.txt 2>&1; python3 step6_stack.py --survey sdss --tag sdss_0.5x0.1_gspphot --nmin 100 --snmin 5 > '$S'/step6_sdss.log 2>&1' > $S/chain_sdss.log 2>&1 &
            echo "chain started (waits for the SDSS pass, then 2b -> 5 -> 6); logs in $S/step2b_sdss.txt step5_sdss.txt step6_sdss.log" ;;
  3b)       nohup nice -n 10 python3 step3b_ebv3d.py >> $S/ebv3d.log 2>&1 & echo "3D extinction started: tail -f $S/ebv3d.log" ;;
  2-3desi)  # only after pass B has finished (external disks); two workers like pass B / pass Q
            if pgrep -f '02_passB_measure.py --all' >/dev/null; then echo "pass B still running - wait"; exit 1; fi
            for i in 0 1; do OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 step2_3_pass_desi.py --all --shard $i/2 >> $S/passS_desi$i.log 2>&1 & done
            for p in $(pgrep -f 'step2_3_pass_desi.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
            echo "DESI pass started (2 workers)" ;;
  5)        python3 step5_hr_cells_check.py --survey sdss --nmin 100 --snmin 5 | tee $S/step5_sdss.txt
            [ -d ~/data/sdssdesi_passS ] && python3 step5_hr_cells_check.py --survey desi --nmin 100 --snmin 5 | tee $S/step5_desi.txt ;;
  6)        nohup nice -n 10 python3 step6_stack.py --survey sdss --tag sdss_0.5x0.1_gspphot --nmin 100 --snmin 5 >> $S/step6_sdss.log 2>&1 &
            echo "SDSS stacks started: tail -f $S/step6_sdss.log" ;;
  monitor)  echo "the stars track is on the main monitor: http://127.0.0.1:8765/monitor.html (programs_monitor/monitor.py)" ;;
  *)        echo "usage: $0 1 | 2-3sdss | 2b | 3b | 2-3desi | 5 | 6 | chain | monitor" ;;
esac
