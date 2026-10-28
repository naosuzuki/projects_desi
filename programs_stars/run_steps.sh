#!/bin/zsh
# Stars thread: the steps in order (Nao, 2026-10-24). Each step is resumable; passes run detached and niced.
# Usage: ./run_steps.sh <step>    steps: 1 | 2-3sdss | 2b | 3b | 2-3desi | 5 | 6 | chain | desi | desi67 | monitor
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
  desi)     # after the DESI walk (2026-10-27): 2b (if not already running/done) -> 5 -> 6 (capped, one read per coadd) -> 7, one after another
            while pgrep -f "[s]tep2b_spcal.py --survey desi" >/dev/null; do sleep 30; done
            [ -f ~/data/sdssdesi_cache/stars/spcal_desi.parquet ] || python3 step2b_spcal.py --survey desi > $S/step2b_desi.txt 2>&1
            python3 step5_hr_cells_check.py --survey desi --nmin 100 --snmin 5 > $S/step5_desi.txt 2>&1
            nice -n 10 python3 step6_stack.py --survey desi --tag desi_0.5x0.1_gspphot --nmin 100 --snmin 5 --maxper ${MAXPER:-300} --onepass > $S/step6_desi.log 2>&1
            nice -n 10 python3 step7_empca.py --survey desi --tag desi_0.5x0.1_gspphot --individuals 5 --plot > $S/step7_desi.txt 2>&1
            echo "DESI chain finished $(date '+%F %T')" ;;
  desi67)   # step 6 from star_hook's store (filled inside the quasar-stack read + fill_store.py) -> step 7. Starts when the
            # 300-per-stack coadds (first 18,875 lines of star_priority_files.txt) are all stored, or at DEADLINE (HHMM, default 0715)
            # with whatever is stored (--store-only: no exdisk2 reads); Nao wants the stars done before 08:30 (2026-10-28)
            P=~/data/sdssdesi_products/stars/star_priority_files.txt; ST=/Volumes/exdisk3/sdssdesi_star_spectra/desi
            until [ $(date +%H%M) -ge ${DEADLINE:-0715} ] && [ $(date +%H%M) -lt 1200 ] || \
                  python3 store_ready.py 18875 >/dev/null; do sleep 120; done
            echo "step 6 starts $(date '+%F %T')"
            nice -n 5 python3 step6_stack.py --survey desi --tag desi_0.5x0.1_gspphot --nmin 100 --snmin 5 --maxper ${MAXPER:-300} --from-store --store-only --onepass > $S/step6_desi.log 2>&1
            nice -n 5 python3 step7_empca.py --survey desi --tag desi_0.5x0.1_gspphot --individuals 5 --plot > $S/step7_desi.txt 2>&1
            echo "DESI steps 6-7 finished $(date '+%F %T')" ;;
  desi-early) # first DESI stacks + EM-PCA with 100 stars per stack as soon as the tier-1 coadds are stored (desi67 redoes them with 300)
            until python3 store_ready.py 11288 >/dev/null; do sleep 120; done
            echo "early step 6 (maxper 100) starts $(date '+%F %T')"
            nice -n 5 python3 step6_stack.py --survey desi --tag desi_0.5x0.1_gspphot --nmin 100 --snmin 5 --maxper 100 --from-store --store-only --onepass > $S/step6_desi_early.log 2>&1
            nice -n 5 python3 step7_empca.py --survey desi --tag desi_0.5x0.1_gspphot --individuals 5 --plot > $S/step7_desi_early.txt 2>&1
            echo "early DESI steps 6-7 finished $(date '+%F %T')" ;;
  monitor)  echo "the stars track is on the main monitor: http://127.0.0.1:8765/monitor.html (programs_monitor/monitor.py)" ;;
  *)        echo "usage: $0 1 | 2-3sdss | 2b | 3b | 2-3desi | 5 | 6 | chain | desi | desi67 | monitor" ;;
esac
