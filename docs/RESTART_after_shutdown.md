# Restart checklist after the 2026-10-26 ~04:30 shutdown

State at shutdown (written 2026-10-25 ~22:15 by the Galaxy Template session; the pre-shutdown script
`programs_templates/preshutdown.sh` runs at 04:15 and appends the final counts to `~/data/sdssdesi_status/shutdown.log`).

## What is finished
- **Pass B** (DESI per-spectrum measurements): COMPLETE — 34,508/34,508 coadds, 53.9 M spectra, 36 GB in
  `~/data/sdssdesi_passB/{bright,dark}/hpx<N>.parquet`. Summary: `~/data/sdssdesi_products/passB_summary.json`,
  compact columns `passB_summary_columns.parquet`, galaxy cell assignment `galaxy_cells.parquet` + `galaxy_cell_table.csv`
  (5,762 cells at floor 300; 6,850 at floor 200).
- SDSS star pass (Stellar thread, internal disk): done. Catalog caches in `~/data/sdssdesi_products/`.

## What was interrupted
- **Pass Q DESI** (Quasar thread's walk: QSO + star rows + EXP_FIBERMAP dump; `programs_quasars/02_passQ_desi.py
  --all --with-stars --expfibermap --shard i/4`): started 2026-10-25 21:51, ETA was Mon 12:45. Per-file atomic
  outputs in `~/data/sdssdesi_passQ`, `~/data/sdssdesi_passS`, `~/data/sdssdesi_expfibermap`; existing files are
  skipped on restart. Resume command: see `docs/COORDINATION.md` (Quasar thread) — or the Quasar Template session.
- Monitor (`programs_monitor/monitor.py --port 8765`) and the keep-awake holds were stopped by the script.

## One-command restart — use this one (it refuses to start on a USB 2 link)
```
zsh ~/data/sdssdesi_status/post_reboot.sh            # Stellar thread: aborts unless exdisk2 reports "Gb/s" and the
                                                      # coadd tree is readable; FORCE=1 to override; then monitor + walk
```
(`resume_after_reboot.sh` from the Quasar thread does the same restart; use it only if post_reboot.sh is missing.)
Stop scripts that ran before the shutdown (all idempotent): `pre_shutdown.sh` 04:10 (Stellar), `stop_before_shutdown.sh`
04:15 (Quasar), `preshutdown.sh` 04:18 (Galaxy). Their logs: `~/data/sdssdesi_status/shutdown.log`, `shutdown_state.txt`.

## Restart order (manual)
1. Power on; wait for both external disks to mount: `mount | grep exdisk` must show `/Volumes/exdisk1` (exFAT) and
   `/Volumes/exdisk2` (HFS+). If exdisk2 is missing, reseat its USB-3 cable (the new one).
2. **Check the USB link of exdisk2** — it must say 5 Gb/s, not 480 Mb/s (480 = USB 2 = 3x slower):
   `system_profiler SPUSBDataType | awk 'BEGIN{RS="\n\n"} /BSD Name: disk5/ {print}' | grep Speed`
   (the BSD number may change after a reboot; look for the 24 TB "Apple_HFS exdisk2" in `diskutil list external`).
3. Start the monitor: `cd ~/github/projects_desi && nohup python3 programs_monitor/monitor.py --port 8765 >> ~/data/sdssdesi_status/monitor.log 2>&1 &`
   and open http://127.0.0.1:8765/monitor.html
4. Resume pass Q DESI (Quasar thread's command in COORDINATION.md; 4 shards, nice 10, OPENBLAS_NUM_THREADS=1) and
   re-attach keep-awake: `for p in $(pgrep -f 'python3 02_passQ_desi.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done`
5. After pass Q finishes: Galaxy pass C (stacking) — written and tested (programs_templates/README.md):
   `cd ~/github/projects_desi/programs_templates && ./run_passC.sh 4 /Volumes/exdisk3/sdssdesi_restframe`
   (refuses to start while pass Q runs or on a USB 2 link; exdisk3 store approved by Nao 2026-10-26),
   then `python3 09_merge_stacks.py --floor 200` and `python3 10_galaxy_empca.py --tag v1`.
   The Stellar thread's `response_coadd_desi.parquet` (Calibration 2) is not used yet; it will enter the next iteration.

## Rules unchanged
Never write to the external disks; keep >= 16 GB RAM free; one external-disk walk at a time
(pass Q, then pass C); ms.tex is the documentation (Nao commits); scripts numbered in programs_templates/.
