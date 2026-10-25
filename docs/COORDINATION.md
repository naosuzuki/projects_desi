# Coordination between the three template threads (Galaxy, Quasar, Stellar)

Written 2026-10-25 17:20 by the Quasar thread at Nao's request. One machine, one coadd disk, one status directory:
these are the rules so that the three sessions do not step on each other. Edit this file when you change a shared
thing; say so in your commit message. Sessions: **Galaxy Template**, **Quasar Template**, **Stellar Template**
(sidebar group projects_desi).

## 1. Shared resources and who owns them
| resource | owner | rule |
|---|---|---|
| **pass B** (`programs_templates/02_passB_measure.py`, outputs `~/data/sdssdesi_passB/`) | Galaxy | the one pass that reads all 34,508 coadds; stars and quasars consume its outputs. **Only the Galaxy thread (or its autoscaler/resume scripts) starts, stops or re-shards it.** Others: hands off. Status 2026-10-25 17:12: restarted with 4 workers on the now USB 3 link (5 Gb/s), 27.5k/34.5k done, ≈25–40 files/min expected |
| **exdisk2 / exdisk1** (coadds, read-only) | all | the bottleneck. Queue order for passes that read the coadds: **1. pass B → 2. pass Q DESI (quasars, 4 shards, auto-started by `~/data/sdssdesi_status/after_passB.sh` when pass B completes) → 3. star passes that need the coadds** (Stellar thread: say in this file if you need the disk before 2 finishes). Never write to the external disks |
| `~/data/sdssdesi_products/` (was `sdssdesi_cache`, symlink kept) | all | rebuildable-but-expensive products: zall Parquet, lookups, spAll/specObj extracts; per-track subdirs `quasar/`, `stars/`; galaxy products go in `galaxy/` |
| `~/data/sdssdesi_status/` | all | one `status.json` (written by the monitor), run logs, `step_<name>.json` progress files (protocol in `programs_monitor/README.md`). The Stellar thread's `~/data/sdssdesi_status_stars/` is being merged into it (commit 089cee7) |
| **monitor** `programs_monitor/monitor.py` on port 8765 | Quasar (started it); the Stellar thread merged its page in | **one instance with the guard**; restart it (`pkill -f 'monitor.py --port 8765'`, relaunch) only when no worker is paused (`n_paused` in status.json); if you add a step, add it to `STEPS` there and say so here |
| `programs_templates/` (pass B, zlookup, galaxy grid) | Galaxy | the Quasar thread imports `02_passB_measure.py` (filters, rebin, extinction, galaxy line list) — **changing function signatures or the LINES table there changes pass Q too**: announce here first |
| `programs_quasars/`, `programs_stars/`, `programs_galaxies/` (to be created) | each thread | nobody edits another thread's directory |
| CPU (14 cores, 64 GB) | all | pass B ≈ 1 core/worker; keep ≥ 16 GB free (guard pauses at 8 GB); nice ≥ 10 for every batch job; anything > 4 cores for more than an hour: note it here |
| **git** (one repo, branch main) | all | commit only your own directory + docs you own; `git pull --rebase` before pushing; never `git add -A` (large files) |

## 2. Cross-thread facts everyone should know
- **Air vs vacuum wavelengths (found 2026-10-24 by the Quasar thread):** the optical line list `LINES` in
  `02_passB_measure.py` (Hα, Hβ, [O III], [O II], [N II], [S II], Hγ, Hδ, [Ne III]) is in AIR wavelengths; SDSS and DESI
  spectra are in VACUUM. Narrow lines are fitted ≈80 km/s off centre → flux loss with the 80 km/s Gaussian.
  Verified with DR16Q ([O III] peak offset +76 km/s → −5 km/s after the fix). Quasars: fixed in `qsolib.py` by
  overriding the table. **Galaxy thread: decide whether to fix the table now (outputs then differ between files done
  before and after) or re-measure the lines in the cell-assignment pass; stars are unaffected.**
- **Step 1 reference for SDSS point sources is the PSF magnitude**, not the fiber magnitude (the pipeline calibrates to
  PSF flux: synthetic − fiber = −0.70 mag BOSS, −0.30 SDSS-I/II). Relevant to the Stellar thread's SDSS pass.
- **Spectrophotometric tilt:** pass B stores the per-object two-band tilt (CAL_B_EST) but applies only the scale; the
  tilt is to come from the stars per exposure/fiber/petal (Stellar thread, calibration 2) and then applied to
  galaxies and quasars. Quasars must not get an imaging-epoch tilt (variability).
- **pass B outputs carry every class**, so the Stellar thread gets CAL_*, S/N, synthetic mags for all DESI stars when
  pass B completes; the Quasar thread gets the pass-B columns through pass Q.
- **DESI S/N prediction** (for counts before a pass is done): log S/N = −0.449 + 0.846 log FIBERFLUX_R + 0.430 log
  TSNR2_QSO (quasars; fit on the pass-B test coadds, rms 0.13 dex). Galaxies/stars: refit on their own class.
- **Systemic redshifts for quasars** come from line peaks ([O III], Mg II), not the pipeline z (Nao's rule); galaxies
  keep the pipeline z.
- **Grids:** galaxies ms.tex Table 3 (5 axes, ≈7,000 stacks); quasars ms.tex tab:qsogrid (6 axes; SDSS alone gives 537
  stacks at floor 300, DESI will add ≈3×); stars HR cells (Stellar thread).

## 3. Timeline (2026-10-25)
| when | what | thread |
|---|---|---|
| now → ~21:00 | pass B finishes (4 workers, USB 3) | Galaxy |
| at pass B end | `after_passB.sh`: deletes passB_test dirs + zlookup_main.parquet (Nao approved), starts pass Q DESI (4 shards, ≈4–6 h on USB 3) | Quasar (automatic) |
| after pass Q DESI | DESI master table, cells, stacks (SDSS stacks can start earlier on the internal disk) | Quasar |
| after pass B | galaxy cell assignment, pass C (reads coadds again → coordinate with pass Q DESI: run after it, or merge the two reads) | Galaxy |
| after pass B | star calibration from pass-B outputs; coadd-reading star passes after pass Q DESI | Stellar |

## 3b. Agreed 2026-10-25 17:30 (Galaxy ↔ Quasar)
- Pass B: Galaxy thread tests 6 workers (keeps them if > 23.4 files/min); `02_passB_measure.py` stays untouched until
  pass B completes. Then the LINES table goes to vacuum (VAC dict in `programs_quasars/qsolib.py`), documented in
  ms.tex + README; all galaxy lines are re-measured in pass C (reads every coadd anyway). Interim cell assignment
  uses the air-wavelength lines (bias ≤ 36 % flux for the narrowest fits; factor-2 EW bins tolerate it).
- Disk after pass B: pass Q DESI first (auto, `after_passB.sh`: triggers only when 34,508 parquet files exist AND no
  pass-B worker is alive — a worker dying early does not trigger it), then galaxy pass C (its cell table is catalog
  work in the meantime). Pass C's per-spectrum stage will be a plug-in function so a later single pass can serve both
  tracks; this round stays sequential.
- Memory: pass Q DESI = 4 shards × ≈4 GB; the guard pauses below 8 GB free.

## 3c. Agreed 2026-10-25 17:50 (Stellar ↔ Quasar): one DESI walk for quasars and stars
`programs_quasars/02_passQ_desi.py --with-stars` reads each coadd once and runs the quasar chain on the QSO rows and
`programs_stars/step2_3_pass_desi.measure_star_chunk` on the STAR rows (outputs `~/data/sdssdesi_passS/`, resumable
per track). Tested on one coadd (112 QSO + 12 STAR rows). `after_passB.sh` starts it this way (4 shards) when pass B
completes. The Stellar thread owns `measure_star_chunk` / `stars_for` / `load_star_lookup` (interface frozen: same
signature as pass B's measure_chunk); the Quasar thread owns the walk. Galaxy pass C follows, as in 3b.

## 4. Open questions for Nao (cross-thread)
1. Pass C (galaxy stacks) and pass Q DESI both re-read all coadds (≈1 day each on USB 3). Should one pass read each
   coadd once and serve both (a "pass C+Q" that writes rest-frame spectra for galaxies and quasars)? Recommended:
   yes — a 1–2 TB SSD for rest-frame spectra was already the plan in programs_templates/README.md.
2. Who fixes the air-wavelength table in pass B and when (see §2).
