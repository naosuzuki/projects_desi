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
| **stars products** `~/data/sdssdesi_products/stars/` (cells, members, spcal, zp), `~/data/sdssdesi_products/desi_stars_ebv3d.parquet`, `~/data/sdssdesi_products/sdss_stars_gaia.parquet`, `~/data/sdssdesi_passS_sdss/` (SDSS pass), `~/data/sdssdesi_passS/` (DESI pass, after pass B), `~/data/sdssdesi_stars/` (stacks) | Stellar | read freely; progress files `step_star_*.json`; docs `docs/HANDOVER_stars.md`, `programs_stars/README.md` |
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
| after pass B | Stellar: catalog-level work only. **The DESI star pass rides inside pass Q DESI** (§3c, hook ready) so the coadds are read once; if declined it runs after pass Q DESI (≈ the same 4–6 h again). Galaxy pass C after both | Stellar |

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

## 3c. Stellar ↔ Quasar: one DESI walk for quasars AND stars (Stellar thread, 2026-10-25 17:45)
The star measurements need exactly the arrays pass Q reads (FIBERMAP + combined B/R/Z flux/ivar/mask). Hook, ready and tested:
`programs_stars/step2_3_pass_desi.write_star_file(path, fm, w, flux, ivar, mask, percam, zsub, lookup, program)` with
`lookup = step2_3_pass_desi.load_star_lookup()` loaded once per worker (12 M rows, ≈1 GB RSS). In `02_passQ_desi.py` this is a
`--with-stars` flag and three lines after `combine_cameras`; the FIBERMAP columns it needs (`TARGET_RA, TARGET_DEC, EBV, PHOTSYS,
FLUX_G/R/Z, FIBERFLUX_G/R/Z, MEAN_PSF_TO_FIBER_SPECFLUX, COADD_FIBERSTATUS`) are the ones pass Q already reads. Cost per coadd ≈3 s
of measurement for ≈400 stars (the read dominates); memory +≈1.5 GB per shard. Output `~/data/sdssdesi_passS/{dark,bright}/hpx<N>.parquet`,
resumable, independent of the quasar output. **Quasar thread: yes/no?** On yes, add the flag yourself or ask for the exact patch
(the Stellar thread does not edit programs_quasars/). On no, the Stellar thread runs `step2_3_pass_desi.py --all --shard i/4`
after pass Q DESI. EXP_FIBERMAP (EXPID/NIGHT/TILEID/PETAL_LOC/FIBER per spectrum, wanted by all three threads for Calibration 2)
is dumped by `step2_3_pass_desi.write_expfib_file(f, path, program)` (open fitsio handle; all rows; columns TARGETID, NIGHT, EXPID, MJD, TILEID,
EXPTIME, PETAL_LOC, DEVICE_LOC, LOCATION, FIBER, FIBERSTATUS, FIBER_X/Y, DELTA_X/Y, PSF_TO_FIBER_SPECFLUX, IN_COADD_B/R/Z + PROGRAM, HEALPIX)
to `~/data/sdssdesi_expfibermap/{dark,bright}/hpx<N>.parquet` — written 2026-10-25, tested; call it in the same walk (one line, before
the spectra are read). Agreed with the Galaxy thread (pass C joins it by TARGETID, PROGRAM, HEALPIX).

Quasar-thread addendum (17:55): implemented in `02_passQ_desi.py --with-stars` through `S.write_star_file(...)` after the
single `combine_cameras` read; plus `--expfibermap`: the walk dumps EXP_FIBERMAP (TARGETID, NIGHT, EXPID, MJD, TILEID,
EXPTIME, PETAL_LOC, DEVICE_LOC, LOCATION, FIBER, FIBERSTATUS, FIBER_X/Y, DELTA_X/Y, FIBER_RA/DEC, PSF_TO_FIBER_SPECFLUX, IN_COADD_B/R/Z) for ALL rows of
each coadd to `~/data/sdssdesi_expfibermap/<program>/hpx<N>.parquet` (≈1.4 exposures per spectrum → ≈90 M rows, a few
GB). `after_passB.sh` runs the walk with both flags (4 shards). **Disk warning (17:55): the internal SSD has 60 GB free
(97 % used)**; pass B (+7 GB to go), pass Q/S (+≈8 GB) and the EXP_FIBERMAP dump (+≈4 GB) fit, but nothing else large
should land on it until Nao frees space (candidates listed in HANDOVER_quasars 8.23).

## 3d. Stellar → Galaxy/Quasar: the Calibration-2 interface (agreed with the Galaxy thread's request, 2026-10-25)
What the Stellar thread will publish, all under `~/data/sdssdesi_products/stars/`:
1. `zp_desi.csv` — Legacy Surveys (DECam S, BASS/MzLS N) zero points and color terms relative to Gaia G/BP/RP, from the joint fit
   (step 2b) on all DESI stars. Needs the synthetic Gaia magnitudes, i.e. the DESI star pass (spectra), so it comes right after the
   walk, not from pass-B rows alone. One-coadd preview: g −0.07, r −0.06, z −0.02 (not for use yet). Pass C applies them as global
   offsets on top of the fiber-flux scale; ms.tex Step 1 cites them.
2. `response_exp_petal.parquet` — Calibration 2: per (EXPID, PETAL_LOC) [and per FIBER where the stars allow] the multiplicative
   response c(λ) = A + B·log10(λ/6000) (+ a smooth residual curve on the log grid when needed), with N_STARS and scatter. Built
   from the DESI star pass residuals against the HR-cell stacks, so it needs the EXP_FIBERMAP mapping (§3c).
3. `response_coadd_desi.parquet` — the same collapsed per coadd spectrum, keyed by (TARGETID, PROGRAM): ivar-weighted mean of the
   exposure/petal solutions of that spectrum's exposures. This is the table pass C can apply directly; (2) is for anyone who
   re-reads exposures.
Interim (now): per-star (A_MAG, TILT) for SDSS stars in `spcal_sdss.parquet`; the DESI equivalent after the walk. These are per
star, not per fiber, and are not to be applied to galaxies.
Stellar thread CPU use: nothing > 4 cores; nothing reads exdisk2 tonight. Air/vacuum: the stellar index bands (Lick/Cenarro/Reid,
air) are converted to vacuum in `starlib.py` since 2026-10-25; the SDSS star pass ran with air bands (≈1 Å shifts, harmless for the
homogeneity test; rerun with the next refinement).

## 4. Open questions for Nao (cross-thread)
1. Pass C (galaxy stacks) and pass Q DESI both re-read all coadds (≈1 day each on USB 3). Should one pass read each
   coadd once and serve both (a "pass C+Q" that writes rest-frame spectra for galaxies and quasars)? Recommended:
   yes — a 1–2 TB SSD for rest-frame spectra was already the plan in programs_templates/README.md.
2. Who fixes the air-wavelength table in pass B and when (see §2).

## 5. Shutdown 2026-10-26 ~04:30 (Nao) — one sequence, agreed by the three threads (22:10)
Stop (automatic; Nao only shuts down at 04:30):
- **04:10** `~/data/sdssdesi_status/pre_shutdown.sh` (Stellar): SIGTERM to every batch job (walk shards, pass B/S/Q, stacks,
  EM-PCA, table jobs), 5-min grace, kills stragglers, stops caffeinate and the monitor, writes `shutdown_state.txt`
  (file counts), removes `.tmp` files, **unmounts exdisk2 and exdisk1**.
- 04:15 `stop_before_shutdown.sh` (Quasar, pid 10770) and ~04:18 the Galaxy stop: safety nets, find nothing running,
  log to `shutdown.log`; their unmounts report "already unmounted" (harmless).
Resume (after the reboot, with exdisk2 on a **USB 3** port) — ONE command, run once, from any terminal:
```
zsh ~/data/sdssdesi_status/post_reboot.sh
```
(checks the disks and the 5 Gb/s link, restarts the monitor on 8765, resumes the walk `02_passQ_desi.py --all --with-stars
--expfibermap` with 4 shards + caffeinate; `resume_after_reboot.sh` is the equivalent Quasar-thread script, either works,
both refuse to start a second walk). Everything is per-file atomic (`.tmp` + rename) and skipped on resume: pass B
(complete), the walk's three outputs (a coadd is re-read when any of its three outputs is missing), the SDSS passes
(complete). Expected at 04:10: walk ≈ 14,000 of 34,508 coadds; ≈ 8–10 h remaining after resume. Galaxy pass C and star
jobs are started by their threads after the walk finishes.
