# programs_templates — galaxy stacking recipe (runbook)

The recipe is documented in `~/github/papers/26_sdssdesi_spectra/ms.tex`, §3.9
("Processing recipe for the galaxy stacks"). This file holds the operational side.
Status 2026-10-24: recipe written; **pass B script tested on 20 random coadds (28,772 spectra), all checks pass;
production run not started.** Scripts: `02_passB_measure.py` (pass B), `03_passB_check.py` (checks).
Test results: Step 1 convention validated on stars (catalog − synthetic = +0.03 r / +0.05 g, σ = 0.05–0.06; after
the scale correction the independent z band agrees to +0.006 ± 0.053 for stars); per-object tilt NOT applied
(its scatter, 0.04 mag in g−r for stars, is independent of spectrum S/N → photometric or per-fiber, to be solved
per exposure/petal from stars); stored as CAL_B_EST. 3.9 s/file → ~37 h single worker, ~20 h with --shard 0/2, 1/2.
Rebin agrees with astrolib to 0.6% rms (differs only at masked pixels, which astrolib includes).

## Order of operations
0. Sample: DESI DR2 loa main dark+bright coadds; SDSS DR17/DR8 spPlates. Cuts: ZCAT_PRIMARY,
   ZWARN=0, COADD_FIBERSTATUS=0, SPECTYPE=GALAXY, DELTACHI2>25, 0<z<1.6.
1. Spectrophotometry: synthetic LS g,r (z as check) vs FIBERFLUX_G/R; DESI spectra first
   multiplied by MEAN_PSF_TO_FIBER_SPECFLUX; 2-parameter (scale, tilt in log λ) correction.
   Validate on stars (0 ± 0.05 mag) before galaxies. SDSS: ugriz vs FIBER2MAG/FIBERMAG, quadratic.
2. Milky Way extinction: EBV (catalog / SFD) × 0.86, Fitzpatrick99, R_V=3.1; same for total fluxes.
3. Rest frame: log10 grid Δ=1e-4 (rebin DESI with astrolib.exec_logrebinning_fluxivarmask),
   integer shift round(log10(1+z)/1e-4), F_λ×(1+z), ivar/(1+z)^2, grid 1400–9900 Å.
4. Absolute magnitude: flat ΛCDM h=0.7, Ωm=0.3, ΩΛ=0.7; m_r = dereddened whole-galaxy
   (FLUX_R / CMODELMAG); K per object from its own calibrated rest-frame spectrum.
   Rest r is on the spectrum only for z<0.40 → band switches by slice:
   M_r (0–0.40), M_g (0.40–0.80), M_u (0.80–1.45), M_NUV (1.45–1.60), all synthesized
   from the spectrum. 1-mag bins apply to this band-switched magnitude.
   Decisions settled 2026-10-24: (1) K from spectrum as above; (2) spectrophotometry vs
   FIBERFLUX with spectrum × MEAN_PSF_TO_FIBER_SPECFLUX, validated on stars first;
   (3) extinction Fitzpatrick99, R_V=3.1, E(B−V)×0.86.
5. Per-spectrum measurements (rest frame): median S/N; EW+flux of Hα Hβ Hγ Hδ [OII] [OIII]
   [NII] [SII] [NeIII] (Gaussian + local linear continuum, ±30 Å sidebands); Dn4000;
   fiber fraction FIBERFLUX_R/FLUX_R. Then cell assignment (grid below).
6. Stack per cell: per-slice normalization window; weighted mean Σwf/Σw, w=ivar·s²;
   median, scatter, N per pixel, member list; one χ² outlier pass; floor 300.
   Output: 2D FITS per z-slice (sdsslib format, rows=cells) + cell table.

## Grid (ms.tex Table 3)
z: Δz=0.05, 0–1.6 (32) × M_r: 1 mag, −24..−16 (8) × color: 0.25 mag, 0–2, g−r (z<0.5) / r−z (8)
× line EW: Hα (z<0.46) / [OII], log2 bins <1 … >64 Å (8) × Dn4000: <1.15, 1.15–1.40, 1.40–1.65,
1.65–1.90, >1.90 (5). Nominal 81,920; measured on 3 axes (≥300): 1,199; target ≈7,000.

## Passes
- Pass A (done): catalog-only cell count → `01_count_grid_cells.py`, `data/grid_cells_z_Mr_color.csv`.
- Pass B (~20 h): steps 1–3, 5 over all coadds → per-spectrum table (Parquet).
- Cell assignment (minutes): catalog + pass-B table.
- Pass C (~20 h, or minutes with the rest-frame spectra stored on a 1–2 TB SSD): step 6.
- SDSS pass (1–2 h): same steps over spPlates on the internal disk.

## Running pass B (production) and the monitor
```
# one-time: compact redshift lookup (1.7 GB npz, 7 s)
python3 00_build_zlookup.py
# two workers, one BLAS thread each, low priority, detached, resumable
cd ~/github/projects_desi/programs_templates
for i in 0 1; do OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passB_measure.py --all --shard $i/2 \
   > ~/data/sdssdesi_status/passB_shard$i.log 2>&1 & done
# monitor + guard + gauge page: now programs_monitor/monitor.py (all tracks, ETA per step); 04_monitor.py is superseded
nohup python3 ../programs_monitor/monitor.py --port 8765 > ~/data/sdssdesi_status/monitor.log 2>&1 &
```
Measured (M4 Pro, 14 cores, 64 GB): 2.7 s/file per worker, peak RSS 4.1 GB per worker (chunked, compact lookup);
two workers read ~100 MB/s from exdisk2 and use ~2.5 cores. Guard (04_monitor.py): SIGSTOP all workers if available
memory < 8 GB or swap > 6 GB (resume at 16 GB); SIGSTOP the highest shard if CPU > 92% for 90 s (resume < 70%).
Overflow coadds: the exFAT disk mounts as /Volumes/exdisk1 or /Volumes/My Book; both are tried.

## Incident 2026-10-24 17:20 — exdisk2 disconnected during pass B
At 17:19:51 the last file was written (10,281 of 34,508 done, 30%); then a FITSIO read error on
coadd-main-bright-30037.fits was followed by `/Volumes/exdisk2` disappearing from `mount` and `diskutil list`
(USB disconnect or drive failure; exdisk1 stayed mounted). Workers walked the rest of the manifest as MISSING
(harmless: no output is written for MISSING files) and were killed. Nothing lost; outputs are atomic.
**Resume once exdisk2 is back** (check `mount | grep exdisk2`, then `ls /Volumes/exdisk2/desiredux/DR2`):
```
cd ~/github/projects_desi/programs_templates
for i in 0 1; do OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passB_measure.py --all --shard $i/2 \
   >> ~/data/sdssdesi_status/passB_shard$i.log 2>&1 & done
for p in $(pgrep -f '02_passB_measure.py --all'); do nohup caffeinate -i -w $p >/dev/null 2>&1 & done
```
Existing outputs are skipped; the file that raised the read error (bright-30037) is re-read and should be
checked (FITS CHECKSUM) if it fails again. Consider `04_monitor.py` adding "disk unmounted → pause workers".

## Operational rules
- Scripts here; caches in `~/data/sdssdesi_cache/`; never write to the external disks.
- Every pass resumable (per-file checkpoint), detached (setsid/nohup), Mac kept awake.
- Env: SPECTRO_REDUX, DESI_REDUX, PYLIB, DESILIB, RR_TEMPLATE_DIR (see ~/.zshrc).
- DESI coadd paths: /Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main/{dark,bright}/...,
  969 overflow files on /Volumes/My Book/desiredux/... (list: desiredux/DR2_overflow_on_exdisk1.txt).

## Checks
Step 1 stars 0±0.05 mag; step 3 line centroids within 1 pixel; step 4 M_r vs luminosity function;
step 6 cell count vs 7,000 (tune floor), stack S/N, visual inspection.

## Pass Q — quasar measurements (`05_passQ_measure.py`), status 2026-10-24
Imports pass B and runs its full chain on the QSO rows only (SPECTYPE=QSO, ZCAT_PRIMARY), then on a 600–9900 Å rest
grid: power-law slope ALPHA_LAM (line-free windows, χ² in linear flux), broad lines Lyα+NV, SiIV, CIV, AlIII+CIII],
MgII, Hβ, Hα and narrow [OIII] (flux, EW, σ, centroid velocity offset), line redshifts Z_MGII/Z_OIII/Z_CIV/Z_HBETA,
forest coverage/transmission, AI-type BAL proxy at CIV/SiIV. Same CLI as pass B; output `~/data/sdssdesi_passQ/`.
Tested on one coadd (112 quasars, 1.5 s); α_λ = −1.60 at S/N>5, CIV blueshift −313 km/s. Full run ≈ 13 h with two
workers (read-limited) — run only after pass B has finished. Design notes and catalog locations:
`docs/HANDOVER_quasars.md` §8.
