# Plan: atmospheric (telluric) eigenspectra from stars, galaxies and quasars

Stars thread, 2026-11-01 (Nao: "check if we have telluric features in eigenspectra; if we do, compose a plan to have just
atmospheric, telluric eigenspectra from galaxies and quasars as well"). Draft for Nao and the Galaxy / Quasar threads.

## 1. What the stellar v2 eigenspectra show

Tests on the v2 bases (`~/data/sdssdesi_stars/v2/empca_*`) and on their cached S/N > 20 individuals (4,382 DESI, 1,575 SDSS).
The telluric template is tau = -ln T of the DECam standard atmosphere (`des_atm_std`, 5 A sampling, in the LSST stack:
`obs_decam_data/.../transmission_atmosphere/2012-09-12T00:00:00.ecsv`). Scripts are in the session scratchpad
(`tell_eig.py`, `tell_rv.py`, `tell_var.py`, `tell_pca.py`, `tell_tmpl.py`).

* **Means and stacks: clean.**
  * Low- and high-RV halves of the stars (median -63 / +62 km/s) line up at lag 0.0 px in every telluric window. A telluric
    feature would sit 1.8 px apart; every feature in the median spectra is stellar (K I 7665/7699, Ca triplet, TiO).
  * The median telluric amplitude per star is |a| < 0.01 (in units of the band optical depth).
  * The Galaxy thread's observed-frame stack check (`programs_templates/25_obs_residuals.py`, `plots_v3_obsframe_residuals.png`)
    agrees: <= 0.3 % in the observed frame.
  * One exception: SDSS 8900-9200 A, with a median a of -0.02 to -0.04 and an RV lag equal to the telluric one (-1.80 vs -1.90 px).
* **Individual stars: small residuals that vary from star to star.**
  * A band, WD/HOT/AF stars: the robust scatter of a is 0.020 against 0.016 for the shifted-template null. That leaves about
    0.012, i.e. about 1 % of the band optical depth, about 1.5 % at the band core.
  * DESI water bands: consistent with noise.
  * SDSS water at 8100-9200 A: in excess (0.13-0.18 vs 0.10-0.14 null).
* **Eigenvectors: weak telluric structure in the higher orders, mostly SDSS.** |corr| with tau_ss, the small-scale part of
  tau (minus a 51-px running median), counted only when it beats the template shifted by every 15-40 px:

| survey | band | eigenvectors with telluric-like structure | note |
|---|---|---|---|
| DESI | A, B, H2O 7200/8200/9000 | WD: none | the cleanest test (no stellar lines); DESI WD basis clean |
| DESI | A band | M e2 (0.74), HOT e2-e4 (0.39-0.51), GK e2 (0.36) | **ambiguous**: the TiO 7590 A head sits on the O2 A head (7595), and HOT = WD + dM composites |
| DESI | H2O 9300-9650 | GK e2/e4, M e2, AF e5/e7, HOT e2 (0.26-0.36) | plausible water |
| SDSS | H2O 8100-8400 | AF e3/e5/e8-e10 (0.37-0.51), M e2/e4/e7/e9 (0.40-0.50) | water |
| SDSS | A, B band | WD e1/e3/e4/e6 (0.20-0.39), AF e6/e9 (0.40-0.50) | O2 |

  * The telluric eigenvectors are high order, so the R^2 they carry is small, < 1 %.
  * But they will put a fake few-percent structure into any fit that uses them.
  * The radial velocities (sigma ~ 60-80 km/s, about 1 px) are too small to separate stellar from telluric features in stars.
    Galaxies and quasars can, because their redshifts spread the rest-frame features over the observed frame.

## 2. Goal

* A small set of **observed-frame, multiplicative** eigenspectra T_k(lambda_obs) per survey (DESI; SDSS-I/II and BOSS separately,
  since calibration pipelines and spectrographs differ). They describe what the pipelines' telluric and flux calibration leave
  behind: O2 A, B, gamma; H2O 7160-7340, 8130-8350, 8950-9800.
* By-products:
  * additive sky-residual eigenspectra S_k (5577, 6300, OH forest);
  * a record of fixed instrumental features (DESI b/r seam 5800 A, r/z seam 7600 A, which overlaps the A band; SDSS b/r
    dichroic 6000 A).
* Use: each class basis is fitted as rest basis(lambda / (1+z)) x (1 + sum_k a_k T_k(lambda_obs)), with a_k free per object.
  The class eigenspectra are then rebuilt from telluric-cleaned inputs (one iteration).

## 3. Method

1. **Per-object residuals, not stacks.** Telluric residuals vary from exposure to exposure (airmass, water vapour), so a stack
   averages them out; the eigenvectors pick them up from individuals and from small stacks.
   * Inputs are high-S/N individual spectra of all three classes.
   * Residual r = f / model - 1 on the observed pixel grid. The model is the object's fit with its own thread's v2 rest-frame
     basis (stars: step 7 bases; galaxies: 10_galaxy_empca; quasars: 06_empca).
2. **Masks.**
   * Galaxy emission lines (±1000 km/s around [O II], Hβ, [O III], Hα + [N II], [S II]).
   * Quasar broad lines: a robust fit, or masks of ±5000 km/s.
   * Stars: only WD/HOT/AF enter the red, with H, Ca triplet and K I masked.
   * The redshift spread is the lever: rest-frame residuals smear out in the observed frame, telluric ones do not.
3. **Multiplicative vs additive.**
   * Telluric scales with the object flux, sky residuals do not.
   * Bright objects (S/N > 20) give T_k from r.
   * Optionally, a joint fit f - model = model x sum a_k T_k + sum b_k S_k on all objects gives S_k as well.
4. **EM-PCA.**
   * Same engine as step 7 (empca with weights ivar x model^2), per survey, on the observed grid (DESI 3600-9824 A).
   * Weight zero where fewer than half the inputs cover a pixel (the min-member-frac lesson from the stellar bases).
   * Expected components: O2, H2O, seams, flux-calibration wiggles.
5. **Validation.**
   * Same-plate / same-tile coherence: objects observed together share the atmosphere, so their a_k must correlate within a
     plate or tile far more than between plates. This is the decisive test against intrinsic features.
   * O2 coefficients against airmass; H2O coefficients against season and humidity; seam terms against petal and camera, not
     airmass.
   * Shape against the des_atm_std template (and an ESO SkyCalc / TelFit model if Nao wants to download one).
   * Hold-out: fitting held-out objects with class basis + T_k must cut chi^2 in the telluric windows and leave it unchanged
     elsewhere. The class coefficients must not change by more than their errors.
6. **Sample size.** About 20k objects per survey (mixed classes, S/N > 20), at least 5 per tile or plate for the coherence test.
   DESI at 7,781 px x (flux, ivar) float32 is about 62 kB per object, about 1.3 GB per 20k. **exdisk3 only** (the home Data
   volume is full).

## 4. Data flow (keeps "one read per coadd")

* **Phase 0, no new reads (quick look, about 1 day per thread).** Each thread takes its existing cached high-S/N individuals,
  maps its residuals back to the observed frame (lambda_obs = lambda_rest (1+z); one extra interpolation, fine for a quick
  look) and writes a common table:
  * Format: npz or parquet with SURVEY, ID, Z, PLATE/TILE, MJD, wave_obs grid, r, w.
  * Stars: done here for 5,957 stars.
  * The stars thread runs the joint observed-frame PCA and the plate/tile coherence test.
  * Go/no-go for Phase 2: do coherent T_k appear, consistent across classes, and correlated within plates?
* **Phase 1, metadata.**
  * DESI: airmass per exposure (DR2 exposures / tiles table; not on disk yet). The coadd-level value is the EXPTIME-weighted
    mean over the coadd's exposures (expfibermap tables: `~/data/sdssdesi_expfibermap`).
  * SDSS: AIRMASS from the spPlate headers (read with the plate in the next pass; spAll / specObj here have no airmass column).
* **Phase 2, observed-frame telluric sample in the next shared DESI walk.**
  * A hook (alongside star_hook / galaxy_hook / native_hook) writes observed-frame flux / ivar / mask at native pixels for the
    selected sample (all classes, ~5-20 per coadd file).
  * No extra read: it rides on the v2 sub-pixel walk. SDSS is the same in the next spPlate pass.
* **Phase 3, the telluric basis.**
  * EM-PCA per survey, validation (section 3.5), products:
    `telluric_<survey>.fits` (T_k, S_k, PIXMASK, header) and per-object coefficient tables.
* **Phase 4, use.**
  * Each thread refits with class basis x (1 + sum a_k T_k) and rebuilds its eigenspectra once from the cleaned inputs.
  * Compare the old and new bases in the telluric windows.

## 4b. Phase 0 status (2026-11-01)

* The plan and the file format went to Galaxy and Quasar (format as in section 4).
* **Stars done:** `programs_stars/telluric_phase0.py` wrote `/Volumes/exdisk3/sdssdesi_telluric/phase0/stars_desi.npz`
  (4,382 stars, 1,871 tiles, 214 with >= 5 stars) and `stars_sdss.npz` (1,575 stars, 614 plates, 63 with >= 5).
  * Fit: each star with its class's v2 basis (mean + 10 eigenvectors), with the telluric windows at zero weight.
  * Per star, the residual is fitted as a tau + linear baseline in each window.
* Coherence test: ANOVA of a across plates/tiles with >= 5 stars, against 500 label shuffles. The control windows use the A-band
  shape moved into telluric-free regions.

| window | DESI scatter / err | DESI F (p) | SDSS scatter / err | SDSS F (p) |
|---|---|---|---|---|
| B band | 0.041 / 0.031 | 1.27 (0.012) | 0.031 / 0.018 | 3.18 (<0.002) |
| H2O 7200 | 0.084 / 0.061 | 1.16 (0.06) | 0.074 / 0.048 | **10.2** (<0.002) |
| A band | 0.011 / 0.007 | 1.27 (0.014) | 0.011 / 0.006 | 5.09 (<0.002) |
| H2O 8200 | 0.063 / 0.041 | 1.27 (0.008) | 0.082 / 0.045 | 3.56 (<0.002) |
| H2O 9000 | 0.046 / 0.032 | 1.44 (<0.002) | 0.086 / 0.057 | 4.52 (<0.002) |
| H2O 9400 | 0.020 / 0.011 | 1.57 (0.006) | - | - |
| control 6100 / 7850 / 8600 | 0.008, 0.005, 0.006 | 1.11, 1.36, 1.56 | 0.007, 0.005, 0.005 | 1.91, 5.17, 1.38 |

* **SDSS: telluric residuals shared within a plate.**
  * Water at 7200 is strongest, then the A band, water at 9000, water at 8200 and the B band.
  * The median water amplitude at 8200-9200 is -0.01 to -0.02, i.e. the pipeline over-corrects slightly on average.
* **DESI: per-star residuals of about 1 % of the band optical depth (A band), about 1.5 % at 9300-9650.**
  * Tile coherence is weak and no stronger than in the control windows.
  * Healpix coadds mix tiles, so the coherence test needs per-exposure weights (expfibermap EXPID), not the dominant tile.
* **Control windows also show plate coherence** (SDSS 7850 A: F = 5.2, at about half the telluric amplitude).
  * Per-plate flux-calibration structure exists apart from telluric.
  * The joint basis must carry such terms, and telluric coherence must always be judged against controls.

## 5. Decisions for Nao

1. Go ahead with Phase 0 in all three threads? Galaxy's `25_obs_residuals.py` already does the stack-level version.
2. Who builds the joint telluric basis? The proposal is the stars thread, in a neutral `programs_telluric/`, with each thread
   exporting its residual table.
3. External telluric model (ESO SkyCalc or TelFit) for validation, or des_atm_std only?
4. DESI airmass: download the DR2 exposures table (small, about tens of MB) to exdisk3?
5. Until a telluric basis exists, should the stellar v2 bases simply mask the telluric windows (as 10_galaxy_empca
   --mask-telluric does)? For stars that removes real stellar features (K I, TiO, Ca triplet near the bands), so the
   recommendation is no; use T_k instead.
