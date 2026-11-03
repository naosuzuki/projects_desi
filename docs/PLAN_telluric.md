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
    SDSS-I/II is already on disk: `~/sdssredux/plates-dr8.fits` (2,880 plates) has AIRMASS, HA, SEEING50 and TAI per plate.
    BOSS needs the DR17 plate list (plates-dr17.fits, small), which is not on disk.
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

* **Quasar SDSS** (`quasar_sdss.npz`: 18,627 quasars, 4,755 plates; native spPlate pixels, v3 window-normalized basis):

| window | scatter / err | plate F (p) |
|---|---|---|
| B band | 0.091 / 0.071 | 1.21 (0.002) |
| H2O 7200 | 0.235 / 0.165 | 1.19 (<0.002) |
| A band | 0.028 / 0.020 | 1.11 (0.04) |
| H2O 8200 | 0.223 / 0.156 | 1.14 (0.02) |
| H2O 9000 | 0.254 / 0.190 | 1.43 (<0.002) |
| H2O 9400 | 0.077 / 0.045 | 1.80 (<0.002) |
| control 6100 / 7850 / 8600 | 0.017 / 0.020 / 0.019 | **1.03 / 1.00** / 1.35 |

  * **The telluric bands are plate-coherent and the controls are not** (6100, 7850). The stellar SDSS controls (7850 A: F = 5.2)
    probably reflect SEGUE plates that target similar stars, so quasars are the cleaner coherence probe.
  * **Mean water under-correction at 8950-9200 A:** the weighted mean R is -1.1, -1.4 and -2.1 % in the S/N tertiles (down to
    -8 % at the deepest pixels). It grows with S/N, so it's multiplicative (telluric), not additive sky. The A-band mean is about 0.
  * **Sky residuals dominate the small-scale variance in the red.** In a PCA of the high-passed residuals of the top-S/N third
    (6000-9200 A), PC1 and PC3 are the OH forest (7250-9200 A). The telluric basis needs OH lines masked or modelled as additive
    S_k, which needs a flux scale per object (requested: FLUX_SCALE, the median observed model flux).
  * **Galactic ISM in the observed frame:** Ca II K/H (3934/3969 A) and Na D (5890/5896 A) appear in the mean residual. An
    observed-frame basis will carry them; the DIB at 6284 A sits on the O2 gamma band. They should be tied to E(B-V), not airmass.

* **Galaxies** (`programs_templates/26_telluric_phase0.py`): `galaxy_sdss.npz` has 19,484 galaxies on 844 plates; `galaxy_desi.npz`
  has 12,287 galaxies on 3,531 dominant tiles, plus the per-exposure EXP_* lists. Both have median z 0.08.

| window | SDSS median a, scatter / err | SDSS plate F | DESI median a, scatter / err | DESI tile F |
|---|---|---|---|---|
| B band | +0.013, 0.067 / 0.047 | 3.1 | +0.011, 0.085 / 0.062 | 1.6 |
| H2O 7200 | -0.038, 0.150 / 0.103 | 5.4 | +0.020, 0.164 / 0.126 | 1.4 |
| A band | -0.002, 0.017 / 0.012 | 4.3 | +0.002, 0.019 / 0.014 | 2.7 |
| H2O 8200 | -0.034, 0.128 / 0.090 | 5.2 | -0.003, 0.113 / 0.083 | 1.7 |
| H2O 9000 | **-0.136**, 0.161 / 0.120 | **12.3** | -0.032, 0.091 / 0.070 | 1.7 |
| H2O 9400 | -0.002, 0.039 / 0.025 | 4.3 | +0.009, 0.036 / 0.024 | 2.2 |
| control 6100 / 7850 / 8600 | 0.011, 0.011, 0.010 scatter | 2.6 / 4.8 / 2.1 | 0.013, 0.010, 0.009 scatter | 1.3 / 1.5 / 1.6 |

  * All p < 0.002 except DESI H2O 7200 (p 0.002).
  * **Galaxy controls are plate-coherent too, as for the SEGUE stars.** Galaxies on a plate share large-scale structure (groups at
    the same z), so their rest-frame model errors land at the same observed wavelength. A one-population ANOVA can't separate
    that from telluric.
  * The SDSS water under-correction at 8920-9200 A is in the galaxies too (median a = -0.14 of the band optical depth). DESI's is
    a quarter of that.

* **Cross-population test** (scratch `tell_cross.py`). It correlates the plate means of a between populations observed on the same
  (PLATE, MJD). The atmosphere is shared; rest-frame model errors are not.

| window | SDSS galaxy x quasar r (654-200 plates) | SDSS galaxy x star r (59-33) | SDSS quasar x star r (373-202) | DESI galaxy x star r, dominant tile (765-623) |
|---|---|---|---|---|
| B band | +0.12 (p 0.006) | +0.20 (0.16) | +0.03 (0.56) | +0.02 (0.59) |
| H2O 7200 | +0.15 (0.001) | +0.12 (0.23) | +0.14 (0.001) | +0.13 (0.006) |
| A band | +0.15 (0.001) | +0.04 (0.81) | +0.14 (0.010) | +0.06 (0.17) |
| H2O 8200 | +0.17 (<0.001) | +0.14 (0.46) | +0.13 (0.015) | +0.02 (0.62) |
| H2O 9000 | **+0.39** (<0.001) | **+0.54** (<0.001) | +0.19 (0.001) | +0.04 (0.29) |
| H2O 9400 | +0.15 (0.04) | - | - | +0.06 (0.12) |
| control 6100 | +0.06 (0.19) | -0.05 (0.76) | +0.09 (0.11) | +0.02 (0.53) |
| control 7850 | +0.18 (<0.001) | +0.35 (0.03) | +0.10 (0.17) | +0.10 (0.014) |
| control 8600 | +0.10 (0.03) | +0.42 (0.005) | +0.02 (0.85) | +0.07 (0.09) |

  * **SDSS water at 8920-9200 A is a real, plate-level telluric residual seen by all three populations.** The galaxy-quasar slope
    is 0.69.
  * The other SDSS bands correlate at r = 0.12-0.17, about as much as the 7850 control. That control sits on the OH forest, so a
    shared plate sky residual (additive) is likely the common term. The joint fit has to model sky as additive (FLUX_SCALE) before
    the weaker bands can be judged.
  * **DESI shows nothing shared by dominant tile.** Healpix coadds mix exposures, so DESI needs the per-exposure weights
    (EXP_* lists; stars_desi.npz does not have them yet). DESI's per-exposure telluric fit from standard stars probably also
    leaves less.

## 4c. Final recipe and tests (2026-11-02; write-up: ms.tex, `sec:telluric`)

Nao's three steps, independently for DR8 (MJD < 55025), DR17 and DESI: telluric-only eigenspectra (`01_wd_residuals.py`, `02_basis.py`, DECam standard atmosphere as guide), correction per spectrum before stacking (`tellapply.py`), restack and EM-PCA.

* **Recipe (accepted by Nao 2026-11-02):** DR8 C = 1 + T0 + sum a_k T_k; DR17 and DESI C = 1 + T0, a_k fitted and recorded only; C = 1 red of 9200 A (DR8) / 9280 A (DR17, DESI). The 9280-9880 A water band is sky-dominated and gets no correction. DR17 H2O 1 is sky spikes and is dropped.
* **Basis tables:** `plots_telluric_<survey>.csv`. **Held-out tests:** `04_star_test.py` -> `plots_telluric_startest_{desi,sdss}.csv`; the table of galaxy, quasar and featureless WD gains is in ms.tex.
* **DR8 stars:** held-out gain 11.34 AF, 11.72 GK, 2.55 M, against -1.34, -1.24, -2.12 for the shifted eigenvectors. Rank correlation of a_k with r_TiO / r_CaT reaches 0.24, in O2 2 and O2 3, about as large in AF and GK as in M; not TiO itself, but a small stellar contribution to C on DR8 is not excluded.
* **Stars:** `programs_stars/tell_hook.py`, `step6_stack.py --telluric`, `step7_empca.py --telluric`, outputs in `/Volumes/exdisk3/sdssdesi_stars/v2_tell`.
* **Step 3, stars (done 2026-11-02):** `06_star_compare.py`, `07_plot_starcompare.py`; tables `csv_files/plots_telluric_starcompare_*.csv`, figure `plots_luminosity/plots_telluric_starcompare.png`, write-up in ms.tex. SDSS / DESI of the same cell improves in 5 of 6 bands with controls unchanged; 6266-6315 A gets slightly worse. SDSS bases lose most water vapour, DESI bases unchanged. T0 amplitude right to about 20 %, sign not settled. On 2026-11-03 `tell_hook.models` was fixed to divide by the members' smeared T0 and SDSS was restacked: water a_k stack medians fell 35-60 %, H2O 1 still follows the old term, cause open; DESI not rerun, old SDSS kept in `v2_tell/sdss_before_smearfix`. Released as `data_eigenspectra/star_eigenspectra_v3_*` the same day, on Nao's OK, with sqrt(N) stack weights as `star_eigenspectra_v4_*`, and without the DESI stars lacking BP or RP as `star_eigenspectra_v5_*`, the current release.

## 5. Decisions for Nao

1. Go ahead with Phase 0 in all three threads? Galaxy's `25_obs_residuals.py` already does the stack-level version.
2. Who builds the joint telluric basis? The proposal is the stars thread, in a neutral `programs_telluric/`, with each thread
   exporting its residual table.
3. External telluric model (ESO SkyCalc or TelFit) for validation, or des_atm_std only?
4. DESI airmass: download the DR2 exposures table (small, about tens of MB) to exdisk3?
5. Until a telluric basis exists, should the stellar v2 bases simply mask the telluric windows (as 10_galaxy_empca
   --mask-telluric does)? For stars that removes real stellar features (K I, TiO, Ca triplet near the bands), so the
   recommendation is no; use T_k instead.
