# Handover: stellar templates (stars thread)

Date: 2026-10-24 (night). From: the "projects_desi" session, after reading `docs/HANDOVER.md`, `docs/HANDOVER_quasars.md`,
`programs_templates/README.md`, `programs_quasars/README.md` and ms.tex §3.5 ("Stars: combining Gaia and DESI"), §3.8 (table
of the three template sets) and §3.9 (processing recipe). Purpose: the steps to generate the **stellar** templates on this Mac,
in the same shape as the galaxy (pass B) and quasar (pass Q) threads. Scripts live in `programs_stars/` (README there).

Companion: `programs_stars/HANDOVER_stars.md` is the brief the galaxy/quasar session wrote for this thread at 20:56 (data table,
XP-fetch and DA-white-dwarf tasks, open decisions); this file is the living status document.

## 1. Where the project stands (2026-10-24 ~21:00)
- **Pass B** (galaxy per-spectrum measurements, all 34,508 DESI coadds): exdisk2 dropped at 17:20; resumed 18:14 with two workers
  (`passB_w{0,1}of2.log`); 12,838 files done (37 %), 13.9 files/min, monitor ETA **Sunday 23:04** (http://127.0.0.1:8765/monitor.html).
  The ~6,900 `MISSING` lines per log in `passB_shard{0,1}.log` are from the pre-incident run; all coadds are present
  (dark 15,656 on exdisk2 + 969 on exdisk1 = 16,625; bright 17,883). Pass B writes **all rows of every coadd, stars included**.
- **Pass Q SDSS** (`programs_quasars/01_passQ_sdss.py`, internal disk) is running (~1,060 of 8,084 plates at 21:00).
  Pass Q DESI waits for pass B (disk-bound).
- Nothing for stars existed before tonight except the ms.tex plan and the Gaia × DESI crossmatch in `~/github/projects_gaia`.
  `programs_stars/` was empty.

## 2. What exists for stars now (built tonight, catalog level, internal disk only)
| file | rows | content | script |
|---|---|---|---|
| `~/data/sdssdesi_cache/zall_stars_main.parquet` | 12,180,554 | every main-survey SPECTYPE=STAR row of zall-pix-loa with SUBTYPE, REF_ID (Gaia source_id), Gaia G/BP/RP, PARALLAX, PM, MWS/SCND/DESI_TARGET, TSNR2, fluxes | `programs_stars/00_build_star_lookup.py` (2 min) |
| `~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet` | 3,919,204 | Gaia quality sample with M_G, M_G0, (BP−RP)_0, plx S/N, has_xp, GSP-Phot Teff/logg/[M/H]/A_G | `programs_stars/01_hr_cells.py` (30 s) |

Census (good = ZCAT_PRIMARY, ZWARN=0, COADD_FIBERSTATUS=0):
- 11,597,017 good main-survey stars: bright 5,993,940, **backup 4,799,617**, dark 803,460. Dark+bright = **6,797,400** (= ms.tex).
  The backup-program coadds are **not on disk** (only main dark+bright were downloaded): the brightest, highest-S/N stars with the
  most XP spectra live there. Decision for Nao: download main/backup healpix coadds (size to be checked on the DESI file listing) or
  build the stellar library from dark+bright only.
- Redrock SUBTYPE of the 6.80 M: K 3,392,606 · G 1,994,968 · M 854,098 · F 449,554 · WD 69,416 · A 33,458 · B 2,756 · CV 544.
- 5,982,399 have a Gaia match in the fibermap (REF_CAT=G2); 2,213,951 have G<17.65 (XP spectra exist). G percentiles 5/50/95 = 16.4/18.0/19.2.
- Radial velocity cz: 1/50/99 percentiles −301/−14/+237 km/s; 1.4 % beyond ±300 km/s; median Redrock error 1.2 km/s.
  One pixel of the log grid is 69 km/s, so **stars need the RV shift** before stacking (integer shift leaves ≤35 km/s; sub-pixel shift
  is cheap on the log grid and recommended for stars — decide).
- Gaia × DESI crossmatch (projects_gaia): 12,623,328 rows (12,105,597 unique TARGETID); quality cuts of ms.tex (plx S/N>10, RUWE<1.4,
  not non-single, not variable) → **3,919,204** (= ms.tex), plx S/N>50: 933,526; with XP: 3,121,305; with GSP-Phot [M/H]: 3,548,478.

HR cells (M_G0 × (BP−RP)_0 [× 4 bins of GSP-Phot [M/H]: <−1, −1…−0.5, −0.5…0, >0, + "unknown"]), stars inside at the floor:
| grid | no [M/H], floor 100 / 300 / 1000 | with [M/H], floor 100 / 300 / 1000 |
|---|---|---|
| 1.0 mag × 0.2 | 175 / 131 / 97 cells (3.90 M stars at 300) | 417 / 315 / 214 |
| **0.5 mag × 0.1** | 456 / **319** / 214 (3.87 M) | 1,044 / **713** / 466 (3.81 M) |
| 0.25 mag × 0.05 | 1,149 / 761 / 459 | 2,468 / 1,600 / 720 |
Subsets at 0.5 × 0.1, floor 300: plx S/N>50 → 168 cells (0.91 M); XP available → 279 cells (3.09 M). Regions: M dwarfs (M_G0>10, BP−RP>2)
216 k; white-dwarf region 43 k; giants (M_G0<2, BP−RP>0.9) 109 k. Extinction here is GSP-Phot A_G where available (not the 3D map yet).

## 3. The steps (Nao's formulation, 2026-10-24 night; SDSS first because the spPlates are on the internal disk)
Goal for now: stacks of the HR-grid cells whose stars are alike; EM-PCA comes after.

| step | what | SDSS (DR17 + DR8) | DESI |
|---|---|---|---|
| 1 | **Catalog**: stars with a Gaia DR3 parallax (crossmatch of projects_gaia, nearest match per fiber), pipeline RV, SUBCLASS/SUBTYPE, PSF mags, ELODIE or GSP-Phot parameters, quality flag (plx S/N>10, RUWE<1.4, single, not variable), M_G, M_G0, (BP-RP)_0 | **done**: `04_sdss_star_catalog.py` → `sdss_stars_gaia.parquet`, 1,064,587 primary stars, 780,269 with parallax (DR17 266,233; DR8 514,036), quality 134,037, XP 255,653 | **done**: `00_build_star_lookup.py` (12.18 M main STAR rows, 5.98 M with Gaia) + `01_hr_cells.py` (3.92 M quality) |
| 2 | **Spectrophotometric calibration**: synthetic g,r vs PSF mags (SDSS; the pipeline calibrates to PSF flux) / fiber flux × MEAN_PSF_TO_FIBER_SPECFLUX (DESI); scale CAL_A, tilt stored; **Gaia G/BP/RP from the raw spectrum vs the catalog** as the independent, more precise check (σ 0.03 vs 0.05) and the tilt handle; Gaia XP spectra as the final anchor (on disk for 257 k SDSS stars; DESI to download) | `05_passS_sdss.py` (written next; reuses the quasar thread's spPlate reader + `sdss_scale`) | `02_passS_measure.py` **tested**: G +0.06 ± 0.03, BP−RP +0.045 |
| 3 | **RV → rest frame**: RV from the pipeline z (ELODIE_Z / spAll Z for SDSS, Redrock Z for DESI; σ68 ≈ 50 km/s, 1 % beyond ±300), integer shift on the log grid (69 km/s pixels), residual stored (sub-pixel shift at stack time if chosen); our own RV by cross-correlation with the cell templates in the iteration of step 6 | in `05_passS_sdss.py` (spPlates are already on the 1e-4 log grid: no rebin) | in `02_passS_measure.py` (`starlib.rv_frame`) |
| 3b | **Dust map → Milky Way extinction**: 3D map at the parallax distance (Bayestar19, later Edenhofer+24) instead of the full SFD column, which over-corrects nearby stars; F99, R_V = 3.1; columns EBV_SFD, EBV_3D, EBV_USED, EBV_SRC, A_G_USED (GSP-Phot A_G when available) | code path in `starlib.distance_columns`; **Bayestar19 download to approve** (dustmaps; only SFD is installed, on exdisk1). Until then EBV_USED = 0.86 × SFD, flagged | same |
| 4 | **Absolute flux**: distance 1/ϖ (ϖ + 0.017 mas; full Lindegren+21 later), DM, ABS_SCALE = (d/10 pc)², dereddened spectrum × ABS_SCALE = the spectrum at 10 pc; check M_G_SYN0 = G_syn − DM − A_G against Gaia M_G0 | columns in pass S; the scaling itself is applied at stack time (step 6) | same |
| 5 | **HR grid**: cells in (M_G0, (BP−RP)_0 [, [M/H]]); candidate 0.5 mag × 0.1 (×4 [M/H] bins) → DESI 319 (713) cells ≥ 300 stars, finer 0.25 × 0.05 → 761 (1,600); SDSS counts after the pass. **Test whether the stars of a cell have similar spectra**: per cell, scatter of the indices (Ca K, Balmer, Mg b, Na D, Ca II triplet, TiO5/CaH) and of the absolute-flux continuum relative to the measurement errors, bimodality tests, fraction of Hα emitters, spread of ELODIE/GSP-Phot [Fe/H]; **report the cells where different stellar groups share the same (M_G, color)** — expected: metal-poor halo turnoff vs disk F/G at the same position, reddened stars still on the sequence after an incomplete extinction correction, unresolved binaries (≈0.75 mag above the sequence), subdwarfs vs dwarfs, active vs inactive M dwarfs, DA/DB/DC white dwarfs at the same (M_G, color), giants vs dwarfs where the branches cross | `06_hr_cells_check.py` on the pass-S tables (SDSS first) | same, after pass S DESI |
| 6 | **Stack** the cells whose stars are alike (RV-shifted, dereddened, at 10 pc: weighted mean + median + N, one outlier pass, member list); split the mixed cells by the extra axis that separates the groups ([M/H], Hα, WD type) before stacking; floor per cell to decide. Output: 2D FITS per slice (sdsslib format) + cell table. **This is the goal for now.** | `07_stack_stars.py` | same |
| 7 | **EM-PCA** on the stacks (+ high-S/N individuals), per temperature class, RV as fit parameter; validation vs Redrock STAR v0.1 and Gaia RVS; ms.tex | later | later |

Where DESI waits for pass B (external disks, ETA Sunday 23:04), steps 2–6 run on SDSS now (internal disk, ≈1 s per plate → 2–3 h for
the 8,100 plates, alongside the quasar thread's SDSS pass).

### Pass S DESI test (coadd-main-bright-0, 381 stars)
Gaia catalog − synthetic: G +0.06 ± 0.03, BP +0.075, RP +0.025, BP−RP +0.045, independent of S/N and SUBTYPE — the DESI pipeline
flux is ≈6 % brighter than Gaia G (the LS fiber check in pass B gives +0.03–0.05 in the same sense) and ≈0.045 mag bluer than Gaia
BP−RP. The Gaia G and LS fiber r residuals of the same star correlate at 0.80, so the per-star scale error is in the spectrum and Gaia G measures it with σ 0.03 against 0.05 for the fiber flux. Step-1 fiber check for stars as before (g +0.04, r +0.02, z −0.01 ± 0.05). RV σ68 48 km/s, integer-shift residual σ68 23 km/s.
Indices follow spectral type (A: Balmer 5–11 Å; K: Ca K 10, Mg b 5.4, Na D 4.3 Å; M: TiO5 0.58, Na D 11 Å; 7 Hα emitters).
Cost 10 s per 381 stars after a 9 s read → the measurement is cheap; the read dominates, hence "inside the pass-Q walk".

## 3a. Step 2 and step 3b clarified (2026-10-24 night; same text condensed into ms.tex §3.1, §3.2, §3.5)

**Step 2 — spectrophotometric calibration of a star.**
1. *What the pipelines deliver.* Both pipelines tie the flux scale to standard stars fitted as point sources: DESI returns the total
   flux of a point source (the fiber loss is corrected with the PSF model), SDSS/BOSS calibrate to the PSF magnitudes of the standards.
   A star's spectrum must therefore reproduce its total (PSF) flux **directly, with no aperture factor** — unlike galaxies, where the
   DESI spectrum is first multiplied by MEAN_PSF_TO_FIBER_SPECFLUX and compared with FIBERFLUX (pass B), and unlike the BOSS fiber
   magnitudes, which are 0.3–0.7 mag fainter than the spectrum (quasar thread test).
2. *Reference photometry.* **Gaia DR3 G, BP, RP** (Vega; zero points Riello+21: AB−Vega = 0.114, 0.016, 0.356), epoch-mean 2014–2017,
   per-star precision of a few millimag to G ≈ 19; available for 5.98 M of the 6.80 M DESI stars and for every star of the SDSS catalog
   (it is the selection). Legacy Surveys g,r,z (DESI) and SDSS PSF ugriz (SDSS) remain as the cross-check of the photometric systems and
   for the 0.8 M DESI stars without a Gaia match. For G < 17.65 the **XP spectrum** (330–1050 nm, absolutely calibrated) fixes the
   continuum shape; broad bands cannot.
3. *Form of the correction.* Per star, multiplicative: the **scale from G** (one number; σ 0.03 mag between stars, dominated by the
   spectrum, see 4), the **tilt from BP−RP** (linear in log λ between the BP and RP pivot wavelengths) applied only where BP and RP are
   precise (errors < 0.02 mag) and otherwise left to Calibration 2, which solves it per exposure, petal and fiber; for XP stars a smooth
   low-order response T(λ) (polynomial in log λ, degree ≤ 3) fitted to the ratio of the spectrum degraded to the XP resolution and the
   XP spectrum. The two-band and XP corrections are stored, not applied blindly: the stacks use the scale per star and the
   exposure/petal/fiber response.
4. *What the test says* (coadd-main-bright-0, 381 stars): catalog − synthetic G = +0.06 (σ68 0.03), BP = +0.075, RP = +0.025,
   BP−RP = +0.045, independent of S/N and spectral type; the G residual of a star correlates at 0.80 with its Legacy Surveys fiber-r
   residual (pass B), so the per-star scale error is in the spectrum, not in either photometric catalog, and G measures it twice as
   precisely. The DESI pipeline flux is ≈6 % brighter than Gaia G and ≈0.045 mag bluer than Gaia BP−RP; both need confirmation at
   scale (and the BP zero point a check) before the tilt is used.
5. *Order.* The calibration compares the observed (reddened) spectrum with the observed (reddened) photometry, so it comes **before**
   the extinction correction (as for galaxies); the Gaia magnitudes are not dereddened for this step.
   Stored: SYN_GAIA_G/BP/RP, GAIA_DG/DBP/DRP/DBPRP, GAIA_EXTRAP_* (band weight in the flat padding beyond the spectrum: 1–2 % DESI,
   more for SDSS-I), CAL_A/CAL_B_EST/CAL_OK from the LS or PSF comparison, SYN_FIBER_*/SYN_G_RAW, SYN_R_RAW.

**Step 3b — Milky Way extinction of a star.**
1. *Galaxies and quasars* lie behind the whole Galactic column: E(B−V) = 0.86 × SFD (Schlafly & Finkbeiner 2011 recalibration),
   Fitzpatrick99 with R_V = 3.1, applied to the spectrum and to the photometry — settled (pass B).
2. *A star* is reddened only by the dust in front of it. At high latitude most dust lies within a few hundred parsecs, so the full SFD
   column over-corrects nearby stars (the WD sequence, the M dwarfs with the best parallaxes) and would be mistaken for a calibration
   or temperature error: an unremoved E(B−V) of 0.03 moves a star by 0.08 in G and 0.04 in BP−RP along the reddening vector.
   Use a **3D map evaluated at the parallax distance**: Bayestar19 (Green+19; Pan-STARRS + 2MASS; Dec > −30°, 3.4′–13.7′, to a few
   kpc) and Edenhofer+24 (Gaia XP based, 69 pc–1.25 kpc, 14′) for the nearby stars; SFD is the far-distance limit and the fallback
   (flag EBV_SRC = 1) where no 3D value exists. Gaia GSP-Phot A_G / E(BP−RP) are kept as a per-star **cross-check** only (model
   dependent, biased at low extinction; Andrae+23), not as the correction.
3. *Application.* The same E(B−V) dereddens the spectrum (F99, R_V = 3.1) and the HR coordinates (A_G, E(BP−RP) from the Gaia
   extinction coefficients at the star's color; the code uses A_G = 2.74 E(B−V) for now). Nearby stars (d < 100 pc) receive ≈0.
4. *Status.* Only SFD is installed (`dustmaps` on exdisk1). Bayestar19 (and Edenhofer+24) must be downloaded — approval needed; until
   then pass S writes EBV_SFD, EBV_3D = NaN, EBV_USED = 0.86 × SFD with EBV_SRC = 1, so the tables can be re-dereddened without re-reading
   the spectra (the per-star correction is a multiplication by 10^{0.4 ΔA(λ)}).
5. *Sequence for a star* (steps 2 → 3 → 3b → 4): calibrate on the reddened spectrum → RV shift → deredden with E(B−V)_3D → scale to
   10 pc with (d/10 pc)², d = 1/(ϖ + 0.017 mas) (full Lindegren+21 zero point later). Only the first two touch the spectra; the last
   two are per-star factors applied at stack time.

## 3b. Grid plan (step 5): size and cell numbers
Quality samples: DESI 3,919,204 (plx S/N > 10, RUWE < 1.4, single, not variable); SDSS 129,378 (same + ZWARN = 0; 230,246 at plx
S/N > 5). Grid in (M_G0, (BP−RP)_0) with GSP-Phot A_G where available (3D map later), four [M/H] bins (< −1, −1…−0.5, −0.5…0, > 0)
plus an "unknown" bin; counts = cells holding ≥ floor stars:
| grid | floor | DESI | SDSS | DESI + SDSS |
|---|---|---|---|---|
| 0.5 mag × 0.1, no [M/H] | 100 / 300 | 456 / 319 | 181 / 88 | 497 / 334 |
| **0.5 mag × 0.1 × [M/H]** | 100 / 300 | 1,044 / 713 | 243 / 93 | **1,096 / 734** |
| 0.25 mag × 0.05, no [M/H] | 100 / 300 | 1,149 / 761 | 267 / 92 | 1,190 / 780 |
| 0.25 mag × 0.05 × [M/H] | 100 / 300 | 2,468 / 1,600 | 274 / 75 | 2,533 / 1,627 |
Plan: **base grid 0.5 × 0.1 × [M/H]4**; the homogeneity test (step 5) runs on every cell with ≥ 100 stars (≈1,100); a cell with
≥ 1,200 stars is split into its four 0.25 × 0.05 quarters when each quarter keeps ≥ 300 (adaptive: fine where the main sequence is
dense, coarse on the giant branch and the WD sequence); templates (step 6) from cells with ≥ 300 (≈730 + splits ≈ 1,000–1,200 stacks).
Stars without a quality parallax (2.9 M DESI, 0.65 M SDSS) go to SUBTYPE/SUBCLASS × S/N fallback stacks. What SDSS adds is not
numbers but coverage: 18 % of its quality stars are metal-poor (ELODIE [Fe/H] < −1) against 5 % in DESI (GSP-Phot), so the [M/H] < −1
turnoff cells (M_G0 3.5–6.6, (BP−RP)_0 0.5–0.8) are filled by SDSS; it also adds 21 late-M cells (M_G0 12.5–13, (BP−RP)_0 2.8–3.3)
with ≥ 100 stars where DESI has < 300. Where the groups mix inside a cell (step 5 report), the cell is split along the separating axis
before stacking.

## 3c. Decisions of 2026-10-24 night (Nao) — supersede the corresponding lines above
- **Gaia XP is not required.** Its few-per-cent, color- and wavelength-dependent systematics (worst below 400 nm and above 950 nm,
  and for faint stars) are of the size we want to correct. XP stays an optional cross-check with the SDSS XP spectra already on disk;
  no DESI XP download.
- **Step 2 = joint fit (step2b_spcal.py):** Gaia G sets the gray scale, the SDSS / Pan-STARRS / Legacy Surveys bands constrain the
  tilt, and each band's zero point (plus a color term) is solved from the stars themselves; Gaia G, BP, RP are the fixed reference
  system. Pan-STARRS joins when a crossmatch is on disk. SDSS filter curves: Doi et al. 2010 (speclite sdss2010, airmass 1.3).
- **Dust maps installed** in `~/data/GAIA` (Bayestar19 0.73 GB, Edenhofer+24 mean map 3.25 GB). Units calibrated on 10.3 M DESI
  stars: E(B−V)_SF11 = 0.86 × E_Bayestar19 (Bayestar is on the raw SFD scale; matches 0.86 SFD at d > 3 kpc to 0.2 %) and
  = 1.00 × E_Edenhofer (matches Bayestar in the 0.3–1.25 kpc overlap to 0.6 %). Stars within 300 pc have 56 % of the SFD column in
  front of them. No star of either survey is south of Dec −30°, so Bayestar covers everything; Edenhofer serves the < 1.25 kpc stars
  (the white dwarfs and M dwarfs).
- **Monitoring:** the stars track of the main monitor (programs_monitor/monitor.py, http://127.0.0.1:8765) now lists Nao's steps
  1, 2–3 SDSS, 2b, 3b, 2–3 DESI, 5, 6, 7 with progress files `~/data/sdssdesi_status/step_star_*.json`; the separate 8767 page was merged into it (colored divider per track).
- **Scripts** renamed after the steps: `step1_*.py`, `step2_3_pass_{sdss,desi}.py`, `step2b_spcal.py`, `step3b_ebv3d.py`,
  `step5_hr_cells_check.py`, `step6_stack.py`, `check_pass.py`, `run_steps.sh` (see programs_stars/README.md).

## 5. Step 5 result on SDSS (2026-10-24 night; `~/data/sdssdesi_status/step5_sdss.txt`, cells in `~/data/sdssdesi_cache/stars/cells_sdss_0.5x0.1_gspphot.csv`)
Sample: 115,444 SDSS quality stars (plx S/N > 10, ZWARN = 0, S/N > 5, CAL_OK); grid 0.5 mag × 0.1 × [M/H] (GSP-Phot; bins
< −1.5, −1.5…−1, −1…−0.5, −0.5…0, > 0, unknown); 216 cells with ≥ 100 stars hold 80 % of the sample. The test works on residuals
after a plane fit in (M_G0, color) inside each cell, with the step-2b tilt applied to the colors; spreads are judged relative to the
median cell. First version (raw spreads, absolute thresholds) flagged every cell: the in-cell temperature gradient and the per-spectrum
tilt dominate the raw scatter, and Na D carries interstellar absorption (now report-only).
**Verdict: 38 cells homogeneous (alike stars), 178 mixed.** Causes of the mixed cells: FEH 66, INDEX 34, COLOR 28, BALMER+COLOR 12,
FEH+COLOR 9, BIMODAL 9, HALPHA(+…) 13, CONT 3.
Cells shared by different stellar groups, by HR region:
- **Metal-poor F turnoff** ([M/H] < −1, M_G0 3.8–5.3, (BP−RP)_0 0.64–0.75; 5–11 cells of 500–1,200 stars): Ca K residual spread 4–6×
  the typical cell, σ[Fe/H] 0.3 inside the cell (ELODIE agrees: 0.31–0.38). Halo stars of different metallicity share the same (M_G,
  color): the [M/H] axis must be finer there (0.5 dex bins below −1 are not enough; use 0.25 or split at the cell median — step 6 does).
- **F/G boundary** (M_G0 4.2–4.8, (BP−RP)_0 0.75–0.77; F 55 % + G 45 %): G-band and Ca K residuals 2–3× typical → two groups at the same
  colour, probably turnoff dwarfs and subgiants / different [Fe/H] (INDEX cause: flagged, stacked whole, to be examined on the stacks).
- **White-dwarf sequence** (M_G0 10.8–11.7, (BP−RP)_0 −0.25…−0.05; subclass W 80–88 % + O 12–15 %): Hδ/Hβ/He I residuals bimodal →
  DA versus non-DA (DB/DC/sdO) at the same position; split on the Balmer EW (BALMER). Colors also spread (COLOR; 3D dereddening of
  WDs at 100–300 pc and the hot-star tilt).
- **M dwarfs** (M_G0 9.7, (BP−RP)_0 2.45, 519 stars): Ca K bimodal, 3 % Hα emitters → active / inactive; split on Hα (HALPHA).
- **K dwarfs** (M_G0 6–7.3, (BP−RP)_0 1.2–1.6): mostly homogeneous; the FEH flags there come from the GSP-Phot [M/H] residual
  distribution (σ 0.16 = GSP-Phot noise) — treat with care; ELODIE [Fe/H] is quantized on its template grid, so its bimodality test is
  meaningless (guarded now: no BC for < 15 distinct values).
Homogeneous cells: the disk F/G/K main sequence at −0.7 < [M/H] < +0.3 (e.g. M_G0 5.2–5.7, (BP−RP)_0 0.87–0.94, 2,000–2,400 stars each;
K dwarfs M_G0 6.3–7.2 with 1,000–1,900 stars). No homogeneous cell among the WDs or M dwarfs at this grid.
Step 6 (running): stacks all 216 cells, splitting the FEH / HALPHA / BALMER cells on their axis; INDEX / COLOR / BIMODAL cells are stacked
whole with the flag in the CELLS table. Next refinement: finer [M/H] below −1, a Balmer-strength axis for the WDs, Hα for the M dwarfs.

## 6. Step 6 result on SDSS (2026-10-24 night) — the goal for now, reached for SDSS
`~/data/sdssdesi_stars/stacks_sdss_0.5x0.1_gspphot.fits` (FLUX/ERR/SCATTER/NSPEC [316 × 4447] + CELLS table; log grid COEFF0/COEFF1,
3550–9881 Å, F_λ at 10 pc), `members_sdss_0.5x0.1_gspphot.csv`, quick-look `stacks_sdss_quicklook.png`. 216 cells → **316 stacks**
after the splits (116 whole, 77 × 2 [M/H] halves, 14 × 2 Balmer halves for the white dwarfs, 9 × 2 Hα halves for the M dwarfs);
91,726 stars (step-2b scale for 82.6 k, tilt applied to 67 k; 3D E(B−V) for 86 %). N per stack 74 / 155 / 446 (16/50/84 %), minimum 12
(a split half; raise --nmin per half later). Stack S/N at 5500 Å 58 / 89 / 148; star-to-star scatter 14 % of the absolute flux (the
per-star scale + distance errors). Absolute level: F(5500 Å) of the M_G0 ≈ 4.7 cells = 4 × 10⁻¹¹ erg/s/cm²/Å, the Sun at 10 pc.
The splits are real spectral differences at equal continuum: metal-poor turnoff low/high-[M/H] halves differ by 11 % at Ca K (0.3 % at
5000–5500 Å); DA/non-DA white-dwarf halves by 20 % at Hβ; M-dwarf Hα-emitter/non-emitter halves by 13 % at Hα. Two passes over the
6,572 plates took 11 min. Next: DESI (after pass B), finer [M/H] bins below −1, Calibration 2 (per exposure/petal/fiber response from
the residuals against these stacks), then EM-PCA (step 7).

## 7. Step 7 plan (EM-PCA), agreed with Nao 2026-10-25
Normalize every stack and individual by its synthetic Gaia G flux (G flux = 1; ivar scaled by the square); keep M_G0 and the HR
coordinates as labels (10^(-0.4 M_G0) restores the absolute-flux stack). Linear flux (the fitter is linear), dereddened and
RV-corrected inputs, mean subtracted (weighted mean = component 0), zero weight on uncovered pixels (EM-PCA). One basis per
temperature class (hot/WD, A-F, G-K, M, CV). Weights: each stack counts once; S/N>20 individuals with ivar scaled so that their
total weight ≈ the stacks'. Outside the basis: amplitude and RV (fit parameters), HR labels and [M/H] (read back from the
best-fitting cell). Written into ms.tex §3.4 (Step 4, "Normalization").
**The normalizing factor is stored, not recomputed** (Nao, 2026-10-25): per stack and per individual, NORM_GFLUX = the synthetic
Gaia-G flux of the input at 10 pc (photon-weighted mean F_λ through the G passband, erg/s/cm²/Å) and the equivalent magnitudes
M_G_SYN (Vega, Riello+21 zero point) and M_G_SYN_AB. In principle the factor follows from M_G alone — F_G = C_G × 10^(−0.4 M_G)
with C_G the passband constant — but that needs the passband constant and the AB/Vega zero point, and the cell's catalog M_G0
differs from the stack's synthetic value by the calibration residuals (0.02–0.05 mag), so the exact number is kept in the
template file: a NORM table (STACK, NORM_GFLUX, M_G_SYN, M_G_SYN_AB, M_G0_CELL, N) next to the eigenspectra, and the same columns
in the individuals' table. Absolute flux = normalized spectrum × NORM_GFLUX, exactly.

## 8. Status 2026-10-25 22:30 and the shutdown of 2026-10-26 ~04:30
- Pass B complete (34,508). The shared DESI walk (`02_passQ_desi.py --all --with-stars --expfibermap`, 4 shards, Quasar thread)
  started automatically and writes the DESI star tables (`~/data/sdssdesi_passS/`) and exposure maps; ≈42 coadds/min → ≈14 h;
  it is stopped cleanly at 04:10 by `~/data/sdssdesi_status/pre_shutdown.sh` and resumed after the reboot by `post_reboot.sh`
  (docs/COORDINATION.md §5). Everything is per-file atomic and resumable.
- SDSS side: the pass is being rerun with the vacuum band edges; a chain then reruns 2b → 5 → 6 → 7 (with 5 individuals per
  stack, S/N > 20); logs step2b_sdss.txt, step5_sdss.txt, step6_sdss.log, step7_sdss.txt in `~/data/sdssdesi_status/`.
- Step 7 (`step7_empca.py`) works: own weighted EM-PCA; first SDSS run on the stacks alone gives, for the M dwarfs, vector 1 =
  temperature slope, 3 = TiO/CaH band structure, 4 = Balmer emission (activity); R² 0.93/0.97/0.98 with 1/2/3 vectors. G/K: 5
  vectors for R² 0.98. The eigenvectors are noisy where only BOSS-era stacks reach (< 3800 Å, > 9200 Å): pixels now need ≥ 50 %
  coverage and a median weight ≥ 0.3 of the typical to enter the basis.
- Next on DESI (after the walk): 2b (zero points vs Gaia at scale → zp_desi.csv for the Galaxy thread), 5, 6, 7; then Calibration 2
  from the star residuals against the stacks with the exposure maps (§3d of COORDINATION.md).

### Step 7 with individuals (SDSS, 2026-10-25 23:30, vacuum-band rerun of steps 2–6 first)
Steps 2b/5/6 rerun on the vacuum-band tables: 216 cells, 39 homogeneous / 177 mixed (unchanged within noise), 315 stacks.
EM-PCA on the stacks + ≤ 5 stars per stack with S/N > 20 (1,575 stars; individuals weighted to match the stacks as a group):
| class | stacks | stars | pixels in basis | R² with 1 / 2 / 3 / 5 / 10 vectors |
|---|---|---|---|---|
| WD | 39 | 195 | 3,290 | 0.62 / 0.74 / 0.79 / 0.88 / 0.97 |
| AF | 29 | 145 | 4,131 | 0.65 / 0.81 / 0.90 / 0.93 / 0.99 |
| GK | 179 | 895 | 3,864 | 0.71 / 0.88 / 0.92 / 0.96 / 0.99 |
| M | 68 | 340 | 2,923 | 0.80 / 0.89 / 0.91 / 0.93 / 0.96 |
Two bugs fixed before this run: the G-flux padding used the grid edge instead of each spectrum's outermost good pixels (dropped
most SDSS-I individuals), and the pixel-coverage criterion counted the low-weight individuals (halved the M-dwarf basis); it is now
evaluated on the stacks only. `empca()`/`explained()` are a stable interface: the Galaxy thread imports them.

## 9. 2026-10-27: the walk also carries the galaxy pass C
Since 00:05 the walk runs `02_passQ_desi.py --all --with-stars --expfibermap --with-galaxies --shard i/4` (Quasar thread, at Nao's
request via the Galaxy thread; COORDINATION.md §6). The first ~17,500 coadds are re-read for galaxies only, so the DESI star tables
pause at 17,462 for ≈7 h, then grow again with the rest of the walk (≈8 h more). The walk stops gracefully on SIGTERM (finishes its
coadd, writes a galaxy checkpoint); the 5-min grace in pre_shutdown.sh is enough. post_reboot.sh launches the same command.
The DESI star steps (2b, 5, 6, 7, Calibration 2) wait for the walk to finish, ≈15 h from 00:05.

## 4. Decisions for Nao (ordered by when they block)
1. Backup-program coadds: download (brightest stars, most XP) or dark+bright only? (blocks the final sample, not the pipeline)
2. Downloads: Bayestar19 (dustmaps), DESI DR1 MWS VAC, Gaia XP coefficients for the DESI quality stars (sizes to be reported before fetching).
3. Pass S runs inside the pass-Q DESI walk (recommended) or as its own pass (another ≈13 h of disk reads)?
4. RV shift: integer pixel (as galaxies) or sub-pixel for stars (recommended)?
5. Grid: 0.5 × 0.1 × [M/H]4 at floor 300 (≈710 cells) or 0.25 × 0.05 (≈1,600)? One EM-PCA basis per temperature class (recommended) or one basis?

## 5. Conventions (unchanged)
Scripts in `programs_stars/`, caches in `~/data/sdssdesi_cache/`, outputs in `~/data/sdssdesi_passS/`; never write to the external
disks; nothing heavy on them while pass B runs; every pass resumable and detached; documentation goes into ms.tex (Nao commits).
