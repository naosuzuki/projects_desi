# Handover: stellar templates and the Gaia-based flux calibration (new thread)

> Written 2026-10-24 20:56 by the galaxy/quasar session as the brief for this thread. The **living status document of the
> stars thread is `docs/HANDOVER_stars.md`** (steps, decisions, census, grid, test results); the runbook is `programs_stars/README.md`.
> This file is kept as the original brief; items 1-3 of §4 are done, 4-6 are open (see docs/HANDOVER_stars.md §4).

Date: 2026-10-24. From: the "projects_desi" session (galaxy/quasar/star template design; pass B and the SDSS quasar pass
running). Purpose: let a fresh session start the **star work** — ms.tex §3.5 "Stars: combining Gaia and DESI" — without
re-deriving the context. Read this, then `~/github/papers/26_sdssdesi_spectra/ms.tex` §3.5 (sec:gaia), §3.9 (common
recipe) and `docs/HANDOVER.md` §1–2 (machine, data). Nao commits ms.tex. This file lives in
`~/github/projects_desi/programs_stars/`, where the star scripts go (numbered, same conventions as
`programs_templates/` and `programs_quasars/`; a README.md there once scripts exist).

## 1. The star plan in one paragraph (ms.tex §3.5)
Stars get two independent flux-calibration handles from Gaia DR3: (1) the **apparent** flux at the epoch — XP spectra
(330–1050 nm, absolutely calibrated) for G<17.65, synthetic PS1/LS/SDSS-u photometry for fainter stars; (2) the
**absolute** flux through the parallax, which places the star on the HR diagram, where stars at the same (M_G, (BP−RP)_0,
[Fe/H]) must have (nearly) the same spectrum. Strategy: quality cuts → calibration 1 (smooth response T_i(λ) per
spectrum against XP / photometry) → 3D extinction at the parallax distance (Bayestar19 / Edenhofer+24, not SFD) →
HR coordinates → calibration 2 (HR-cell median templates, residual response solved jointly per exposure / fiber /
petal, iterated; anchors: XP and DA white dwarfs) → validation (repeats, held-out XP, synthetic − observed photometry).
Products: (a) flux-calibrated DESI stellar spectra; (b) an empirical stellar library indexed by HR position and
[Fe/H] = the **star templates**; (c) per-exposure / fiber / petal response corrections that feed Step 1 of the galaxies
and quasars observed in the same exposures (this is why pass B stores the per-object tilt CAL_B_EST but does not apply
it: the tilt is to be solved from stars per exposure/petal — see programs_templates/README.md).

## 2. Data on this Mac (Naos-Mac-mini; see HANDOVER.md §2 for the general layout)
| data | where | notes |
|---|---|---|
| Gaia DR3 × DESI DR2 crossmatch (12.6 M stars) | `~/github/projects_gaia/csvfiles/gaiadr3_desidr2pix_star{,1..5}.csv` (7.4 GB) | columns: source_id… phot_g_mean_mag, bp_rp, parallax, parallax_over_error, non_single_star, phot_variable_flag, ruwe (+ TARGETID, astrometry, GSP-Phot in the full files — check headers) |
| Gaia × DESI DR1, × SDSS DR8/DR17 stars | `gaiadr3_desidr1_star*.csv`, `gaiadr3_sdssdr8_star.csv` (953 MB), `gaiadr3_sdssdr17_star.csv` (537 MB), `*_combined.csv` | SDSS crossmatches include PLATE/MJD/FIBERID |
| Gaia XP spectra, SDSS stars | `~/github/projects_gaia/data/gaiadr3_xpspec_sdssdr8_star.fits` (209,471 × 251 px, flux, error, ids), `..._sdssdr17_star.fits` | sampled 10^(2.5186+0.002 i) nm, i=0..250 (330–1049 nm), flux ×1e17; made with `programs/01_read_xpspecs.py` (gaiaxpy.calibrate) |
| Gaia XP spectra, DESI stars | **not downloaded** (6.6 M stars with G<17.65) | fetch via Gaia archive DataLink by source_id in chunks (gaiaxpy `calibrate` on the continuous XP coefficients); needs `pip install gaiaxpy` (not installed on this python3.13) |
| DESI DR2 loa coadds | `/Volumes/exdisk2/desiredux/...` (+969 overflow on exdisk1) | read-only; **external disks are the bottleneck — do not start a heavy pass while pass B runs** |
| zall-pix-loa cache | `~/data/sdssdesi_cache/zall_loa.parquet`, `zlookup_compact.npz` | SPECTYPE=STAR rows: 6.8 M after cuts (main dark+bright) |
| pass B output (per spectrum, incl. **all stars**) | `~/data/sdssdesi_passB/{dark,bright}/hpx<N>.parquet` (complete Sunday 2026-10-25 ~07:00) | SNR_*, synthetic fiber g/r/z vs FIBERFLUX (CAL_DG, CAL_DR, CAL_A, CAL_B_EST tilt, CAL_DZ_CHECK), EBV, PHOTSYS, TARGET_RA/DEC, Z (= radial velocity for stars), ZWARN |
| SDSS spPlates DR17 / DR8 | `~/sdssredux/v5_13_2/`, `~/sdssredux/DR8/{26,103,104}/` | internal disk, free to use any time |
| spAll / specObj | `~/sdssredux/spAll-v5_13_2.fits` (CLASS, PSFMAG, FIBER2MAG, EXTINCTION…), `specObj-dr8.fits` (no magnitudes; plugmap MAG = fiber mags) | |
| dust maps | `~/.dustmapsrc` → `/Volumes/exdisk1/data/photoz/dustmaps` (only **sfd/** present) | Bayestar19 (`dustmaps.bayestar`, ~1 GB) and Edenhofer+24 not fetched — needed for item 3 |
| DESI MWS stellar parameters (RVS/SP, [Fe/H]) | not on disk | DR1 VAC: https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/ (needs Nao's approval to download; DR2 not public) |
| redrock star templates | `~/github/redrock-templates/rrtemplate-STAR-{A,B,CV,F,G,K,M,WD}-v0.1.fits` | reference / validation |
| HR-diagram figures | `~/github/papers/26_sdssdesi_spectra/figures/hr/`, made by `projects_gaia/programs_plot/01_plot_HRdiagram_SNR.py` | numbers in ms.tex §3.5: 3.9 M DESI stars at parallax S/N>10, 934 k at >50, 240 k Gaia non-single excluded |

## 3. What pass B already gives for stars (useful immediately when it finishes)
Pass B ran on *all* spectra, so every DESI star has: observed S/N, synthetic DECam/BASS g, r, z of the spectrum vs the
LS fiber fluxes (CAL_DG, CAL_DR; the independent z-band check CAL_DZ_CHECK), the scale CAL_A and the tilt estimate
CAL_B_EST. Validated on 20 test coadds: catalog − synthetic = +0.03 (r), +0.05 (g), σ = 0.05–0.06 mag; after the scale
the z band agrees to +0.006 ± 0.053; the tilt scatter (0.04 mag in g−r) does not depend on S/N → it is photometric or
per-fiber, i.e. the thing calibration 2 is meant to solve. First star task once pass B is done: map CAL_B_EST (and
CAL_DG−CAL_DR) against petal, fiber, exposure (from the coadd FIBERMAP / EXP_FIBERMAP) and against BP−RP to separate
the per-fiber response from the color term of the filters.

## 4. Suggested first tasks (catalog-level, safe while pass B and the SDSS quasar pass run)
1. **Sample table**: from the Gaia × DESI DR2 csvs build one Parquet (TARGETID, source_id, G, BP, RP, parallax, σ,
   RUWE, non_single, variable, GSP-Phot Teff/logg/[M/H] if present) with the Lindegren+21 parallax zero-point applied;
   apply the §3.5 cuts (parallax S/N>10, RUWE<1.4, not non-single, not variable, ZWARN=0) → ≈3.9 M stars. Cache in
   `~/data/sdssdesi_cache/stars/`.
2. **HR grid design** (the star analogue of the galaxy/quasar grids): bins in M_G (0.25 mag), (BP−RP)_0 (0.05 mag),
   [Fe/H] (3–5 bins) → count stars per cell at S/N>5 (predict S/N from FIBERFLUX_R and TSNR2 as for quasars, or wait
   for pass B); floor 300 → number of templates. White dwarfs and giants separately. Write the counts into ms.tex.
3. **3D extinction**: install Bayestar19 into the dustmaps dir (~1 GB; needs approval); A_G, E(BP−RP) at the parallax
   distance for the sample; compare with SFD for nearby stars (expect SFD to over-correct).
4. **Gaia XP for DESI stars**: write the DataLink fetch (chunks of 5000 source_ids), store as `gaiadr3_xpspec_desidr2_star_<k>.fits`
   in the projects_gaia layout; start with the ≈100 k brightest (G<15) as the calibration-1 anchors. Needs gaiaxpy.
5. **Calibration 1 prototype on SDSS**: the XP spectra for SDSS DR8 stars are already on disk → fit a smooth response
   T(λ) (low-order polynomial in log λ, or spline) of the SDSS spectrum to the XP spectrum for the 209 k DR8 stars;
   this tests the method, the XP systematics (XP is 330–1050 nm at R≈50) and gives the SDSS-I/II spectrophotometry
   check. spPlates are on the internal disk, so this can run now (`sdsslib/pylib/sdss_db.py` reads them).
6. **DA white dwarfs**: assemble the anchor list (Gaia WD candidates × DESI; DESI DR1 MWS WD VAC if downloaded) and
   the model-atmosphere grid to use (e.g. Koester models) for the absolute zero point.
7. When pass B is done (Sunday): task of §3 above, then calibration 2 design on one petal/exposure as a test case.

## 5. Open decisions for Nao
(a) the HR grid bin sizes and the [Fe/H] source (DESI MWS SP vs Gaia GSP-Phot); (b) the response model T(λ) form and
order; (c) whether calibration 2 is solved per exposure or per (night, petal); (d) downloads: Bayestar19, DESI MWS VAC,
Gaia XP for DESI stars (~GBs); (e) nothing on the directory: `programs_stars/` is set.

## 6. Conventions and rules
Same as HANDOVER_quasars.md §3 and §7: never write to the external disks; nothing heavy on them while pass B runs
(ETA Sunday 07:00); keep ≥16 GB free; scripts in the repo, caches in `~/data/`; document in ms.tex (Nao commits);
every pass resumable and detached. Session filed under the sidebar group **projects_desi**.
