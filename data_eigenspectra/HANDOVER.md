# Galaxy eigenspectra: how they are built (handover)

Galaxy Template thread, 2026-10-28. Rest-frame galaxy eigenspectra for 0 < z < 1.6 from DESI DR2 (loa) and SDSS
DR8 + DR17, built by a weighted EM-PCA of 5-axis cell stacks. Paper: `~/github/papers/26_sdssdesi_spectra/ms.tex`
(Sec. 3.6 galaxies, Sec. 3.9 recipe). Code: `programs_templates/` (scripts 02, 07–12, `gridcells.py`, `normwin.py`,
`galaxy_hook.py`); runbook in `programs_templates/README.md`.

## Files in this directory

| file | input | tilt corrected | status |
|---|---|---|---|
| `galaxy_empca_desi_v1.fits` | 7,641 DESI stacks | no | superseded: normalization-window jump at z = 0.60, noisy blue end |
| `galaxy_empca_desi_v2.fits` | 7,641 DESI stacks | no | provisional (DESI-only comparison) |
| `galaxy_empca_desisdss_v2.fits` | 7,641 DESI + 1,231 SDSS stacks | no | provisional |
| `galaxy_empca_desisdss_v3.fits` | 7,641 DESI (tilt-corrected) + 1,231 SDSS stacks | DESI: yes | **current** |
| `galaxy_eigenspectra_v1.fits` | = `galaxy_empca_desisdss_v3.fits` (byte-identical copy) | DESI: yes | **release v1** (Nao, 2026-10-29) |

Each `.png` shows e0–e5 of the matching file.

Release naming: `galaxy_eigenspectra_v<N>` counts released bases; `galaxy_empca_<inputs>_v<N>` are the build iterations.
Release v1 de-redshifts by a whole-pixel shift (round(log10(1+z)/1e-4)). The next release (v2) redoes everything from the
original spectra with a sub-pixel shift on the same 1e-4 dex grid (`programs_templates/restframe.py`).

FITS extensions (all files):

| HDU | content |
|---|---|
| `EIGVEC` | [11, 8496] float32. Row 0 = mean shape (no mean subtraction), rows 1–10 = eigenvectors. Rows are orthonormal over the basis pixels and 0 outside them. The sign of each row is arbitrary. |
| `COEFF` | [n_stack, 11] float32, the coefficients of each input stack (in the renormalized units, see step 8) |
| `PIXMASK` | [8496] int16, 1 = pixel in the basis |
| `VARIANCE` | weighted R² of the data reconstructed with the first k vectors, k = 1…11 |
| `WAVE` | [8496] rest-frame **vacuum** wavelength in Å: 10^(i × 1e-4), log-uniform, Δlog λ = 1e-4 (69 km/s), 1400–9900 Å |
| `CELLS` | one row per stack: cell indices ZB, MB, CB, EB, DB; N; Z_MED, MABS_MED, …; STACK_SNR; NPIX; SOURCE (`desi`/`sdss`, v2 DESI+SDSS and later); RENORM (v2 and later) |

Basis wavelength range: v1 1400–9799 Å; v2 DESI 1470–9703 Å; v2 and v3 DESI+SDSS 1470–9790 Å (8,237 pixels).

**Using them:** model a rest-frame galaxy spectrum as f(λ) = Σ_k a_k EIGVEC[k](λ) with a_k free (including a_0,
which carries the amplitude). Fit only pixels with PIXMASK = 1. The vectors are dimensionless shapes; overall flux
scale and sign are absorbed in the a_k.

## Construction, step by step

### 0. Samples
- **DESI DR2 loa**, healpix coadds, main survey bright + dark (34,508 coadds).
  - Selection: SPECTYPE = GALAXY, ZWARN = 0, ZCAT_PRIMARY, DELTACHI2 > 25, 0.001 < z < 1.6, COADD_FIBERSTATUS = 0.
- **SDSS DR8** (SDSS-I/II spPlates, 3,953 plates) and **DR17** (BOSS/eBOSS spPlates v5_13_2, 4,131 plates).
  - Selection: CLASS = GALAXY, SPECPRIMARY, ZWARNING = 0, good plate quality, 0.001 < z < 1.6.
- Read counts:

  | survey | galaxies read |
  |---|---|
  | DESI | 26.0 M |
  | DR8 | 616,630 |
  | DR17 | 1,774,533 |

### 1. Spectrophotometry: fiber-aperture flux, scale only
The aperture is the fiber because that is the light the spectrum sees.
- **DESI:**
  1. Pipeline flux × MEAN_PSF_TO_FIBER_SPECFLUX, which undoes the point-source calibration.
  2. Then one scale factor so that the synthetic g, r through the DECam/BASS-MzLS curves match FIBERFLUX_G/R. This is ivar-weighted; speclite bands are cut to 0.5–99.5 % of the cumulative response.
  3. A per-object tilt was fitted, but was too noisy to apply. It is stored as CAL_B_EST in the pass-C tables only.
  4. Spectra failing the ±0.5 mag check keep scale 1 (CALC = 0) and are still stacked.
  5. 94.7 % pass.
- **SDSS:** reference = fiber magnitude (DR8 FIBERMAG 3″ from photoPlate-dr8; DR17 FIBER2MAG 2″ from spAll) + the plate's point-source offset median(PSFMAG − FIBERMAG) of its stars (DR17 −0.82, DR8 −0.33 mag). Scale only.

### 2. Milky Way extinction
Fitzpatrick (1999), R_V = 3.1, E(B−V) = 0.86 × SFD. Photometry uses the SF11 coefficients.

### 3. De-redshift
λ_rest = λ_obs / (1+z) on the common log grid (integer-pixel shift; both survey grids have Δlog λ = 1e-4), and
F_λ,rest = F_λ,obs × (1+z). Redshifts are the catalogue spec-z.

### 4. Absolute magnitude
- Whole-galaxy luminosity (not fiber) in flat ΛCDM h = 0.7, Ωm = 0.3, ΩΛ = 0.7.
- The K-correction comes from the galaxy's own spectrum.
- The rest band switches with z so that it stays inside the observed spectrum:

  | z | band |
  |---|---|
  | < 0.35 | M_R |
  | 0.35–0.72 | M_G |
  | 0.72–1.45 | M_U |
  | 1.45–1.6 | M_NUV |

### 5. Cell assignment (`gridcells.py`)
Five axes:

| axis | bins |
|---|---|
| z | Δz = 0.05, 32 slices |
| M | 1 mag steps, −24 to −16 (8) |
| colour | dereddened photometric g−r (z < 0.5) or r−z (z ≥ 0.5), 0.25 mag steps (8) |
| emission-line EW | Hα for z < 0.46, [O II] 3726+3729 above, edges 1, 2, 4, 8, 16, 32, 64 Å (8) |
| Dn4000 | edges 1.15, 1.40, 1.65, 1.90 (5) |

- Line fluxes and EWs use **vacuum** wavelengths.
- Cells with ≥ 100 galaxies in the pass-B census are accumulated; 9,117 cells.

### 6. Normalization of each galaxy (`normwin.py`)
- Before stacking, each rest-frame spectrum is scaled so that its mean F_λ over good pixels equals 1 in a line-free window.
- The window is chosen per z slice: inside the observed coverage with a 30 Å margin, off the A band (7580–7700 Å observed), and as close as possible to 6400 Å observed.
- Galaxies with < 70 % usable pixels in the window are not stacked.

| z | window (rest frame, vacuum Å) |
|---|---|
| < 0.10 | 6000–6250 |
| 0.10–0.25 | 5400–5500 |
| 0.25–0.60 | 4600–4800 |
| 0.60–1.45 | 3000–3100 |
| 1.45–1.60 | 2100–2200 |

### 7. Stacking (`08_passC_stack.py`, `galaxy_hook.py`, `11_passG_sdss.py`, merge `09_merge_stacks.py`)
- Per cell: inverse-variance-weighted mean of the normalized spectra. The accumulators store Σwf, Σw, Σwf² and N per pixel, which gives FLUX, IVAR, SCATTER and NSPEC.
- Only cells with ≥ 200 galaxies are written.
- DESI and SDSS are stacked separately; the same cell can appear once from each survey.

| set | cells written | galaxies | median stack S/N |
|---|---|---|---|
| DESI | 7,643 | 26.0 M | 47.7 |
| SDSS DR8 + DR17 | 1,231 | 1.60 M | 66 |

### 7b. Tilt correction (DESI only; used from v3; `12_restack_response.py`)
- Source: Calibration 2 of the Stellar Template thread (`programs_stars/step8_response.py`), solved from 6.2 M DESI stars with Gaia G/BP/RP as fixed anchors.
- Per (EXPID, PETAL_LOC), averaged per coadd: `~/data/sdssdesi_products/stars/response_coadd_desi.parquet`, key (TARGETID, PROGRAM).
- Correction, multiplying the observed-frame flux:
  `CORR(λ_obs) = 10^(−0.4 B_MAG_DEX · log10(λ_obs/6000))`
  - Only the shape term B is applied. The scale is already set by the fiber flux, and every spectrum is renormalized anyway.
  - It is applied to each stored rest-frame galaxy at λ_obs = λ_rest (1+z). Then the galaxy is renormalized in its window (step 6) and re-accumulated into the **same** cell. The tilt moves Dn4000, the EWs and M by much less than a bin.
- What it contains:
  - A global B ≈ −0.207 mag/dex, the same on all petals; about 0.09 mag of colour across the DESI range.
  - A real exposure-petal spread of about 0.06 mag/dex.
  - Coadd-level B, 5/50/95 %: −0.38 / −0.21 / −0.04 mag/dex.
- Coverage: 0.35 % of galaxies had no solved B and entered uncorrected.
- About 0.21 mag/dex of per-spectrum tilt is not predictable from neighbouring stars (likely per-fiber chromatic loss). Within a cell it averages down as 1/√N.
- Validation:
  - With no correction, the re-stack reproduces pass C (same 3,349 galaxies in a 40-coadd test; fluxes within 0.1 %, the float16 store precision).
  - A uniform −0.21 mag/dex gives the expected 4500→7000 Å change (1.0375 measured vs 1.0378 expected).
- SDSS spectra have **no** tilt correction yet: SDSS has no Calibration 2 equivalent so far.

### 8. EM-PCA (`10_galaxy_empca.py`, using `programs_stars/step7_empca.empca`, Bailey 2012)
- **Input:** the cell stacks with median S/N ≥ 5.
- **No mean subtraction.** With z-dependent rest-frame coverage, a per-pixel mean would be biased by which cells cover which pixel. The mean shape is fitted jointly as vector 0, so 10 + 1 vectors in total.
- **Weights:** each stack's IVAR is rescaled to median 1 over its good pixels, so each stack counts once and the many high-S/N low-z cells do not dominate the basis. Relative per-pixel weights within a stack are kept.
- **Common amplitude (v2 and later, `--renorm 1`):**
  - The step-6 windows give each z slice its own flux scale. In v1 this leaked into the coefficients: c0 jumped by ×3 at z = 0.60, where the window changes from 4600–4800 Å to 3000–3100 Å.
  - Fix: after a first decomposition, each stack is divided by its weighted projection on e0, so all stacks share the mean-shape amplitude. Then the decomposition is redone.
  - The factor is stored as RENORM in CELLS; it carries the window change (median ×59 → ×168 across z = 0.60).
  - After the fix the coefficients are continuous across all window switches.
- **Pixel criterion:** a pixel enters the basis if ≥ 200 stacks cover it (v1: 20). Below about 1550 Å only 38–280 faint high-z stacks contribute, and their noise dominated e2 in v1.
- **Settings:** 25 iterations, float64.
- **Runtimes:** v1 16 min; v2 DESI 32 min; v2 DESI+SDSS 40 min; v3 29 min (6 threads, 6.6 GB).

## Results so far

Weighted R² vs number of vectors:

| file | 1 | 2 | 3 | 4 | 5 | 7 | 11 |
|---|---|---|---|---|---|---|---|
| desi_v1 | 0.904 | 0.942 | 0.960 | 0.968 | 0.986 | 0.996 | 0.9992 |
| desi_v2 | 0.784 | 0.934 | 0.954 | 0.983 | 0.993 | 0.997 | 0.9988 |
| desisdss_v2 | 0.861 | 0.911 | 0.978 | 0.984 | 0.993 | 0.997 | 0.9988 |
| **desisdss_v3** | 0.895 | 0.951 | 0.983 | 0.987 | 0.994 | 0.998 | 0.9989 |

v1's high first term is mostly the arbitrary window scale; from v2 on the first term is shape only. With the tilt
correction (v3) the basis is more compact: more variance in the first 2–4 vectors than v2 with the same inputs.

Physical reading of desisdss_v3 (the signs are arbitrary):

| vector | content |
|---|---|
| e0 | mean galaxy (continuum, 4000 Å break, [O II], Hβ, [O III], Hα) |
| e1 | smooth UV-to-red continuum slope with [O II] (age / star-formation sequence) |
| e2 | red continuum (> 7000 Å) vs [O III] |
| e3 | [O III] vs [O II] and Hα ([O III]/[O II]: ionization) |
| e4 | [O III] + Hβ + Hα strength |
| e5 | [O II] and Hα+[N II] vs [O III] (metallicity / ionization) |

(desi_v2: e0 mean, e1 continuum slope, e2–e4 [O III]/Hβ, e5 Hα+[N II] vs [O III].)

**DESI vs SDSS in the 1,168 cells stacked by both surveys** (c_k/c0, SDSS − DESI):

| | v2 (no tilt) median / NMAD | v3 (DESI tilt-corrected) median / NMAD | v3 spread across DESI cells |
|---|---|---|---|
| c1/c0 | −0.047 / 0.085 | +0.033 / 0.061 | 0.577 |
| c2/c0 | +0.039 / 0.066 | −0.032 / 0.055 | 0.246 |
| c3/c0 | −0.037 / 0.080 | −0.008 / 0.050 | 0.151 |
| c4/c0 | −0.002 / 0.040 | +0.004 / 0.023 | 0.102 |

The tilt correction cuts the survey-to-survey scatter by 25–45 %. The offsets are now ≤ 13 % of the real spread across
cells (c1: 6 %, c3: 5 %). The vectors are reordered between v2 and v3, so compare row by row only loosely.

## Known issues / to do
- The blue end (< 1600 Å) is still somewhat noisy in e1, e3, e5 (few faint z > 1.4 stacks), less so in v3.
- Features > 9300 Å, probably sky residuals; consider a red cut at ~9300 Å.
- The lowest z slice (z < 0.05) has a c0 outlier (−1.5 vs about −0.8); check those cells (aperture effects in nearby galaxies).
- DR8 M_R is missing for 0.28 < z < 0.35 (DR8 spectra end at 9200 Å, before the band switch at 0.35). BOSS stacks only 62 % of its galaxies (cells missing from the DESI-based cell table).
- A smooth-curve refinement of the DESI response from star-stack residuals, and possibly a per-fiber term, are expected from the Stellar thread. Same files and keys, so rerun `12_restack_response.py` → `09_merge_stacks.py` → `10_galaxy_empca.py`.

## Reproduce

```bash
cd ~/github/projects_desi/programs_templates
python3 12_restack_response.py --workers 4                                   # tilt-corrected DESI accumulators (~4.8 h, reads exdisk3)
python3 09_merge_stacks.py --passc ~/data/sdssdesi_passC_resp --floor 200    # DESI stacks
python3 09_merge_stacks.py --passc ~/data/sdssdesi_passG_sdss --acc 'acc_dr*.npz' --out stacks_sdss --floor 200
OPENBLAS_NUM_THREADS=6 python3 10_galaxy_empca.py --stacks ~/data/sdssdesi_passC_resp/stacks ~/data/sdssdesi_passG_sdss/stacks_sdss \
    --out ~/data/sdssdesi_passC_resp/empca --nvec 10 --niter 25 --min-cover 200 --renorm 1 --tag desisdss_v3
```

## Data locations

| what | path |
|---|---|
| DESI stacks, tilt-corrected | `~/data/sdssdesi_passC_resp/stacks/` (996 MB) |
| DESI stacks, uncorrected | `~/data/sdssdesi_passC/stacks/` |
| DESI accumulators, uncorrected | archived on `/Volumes/exdisk3/sdssdesi_archive/passC_uncorrected/acc` (symlink at the old place) |
| SDSS stacks | `~/data/sdssdesi_passG_sdss/stacks_{sdss,dr8,dr17}/` |
| Per-galaxy normalized rest-frame spectra | `/Volumes/exdisk3/sdssdesi_restframe/` (DESI), `/Volumes/exdisk3/sdssdesi_restframe_sdss/` |
| Per-galaxy tables (calibration, M, colour, EW, Dn4000, cell) | `~/data/sdssdesi_passC/tables/`, `~/data/sdssdesi_passG_sdss/tables/` |
| Eigenspectra originals | `~/data/sdssdesi_passC/empca/`, `~/data/sdssdesi_passC_resp/empca/` |
