# Stellar eigenspectra (Stellar Template thread)

Release v1 (Nao, 2026-10-29): byte-identical copies of `~/data/sdssdesi_stars/v1_20261029/empca_<survey>_0.5x0.1_gspphot_<class>.fits`
(code at git tag `stars-v1`). Construction: `docs/HANDOVER_stars.md`; scripts `programs_stars/step6_stack.py` (stacks) and
`step7_empca.py` (EM-PCA).

## Files

`star_eigenspectra_v1_<survey>_<class>.fits` (+ `.png`: mean and eigenvectors 1-4), one basis per survey and temperature class:

| class | HR region of the cell | DESI inputs | DESI R2 at 1/2/3/10 vectors | SDSS inputs | SDSS R2 at 1/2/3/10 |
|---|---|---|---|---|---|
| WD | M_G0 > 9 and (BP-RP)_0 < 0.8 | 87 stacks + 432 stars | .464 / .965 / .988 / .999 | 39 + 195 | .514 / .731 / .807 / .978 |
| HOT | (BP-RP)_0 < 0.3 (not WD) | 5 + 25 | .866 / .961 / .974 / .996 | (none) | |
| AF | 0.3 - 0.7 | 12 + 60 | .695 / .886 / .920 / .971 | 29 + 145 | .621 / .763 / .883 / .990 |
| GK | 0.7 - 2.0 | 369 + 1,832 | .807 / .940 / .971 / .988 | 179 + 895 | .713 / .893 / .930 / .987 |
| M | > 2.0 | 425 + 2,033 | .888 / .935 / .947 / .972 | 68 + 340 | .786 / .833 / .847 / .962 |

Inputs: stacks of the HR cells (0.5 mag in M_G0 x 0.1 mag in (BP-RP)_0, split by [M/H], H-alpha or Balmer where the cell is
mixed; DESI <= 300 stars per stack, median 274; SDSS median 155) plus up to 5 individual stars per stack with S/N > 20.

## Format

| HDU | content |
|---|---|
| `MEAN` | [4447] weighted mean spectrum (component 0; subtracted before the PCA). Header: `COEFF0` = 3.5502, `COEFF1` = 1e-4 |
| `EIGVEC` | [10, 4447] eigenvectors, orthonormal over the basis pixels, 0 outside them; sign arbitrary |
| `COEFF` | [n_inputs, 10] coefficients of each input (stacks first, then individual stars) |
| `VARIANCE` | weighted R2 of the inputs reconstructed with the first k vectors, k = 1..10 |
| `PIXMASK` | [4447] 1 = pixel in the basis |
| `NORM` | per input: ID, KIND (stack/star), NORM_GFLUX, M_G_SYN (Vega), M_G_SYN_AB, M_G0_CELL, COL0_CELL, FEH_CELL, N |

Wavelength: rest-frame **vacuum**, log10 lambda = COEFF0 + i x COEFF1 (i = 0..4446), 3550-9881 A, 69 km/s per pixel.
Units: every input is normalized to synthetic Gaia G flux = 1; the factor is stored, so absolute F_lambda at 10 pc
= normalized spectrum x NORM_GFLUX. Model a star as f = MEAN + sum_k a_k EIGVEC[k] (fit PIXMASK = 1 pixels), times a scale.

## How v1 was built (short)

Per star before stacking: step-2b joint spectrophotometric fit (Gaia-G scale, Gaia/SDSS/LS tilt), RV shift to the rest frame by
the pipeline RV with a sub-pixel Lanczos-3 kernel (DESI: after the flux-conserving rebin onto the observed log grid; SDSS: from the
native spPlate pixels), 3D dereddening (Edenhofer+24 / Bayestar19, F99), x (d / 10 pc)^2 from the Gaia parallax. Stacks: weighted
mean with one 4-sigma outlier pass. EM-PCA per class, 10 vectors; pixel cut `--min-wfrac` 0 for DESI, 0.3 for SDSS (0.1 for M).

Next release (v2): DESI spectra resampled once from the native pixels straight onto the rest-frame 1e-4 dex grid at the
pipeline RV (`programs_templates/restframe.desi_rest`), from the shared v2 walk; SDSS unchanged at 1e-4.
