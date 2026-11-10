# Continuous normalization across redshift: recipe for all three threads

Nao, 2026-10-31. Written by the quasar thread for Galaxy and Stellar. The goal is spectra that stay continuous when different
redshift ranges are normalized in different rest-frame windows. Every normalization factor is recorded, so absolute
flux can always be restored.

## 1. The problem

A rest-frame window that is clean and covered at low z drops out of the observed range at high z, so each z range
uses its own window. If every spectrum is simply scaled to 1 in its own window, two things go wrong:

* Stacks from neighbouring z ranges do not join. One is 1 at window A, the next is 1 at window B, and the true
  ratio F(B)/F(A) is lost.
* Any axis measured in "the window" (luminosity bins, for example) jumps at the boundary. Quasars use log F_lam,abs in
  the window; at alpha_lam = -1.5 a quasar moves up about 1 bin (0.34 dex) at z = 0.8 and about 1.4 bins (0.56 dex) at z = 2.0.

PCA on such stacks spends its first eigenvectors on the normalization pattern, not on spectral variation.

## 2. Recipe

1. **Windows.** Pick line-free continuum windows and the z range each serves. Quasars (qsolib.NORM_WINDOWS):
   5100-5130 A for z < 0.8, 3000-3050 A for 0.8 <= z < 2.0, 1275-1290 A for z >= 2.0. Put the boundaries on z-slice edges.
2. **Normalize and record.** Divide each spectrum by its ivar-weighted window mean (NORM), with a fallback (power-law fit)
   when the window is less than half covered. Record per spectrum:
   * NORM and NORM_SRC (0 window mean, 1 fallback, 2 none);
   * the mean of the normalized spectrum in **every** window, plus the good-pixel fraction there (W*, COV*).
   The own window reads 1. Any other covered window gives that spectrum's flux ratio directly.
   Quasars: the stacker writes NORM/NORM_SRC into each rest-frame store file. `programs_quasars/09_norm_table.py` gathers
   them with the window means into `quasar/v2/qso_norm_table.parquet` (one row per stacked quasar).
3. **Measure the ratios between adjacent windows** on spectra that cover both (COV >= 0.8 in each, NORM_SRC 0).
   * As a function of absolute magnitude at the reference window: 0.5 mag bins where every bin of the populated range
     holds >= 100 spectra, otherwise 1.0 mag.
   * Split by survey and by z inside the overlap range. A real ratio does not depend on the survey or jump within a few
     hundredths in z; an instrumental one does (see section 4).
   * Also check the continuum slope / colour. For quasars R1 = F(3025)/F(5115) follows alpha_lam almost exactly
     (log R1 = -0.228 alpha for a power law), and at fixed slope it barely depends on luminosity.
   Quasars: `programs_quasars/10_ratio_vs_mag.py` -> `checks/qso_ratio_vs_mag.png/.csv`.
4. **Chain to one reference window.** Choose the middle window (quasars: 3000-3050 A). Multiply every stack normalized
   elsewhere by 1/R, where R is measured for that stack's own population: its magnitude bin and slope bin, taken from
   the overlap spectra. Every stack is then 1 at the reference and keeps the true level at the other windows, so
   neighbouring z ranges join.
5. **Bin luminosity at the reference wavelength too**, not at each range's own window. Otherwise the luminosity bins
   still jump at the boundaries even after the fluxes are chained.
6. **Trim stack edges.** A stack's mean jumps where the subset of members covering a wavelength changes. Drop pixels
   covered by fewer than ~50 % of the stack's members and let the PCA fill them.

## 3. Bad pixels (applies to all three threads)

SDSS DR17 6322.66-6335.78 A observed (a sky-subtraction residual peaking at 6331.40 A): `restframe.SDSS_BAD_RID`,
`restframe.sdss_badpix_ok(coeff0, npix, survey)`. DR8: none. The 6565 A candidate was stellar H-alpha absorption
from stars on the plates, not an instrument defect.

## 4. Quasar results so far (2026-10-31)

Overlap populations (instrument ranges; training set in brackets): 5100-5130 + 3000-3050 A, 614k (255k), z 0.2-1.0;
3000-3050 + 1275-1290 A, 753k (209k), z 1.8-2.4.

Stack-level (v1 stacks, `08_window_ratio.py`, `checks/qso_window_ratio*.{csv,json,png}`):

* **R1 = F(3025)/F(5115).** Overall log R1 = +0.23 (R1 = 1.69), scatter 0.23 dex. It depends strongly on luminosity at the
  faint end (log R1 = -0.31 at LBIN -2, rising to +0.31..+0.37 for LBIN >= 1) because faint, low-z quasars are red
  (host-galaxy light at 5115 A). At fixed slope bin the luminosity trend is only about +-0.05 dex. SDSS and DESI agree
  (z 0.5-0.6: +0.334 / +0.324; z 0.7-0.8: +0.382 / +0.372).
* **R2 = F(3025)/F(1282).** The two sides of z = 2.0 disagree:
  * z 1.8-2.0, where the 1282 A window falls at ~3700 A observed: SDSS -0.534, DESI -0.541 (agree);
  * z 2.0-2.2, where the 3025 A window falls at ~9400 A observed: SDSS -0.384, **DESI -0.104**.
  A power law with alpha_lam ~ -1.5 gives -0.56. The z 2.0-2.2 values, the DESI one above all, are biased by excess flux
  at the red end of the observed range (also visible as a x1.5-2 step at ~3400 A rest in the z 2-3 stacks).
  Per quasar (below) this turned out to be an outlier effect on the stack MEAN, not a shift of the typical spectrum.

Per quasar (v2 spectra, `09_norm_table.py` -> `quasar/v2/qso_norm_table.parquet`, 1,316,874 quasars;
`10_ratio_vs_mag.py` -> `checks/qso_ratio_vs_mag.{png,csv}` and the pooled lookup `checks/qso_ratio_lookup.csv`).
Median per 0.5 mag bin of M_AB(3025) (every bin has >= 100 quasars, so 0.5 mag throughout):

| M_AB(3025) | log R1 | n | | M_AB(3025) | log R2 | n |
|---|---|---|---|---|---|---|
| -25.25 | +0.478 | 1,542 | | -27.25 | -0.499 | 1,414 |
| -24.25 | +0.478 | 7,066 | | -26.25 | -0.522 | 9,279 |
| -23.25 | +0.437 | 23,289 | | -25.25 | -0.530 | 32,799 |
| -22.25 | +0.353 | 36,597 | | -24.75 | -0.535 | 44,000 |
| -21.25 | +0.234 | 20,470 | | -24.25 | -0.550 | 43,202 |
| -20.25 | +0.040 | 9,867 | | -23.75 | -0.562 | 19,742 |
| -19.25 | -0.247 | 8,328 | | -23.25 | -0.581 | 3,952 |
| -18.25 | -0.469 | 6,189 | | | | |
| -17.25 | -0.754 | 1,540 | | | | |

(full table in the csv; the 16-84 % spread is about +-0.1 dex.)

* **R1 depends strongly on luminosity**: flat at +0.48 brighter than M = -24, then falling by 1.2 dex to M = -17
  (host galaxy at 5115 A). At fixed magnitude, all three surveys and all three z ranges agree to <~0.05 dex. The survey
  differences seen at fixed z were population differences: DESI reaches fainter.
* **R2 is nearly constant**: -0.50 to -0.58 over 4 mag, the same on both sides of z = 2.0 and in all surveys.
* **Lesson: medians of individual ratios, not stack means.** The stack-level R2 at z 2.0-2.2 (DESI -0.10) was pulled up
  by a minority of spectra with red-end excess. The per-quasar median there is -0.536. Stack combination needs
  clipping, or medians, at the red and blue ends.

## 5. For Galaxy and Stellar

* Galaxies: same steps with your windows. Host light makes the ratios depend on luminosity and type, so bin the ratio by
  your classification axes as well as by magnitude.
* Stars: the rest frame is the observed frame. Ratios between windows are colours, set by temperature, so measure them
  per spectral-type cell and check them per survey in the same way. The red-end flux excess seen in DESI quasars at
  ~9400 A observed should show up directly in your stars if it is instrumental.
* Please record NORM / NORM_SRC / per-window means for every spectrum, as in step 2. That keeps the factors needed to
  restore absolute flux, or to re-chain later.
