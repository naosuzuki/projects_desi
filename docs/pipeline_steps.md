# Spectral preprocessing pipeline (SDSS + DESI)

Date: 2026-10-22. Steps as defined by Nao; notes on local data from the
inventory done the same day on this Mac (dhcp-124.hep.fsu.edu), which holds both
datasets (see `docs/HANDOVER.md`).

## Step 1 — Spectrophotometric calibration

Make every SDSS and DESI spectrum spectrophotometric by tying it to broad-band
photometry: Pan-STARRS (grizy), DESI Legacy Survey (g r z), and SDSS for the
u band.

Photometry available locally:
- SDSS spectra: `spAll-v5_13_2.fits` / `specObj-dr8.fits` carry SDSS ugriz
  (PSFMAG, FIBER2MAG, CMODELMAG, EXTINCTION) for every spectrum, so u is already here.
- DESI spectra: `zall-pix-loa.fits` (64.1 M rows) carries LS photometry for
  every target: FLUX_G/R/Z (+ IVAR), FIBERFLUX_*, FIBERTOTFLUX_*, W1/W2,
  PHOTSYS (N/S), Gaia G/BP/RP, MEAN_PSF_TO_FIBER_SPECFLUX.
- Pan-STARRS DR2 and LS DR10 catalogs exist locally only for COSMOS
  (`/Volumes/My Book/data/catalog/{PanSTARRS,DESI_Legacy,SDSS}/COSMOS/`).
  All-sky use needs a crossmatch (PS1 DR2 via MAST CasJobs; SDSS u for DESI
  targets via SDSS photoObj; LS for SDSS targets via the LS DR10 tractor/sweeps).

Open: fiber vs. total magnitudes as the reference (aperture correction for
extended sources); form of the correction (scale-only vs. low-order polynomial
in λ); handling of the SDSS/PS1/LS filter-system offsets.

## Step 2 — Milky Way extinction correction

- `dustmaps` + SFD maps work (`~/.dustmapsrc` → `/Volumes/exdisk1/data/photoz/dustmaps`).
- `zall` has an EBV column per target; `sdsslib/pylib/sdss_db.py` already does SFD + CCM89.
Open: extinction law (CCM89 vs. Fitzpatrick99) and whether to use the
Schlafly & Finkbeiner (2011) 0.86 rescaling of SFD.

## Step 3 — De-redshift extragalactic spectra with known spec-z

- SDSS: Z, Z_ERR, ZWARNING, CLASS from spAll / specObj.
- DESI: redrock Z, ZERR, ZWARN, SPECTYPE from zall-pix-loa.
- `SDSSspec.deredshift()` already does rest-frame log rebinning (dlog = 1e-4,
  ≈ 69 km/s) with flux/ivar/mask propagation.
Open: which spec-z counts as "known" (ZWARNING == 0 only? which CLASS?).

## Step 4 — High-S/N selection, PCA, eigenspectra

Select high-S/N rest-frame spectra from Step 3 and run PCA to obtain eigenspectra.

- Input container: sdsslib's 2D stacks (`create_2dspec` / `SDSSspec2d`) already
  put flux/ivar/mask on a common rest-frame log-λ grid — a natural PCA input.
- S/N measures available: SDSS SN_MEDIAN_ALL (spAll); DESI TSNR2_* in zall,
  or a median S/N computed from the coadd ivar.
- De-redshifted spectra cover different rest-frame windows and have masked
  pixels, so the PCA must handle missing data and per-pixel weights
  (weighted/EM PCA, e.g. Bailey 2012 `empca`, or gappy PCA), not plain SVD.
Open: S/N cut; separate PCA per class (GALAXY / QSO / STAR, maybe redshift
bins); rest-frame λ range; normalization before PCA; number of eigenspectra
kept; one joint SDSS+DESI basis or one per survey (resolution and λ coverage differ).
