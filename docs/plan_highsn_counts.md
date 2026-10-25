# Plan: how many high-S/N spectra do we have?

Date: 2026-10-22. Status: **draft, awaiting Nao's OK before running anything.**
Feeds Step 4 of `docs/pipeline_steps.md` (high-S/N selection → PCA → eigenspectra).

## Goal

A count table: survey × class × S/N threshold (and z bin), for spectra that are
(a) in the catalogs and (b) actually on disk on this Mac, so we can choose the
Step 4 S/N cut knowing the sample size it gives.

## Data

| Survey | Catalog | Rows | Spectra on disk |
|---|---|---|---|
| SDSS DR17 (BOSS/eBOSS, v5_13_2) | `~/sdssredux/spAll-v5_13_2.fits` (16 GB) | 3,958,000 | 4,131 spPlate files (3,853 plates) in `~/sdssredux/v5_13_2/` |
| SDSS DR8 (legacy) | `~/sdssredux/specObj-dr8.fits` | ~1.8 M | 3 spPlate files (plate 0416) only |
| DESI DR2 loa | `~/desiredux/zall-pix-loa.fits` (51 GB) | 64,100,684 | all 34,508 main dark+bright healpix coadds (exdisk2 + My Book) |

## S/N definition — the main issue

- **SDSS:** spAll/specObj give `SN_MEDIAN_ALL` = median S/N per pixel over the
  whole spectrum (also `SN_MEDIAN` per ugriz band). Ready to use.
- **DESI:** zall has **no per-spectrum median S/N**. It has `TSNR2_*`
  (template S/N² for a reference target of each class: an effective-depth
  measure, not this object's S/N), `DELTACHI2` (redshift confidence), and
  fiber fluxes. Three options:
  1. **Proxy (fast, minutes):** predict S/N from FIBERFLUX_R (and TSNR2 /
     COADD_EXPTIME), calibrated on a sample of real coadds (~200 files, reading
     only the FLUX/IVAR HDUs, ~150 MB each → ~30 GB, < 1 h).
  2. **Exact (slow, ~10–15 h):** compute median S/N for every spectrum from all
     34,508 coadds, reading FLUX/IVAR/MASK only (~22 % of 24 TB ≈ 5 TB from the
     USB disks). Steps 1–3 will need this pass anyway, so its output (a per-target
     S/N table) is reusable.
  3. Both: proxy now for a first answer, exact pass overnight.
  **Recommended: 3.**
- **Comparability:** SDSS pixels are 69 km/s (≈1.15 Å at 5000 Å); DESI pixels
  are 0.8 Å. S/N per pixel is not comparable across surveys, so the table will
  also give S/N per Å (S/N_pix / √Δλ). Proposed common measure: median S/N per Å
  over observed 4000–8500 Å, the range both surveys cover well.

## Quality cuts (baseline; each one's effect reported)

- SDSS: SPECPRIMARY == 1, PLATEQUALITY == 'good', ZWARNING == 0; class from CLASS.
- DESI: SURVEY == 'main', PROGRAM in {dark, bright} (what is on disk),
  ZCAT_PRIMARY, ZWARN == 0, COADD_FIBERSTATUS == 0; class from SPECTYPE.

## Counts produced

- S/N thresholds: median S/N ≥ 5, 10, 20, 30, 50 (per pixel and per Å).
- Split by class (GALAXY / QSO / STAR) and by z bin (galaxies 0–0.3–0.6–1.0–1.6;
  QSOs 0–1–2–3–5), since the PCA will probably be per class and z range.
- Two columns per cell: in catalog, and on disk here.

## Execution

1. Read the needed columns from spAll, specObj-dr8 and zall once, and cache them
   as Parquet in the scratchpad (zall is a full 51 GB read, ~10 min).
2. SDSS counts directly; mark which plate-MJDs have a local spPlate.
3. DESI: proxy calibration sample → proxy S/N for all rows → counts.
4. Write `csvfiles/highsn_counts.csv` and a short summary added to this doc.
5. (Option 2/3) start the exact DESI S/N pass in the background, writing
   `data/desi_loa_main_snr.parquet` (TARGETID, per-band and overall median S/N);
   recount when it finishes.

Nothing is modified on the data disks; outputs go only to this repo and the scratchpad.

## Questions for Nao

1. DESI S/N: option 1, 2 or 3 (recommended 3)?
2. Is "median S/N per Å over 4000–8500 Å" the right common definition, or
   would you prefer per-pixel, or S/N in a specific band/continuum window?
3. Count SDSS at catalog level too (all 3.96 M + DR8), or only what is on disk?
   Getting more SDSS spPlates would be a separate download.
