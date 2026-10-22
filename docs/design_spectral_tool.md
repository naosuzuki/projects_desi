# Design: Spectral Classification & Redshifting with Proper Errors

**Working title:** `zport` (name TBD)
**Status:** draft v0.1 — 2026-10-22
**Vision:** Take the good parts of [Redrock](https://github.com/desihub/redrock)
(DESI's production redshift fitter) and [Marz](https://github.com/samreay/Marz)
(OzDES's manual/automatic redshifting tool) and build a tool that does spectral
classification and redshift measurement **with statistically proper errors**,
for both **DESI and SDSS** spectra.

---

## 1. What we take from each parent

### From Redrock (BSD-3-Clause, Python)
Redrock is the physics-grounded forward-modeling approach:

- **χ² template fitting in flux space.** PCA (galaxy: 10 eigenvectors from
  stellar-population-synthesis + emission-line models; QSO: 4 eigenvectors;
  stars: 5 per subtype B/A/F/G/K/M/WD) fit at each point of a
  class-dependent redshift grid; linear coefficients solved analytically.
- **Instrument forward model.** Per-camera resolution matrix `R` applied to
  the redshifted template before comparison — correct treatment of the
  line-spread function, essential for narrow lines.
- **Inverse-variance weighting** with proper mask-bit handling.
- **Archetype mode** (Anand et al. 2024, [arXiv:2405.19288](https://arxiv.org/abs/2405.19288)):
  physical spectra + Legendre polynomials, bounded least squares — avoids
  unphysical PCA combinations (e.g. negative emission lines).
- **ZWARN flag logic** and **DELTACHI2** (χ² gap to the second minimum) as a
  confidence statistic.
- **Multi-exposure fitting:** fits all exposures/cameras simultaneously
  rather than a merged coadd.

**Known weakness (our main target):** `ZERR` comes from a parabola fit to the
χ²(z) minimum — a Gaussian approximation at a single mode. It ignores
multi-modality (line confusion), template mismatch, and correlated residuals;
DESI validates empirically that ZERR is underestimated for some classes.
Classification is a hard arg-min over classes, with no probability attached.

### From Marz (MIT, JavaScript — Hinton et al. 2016, [arXiv:1603.09438](https://arxiv.org/abs/1603.09438))
Marz is the robust, fast, human-in-the-loop approach:

- **Preprocessing pipeline that makes real data tractable:** cosmic-ray /
  outlier sigma-clipping, continuum subtraction (polynomial + smoothing),
  tapering — the fit sees features, not flux calibration errors.
- **FFT cross-correlation** (modified Autoz; Tonry & Davis 1979 heritage):
  orders of magnitude faster than a χ² scan — ideal as a *first pass* to
  localize candidate redshifts.
- **Combining peaks across all templates** into a ranked candidate list,
  not just a single best fit.
- **The inspection UI** — the best part of Marz: keyboard-driven visual
  inspection, overlay of template at candidate z, marked spectral features,
  cycling through the top-N automatic results, per-object comments,
  **QOP quality flags**, versioned results, and **merging of results from
  multiple inspectors**.

**Known weakness:** no formal redshift errors at all; QOP is a human/heuristic
confidence, not a statistical one.

---

## 2. The new part: proper errors

The deliverable that neither parent provides: a **posterior over (class, z)**
per object, not a point estimate with a curvature error.

1. **Likelihood from the χ² scan.** For each class *c* and grid redshift *z*,
   the template coefficients are linear → marginalize them analytically
   (Gaussian integral), giving a marginal likelihood L(z, c) rather than a
   profile χ² at the best coefficients. This is cheap: it reuses the same
   matrices Redrock already builds.
2. **Full redshift PDF.** P(z | c, data) ∝ L(z, c)·π(z, c) over the whole
   grid, refined around each mode. Report:
   - MAP redshift and **asymmetric credible intervals** (68/95%),
   - **multi-modality explicitly**: all modes with their posterior mass
     (this replaces DELTACHI2 — "P(second mode) = 3%" is interpretable;
     ΔΧ²=25 is not),
   - the full PDF saved compactly (quantiles or samples) for downstream use.
3. **Classification probabilities.** P(c | data) by marginalizing z within
   each class, with class priors (optionally magnitude/color-dependent).
   Output `P(GALAXY), P(QSO), P(STAR)` instead of a hard SPECTYPE.
4. **Empirical calibration — the crucial step.** Analytic posteriors are
   only as good as the noise model and templates. Calibrate with:
   - **repeat observations**: DESI repeats and SDSS duplicates give the
     true error distribution of Δz/(σ1⊕σ2); fit a per-class error-scaling
     (or a rank/quantile recalibration of the PDF),
   - **truth tables**: deep VI catalogs (DESI EDR VI, OzDES QOP≥4),
   - report both *raw* and *calibrated* errors.
5. **Honest failure modes.** Keep a ZWARN-like bitmask, but driven by the
   posterior: NO_DOMINANT_MODE, CLASS_AMBIGUOUS, EDGE_OF_GRID, BAD_SKY_OVERLAP.

## 3. Two-stage fitting engine

```
stage 1 (Marz-style, fast):   continuum-subtract → FFT cross-correlation
                              against all templates → top-N candidate
                              (z, class) regions per object
stage 2 (Redrock-style, exact): forward-model χ²/likelihood scan only in
                              windows around candidates + a sparse global
                              grid (to catch what stage 1 missed) →
                              posterior construction
```

This gives Marz's speed with Redrock's correctness, and the sparse global
scan bounds the risk of stage-1 misses. A `--full-scan` mode reproduces the
dense Redrock grid for validation.

## 4. Data layer: DESI + SDSS

One internal representation; two (extensible) readers.

```python
@dataclass
class Spectrum:            # one object, possibly many bands/exposures
    bands: list[Band]      # per camera/arm
    meta: dict             # targetid/plate-mjd-fiber, ra, dec, photometry

@dataclass
class Band:
    wave: np.ndarray       # vacuum, barycentric-corrected [Å]
    flux: np.ndarray
    ivar: np.ndarray
    mask: np.ndarray
    R: scipy.sparse | None # resolution matrix (DESI) or
    wdisp: np.ndarray|None # Gaussian dispersion (SDSS)
```

- **DESI reader:** healpix `coadd-*.fits` / `spectra-*.fits` (b/r/z cameras,
  resolution matrix, DESI mask bits); works against DR2/loa files already
  used by `desilib/desi_db.py`.
- **SDSS reader: built on the user's `sdsslib`** (`~/github/sdsslib`),
  developed for SDSS-III/IV. Assets we reuse directly:
  - `pylib/sdss_db.py` → `SDSSspec`: spPlate reader for DR17 (`v5_13_2`)
    and DR8 (plate>3509 switch), log10-λ grid (`coeff0/coeff1`, wid pixel
    IDs), flux/ivar/mask, SFD + CCM89 Milky-Way dereddening, deredshift +
    log rebinning.
  - `pylib/sdss_catalog.py` → `spall` / `DR8specobj`: catalog readers with
    Z, Z_ERR, ZWARNING, CLASS/SUBCLASS, **THING_ID** (duplicate spectra →
    error calibration pairs), ELODIE stellar parameters (Teff/logg/[Fe/H]
    → stellar truth table), photometry, S/N.
  - `pylib/astrolib.py`: log/velocity rebinning with proper flux+ivar+mask
    propagation (`exec_logrebinning_fluxivarmask`, dlog=1e-4 ≈ 69 km/s
    pixels) — the resampling core for stage 1.
  - `create_2dspec` / `SDSSspec2d`: **rectangular 2D flux/ivar/mask stacks
    on a common wavelength grid** (observed or rest frame) — adopted as the
    tool's batch container format for fitting, composites, and PCA template
    building; extend the same format to DESI spectra.
  - `pylib/cosmology.py`: distance modulus utilities for rest-frame work.
  SDSS wavelengths are vacuum, heliocentric; harmonize frames in the reader
  so the engine never cares. Note: spPlate/spAll files live at
  NERSC/LBNL (`$SPECTRO_REDUX`), not on this machine — for local testing,
  fetch a handful of per-object `spec-PLATE-MJD-FIBER.fits` "lite" files
  from the SAS.
- **Templates:** start from redrock's template files (`rrtemplate-*.fits`)
  and redrock-archetypes; same format so templates stay interchangeable.
  Later, build our own templates from sdsslib's rest-frame 2D stacks +
  PCA machinery (composites per class/subclass).

## 5. Package layout (proposed)

```
zport/
  io/         desi.py, sdss.py, templates.py, results.py
  prep/       continuum.py, masking.py, resample.py
  engine/     xcorr.py (stage 1), likelihood.py (stage 2, R-matrix fwd model)
  posterior/  pdf.py, classify.py, calibrate.py, flags.py
  viewer/     (phase 4) web inspection UI, Marz-inspired
  scripts/    zport-run, zport-inspect, zport-calibrate
  tests/
```

Licensing: BSD-3 (redrock) and MIT (Marz) are both permissive and
compatible — we can port/adapt code with attribution. Suggest BSD-3 for the
new tool.

## 6. Roadmap

| Phase | Deliverable | Validates against |
|------|--------------|-------------------|
| 0 | `Spectrum` + readers: SDSS via `sdsslib` wrap, DESI coadd; quick-look plotting | visual |
| 1 | Stage-2 engine with redrock templates; reproduce redrock Z/ZERR/SPECTYPE | `zall-pix` (COSMOS/SXDS extracts already in `data/`) |
| 2 | Analytic posteriors: z-PDF, credible intervals, P(class) | internal consistency |
| 3 | Stage-1 xcorr accelerator + two-stage pipeline | phase-1 results, timing |
| 4 | Error calibration from DESI repeats + SDSS duplicates | Δz pull distributions |
| 5 | Inspection viewer (QOP workflow, multi-inspector merge) | user |

Phase 1 is the keystone: once we match redrock on the user's existing
COSMOS/SXDS zall extracts, everything after is additive.

## 7. Open questions

- New repo vs. subpackage of `projects_desi`? (Suggest a new repo once
  phase 0 works; prototype here under `zport/`.)
- Viewer technology: reuse/port Marz's JS front end, or Python-native
  (bokeh/solara) to keep one language?
- GPU support: redrock has CuPy paths; defer until profiling says needed?
- Template upgrade path: NMF galaxy templates / newer QSO templates vs.
  frozen redrock DR2 templates for reproducibility.

## References

- Redrock: https://github.com/desihub/redrock (Bailey et al., in prep.)
- Archetypes: Anand et al. 2024, https://arxiv.org/abs/2405.19288
- DESI pipeline: Guy et al. 2023, https://arxiv.org/abs/2209.14482
- Marz: Hinton et al. 2016, A&C 15, 61, https://arxiv.org/abs/1603.09438
- Autoz: Baldry et al. 2014; cross-correlation: Tonry & Davis 1979
