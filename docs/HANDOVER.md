# Handover: Spectral Classification & Redshifting Tool

Date: 2026-10-22
From: Claude Code session on Mac (desktop), with Nao Suzuki
Purpose: transfer this project to the computer that hosts the SDSS/DESI data
disks, so work can continue there (human or a fresh Claude session).

---

## 1. The project in one paragraph

Build a new tool that combines the good parts of **Redrock**
(https://github.com/desihub/redrock — DESI's forward-model χ² template
fitter; BSD-3) and **Marz** (https://github.com/samreay/Marz — OzDES's fast
cross-correlation + inspection UI; MIT) to do **spectral classification and
redshifting with proper errors**, for **both DESI and SDSS spectra**. "Proper
errors" means full posteriors over (class, z) — asymmetric credible
intervals, explicit multi-modality, P(GALAXY)/P(QSO)/P(STAR) — empirically
calibrated with repeat observations, instead of Redrock's χ²-curvature ZERR
and hard SPECTYPE.

The full design (architecture, two-stage engine, data layer, roadmap, open
questions) is in **`docs/design_spectral_tool.md`** in this repo. Read that
first; this handover only adds machine/transfer specifics.

## 2. Repositories involved

| Repo | Where | Role |
|---|---|---|
| `projects_desi` (this repo) | github.com/… (user's) | design docs, DESI catalog work; prototype home for the new tool |
| `sdsslib` | `~/github/sdsslib` on the Mac (user's own, has .git) | **the SDSS foundation** — see §4 |
| `desilib/desi_db.py` | inside this repo | DESI zall-pix catalog reader (DR2 "loa") |
| desihub/redrock | GitHub | reference implementation + template files |
| samreay/Marz | GitHub | reference for xcorr stage + inspection UI |

Transfer: clone/copy `projects_desi` AND `sdsslib` to the data machine.

## 3. State as of handover

- `docs/design_spectral_tool.md` written (v0.1 draft). **No tool code exists yet.**
- Nothing committed yet on the Mac beyond prior catalog work (branch `main`;
  untracked: `csvfiles/`, `docs/`).
- Repo already contains validation catalogs: `data/zall-pix-loa-cosmos.csv`,
  `data/zall-pix-loa-sxds.csv` (DESI DR2 zall extracts; made by
  `programs/03_read_zall_extract.py`), and `csvfiles/desidr2star*.csv`.
- The Mac has **no spectra files locally** (internal disk nearly full);
  SDSS/DESI spectra live on the data machine — hence this handover.

## 4. User's existing SDSS resources (sdsslib) — reuse, don't rewrite

- `pylib/sdss_db.py` — `SDSSspec(plate,mjd,fiber)`: reads spPlate for
  **DR17 (`v5_13_2`, plate>3509)** and **DR8** (with run2d 103/104/26
  fallbacks); log10-λ grid via `coeff0/coeff1` + integer pixel IDs (`wid =
  coeff0/1e-4 + i`); SFD+CCM89 MW dereddening; `deredshift()` with
  log-rebinning. Also `create_2dspec`/`SDSSspec2d`: rectangular 2D
  flux/ivar/mask FITS stacks on a common grid (observed or rest frame) —
  **adopted as the tool's batch container format**.
- `pylib/sdss_catalog.py` — `spall` (DR17) and `DR8specobj` readers:
  Z, Z_ERR, ZWARNING, CLASS/SUBCLASS, **THING_ID** (→ repeat-observation
  pairs for error calibration), ELODIE Teff/logg/[Fe/H] (→ stellar truth),
  photometry, S/N.
- `pylib/astrolib.py` — `exec_logrebinning_fluxivarmask` (dlog=1e-4 ≈
  69 km/s) with flux+ivar+mask propagation → resampling core for the
  cross-correlation stage. Also CCM extinction.
- `pylib/cosmology.py` — distance modulus utilities.
- `programs/` — 01a/01b build spAll→CSV; 08/09 build 2D spec stacks (qsub
  scripts exist for batch runs); 02/03 are read/plot examples.

## 5. Environment the code expects (set on the data machine)

```bash
export SPECTRO_REDUX=/path/to/sdss/spectro/redux/   # contains v5_13_2/, spAll-v5_13_2.fits, specObj-dr8.fits …
export DESI_REDUX=/path/to/desi/spectro/redux/loa/  # contains zall-pix-loa.fits (+ healpix/ coadds)
export PYLIB=~/github/sdsslib/pylib                 # sdss_catalog etc.
export DESILIB=~/github/projects_desi/desilib
```
Note the trailing-slash convention: both `sdss_db.py` and `desi_db.py`
concatenate paths as `os.environ[VAR] + 'file.fits'`.
Python deps: numpy, scipy, pandas, fitsio, astropy, dustmaps (SFD maps
downloaded), extinction, matplotlib.

## 6. Data inventory needed on the data machine (verify, then record paths)

- [ ] SDSS DR17 spPlate files (`v5_13_2/PLATE/spPlate-PLATE-MJD.fits`)
- [ ] `spAll-v5_13_2.fits`; DR8: `specObj-dr8.fits`, `photoPlate/photoMatchPlate-dr8.fits`
- [ ] DESI DR2 (loa): `zall-pix-loa.fits`; healpix `coadd-*.fits` for at
      least COSMOS (RA 148–152, Dec 0–4) and SXDS (RA 33.5–38, Dec −6.2…−3)
- [ ] redrock templates (`rrtemplate-*.fits`; pip-install redrock or clone
      desihub/redrock + desihub/redrock-archetypes)

## 7. First task on the data machine: Phase 0

1. Set env vars (§5); smoke-test `sdsslib`: read one spPlate spectrum and
   plot it (`programs/02_readspec.py` is the example).
2. Smoke-test DESI: read one healpix coadd for a COSMOS target from
   `data/zall-pix-loa-cosmos.csv` (targetid → healpix via the zall catalog).
3. Build the unified data model in a new package dir (working name `zport/`,
   see design doc §5): `Spectrum`/`Band` dataclasses; SDSS backend wrapping
   `SDSSspec`; DESI backend reading coadd b/r/z + resolution matrix;
   one quick-look plot function for both.
4. Then Phase 1 (keystone): reproduce redrock Z/ZERR/SPECTYPE on the
   COSMOS/SXDS extracts before adding anything new.

## 8. Decisions already made / open

**Decided:** reuse sdsslib rather than rewrite; 2D stack as batch format;
two-stage engine (xcorr → forward-model likelihood); DESI+SDSS from day one;
calibrate errors with THING_ID / DESI repeats; start from redrock templates.
**Open (ask Nao):** tool name (`zport` is a placeholder); new repo vs.
subdir here; viewer tech (port Marz JS vs. Python-native bokeh/solara);
GPU support timing; template freeze vs. upgrade.
