#!/usr/bin/env python3
"""Spectrophotometric quasar spectra and their eigencoefficients (Nao 2026-11-07; ms.tex sec:qsocoeff).
Every pass-Q quasar after the quality cut (no S/N cut): DESI ZCAT_PRIMARY & ZWARN 0 & fiber ok, SDSS ZWARNING 0 on a
good plate, 0 < z < 5.  Per spectrum:
  1. photometric correction with a smooth curve: synthetic AB mags of the observed (reddened, uncorrected) spectrum vs the
     broad-band photometry (DESI: Legacy total FLUX_G/R/Z by PHOTSYS; DR17: spAll PSFFLUX g r i; DR8: photoPlate PSFFLUX
     g r i); a band counts when the spectrum covers 95 % of it.  dm_b = m_phot - m_syn fitted by a0 + a1 x,
     x = log10(lambda_eff / 6000), weights 1 / (s_phot^2 + s_syn^2 + 0.01^2); spectrum x 10^(-0.4 (a0 + a1 x)), x held
     at the outermost effective wavelengths; the tilt a1 is kept only when |a1| > 2 sigma(a1) and |a1| <= 3 mag/dex, else the
     weighted mean gray offset alone.  All |dm_b| < 0.032 mag (3 %): no correction (PCORR_MODE 0; 1 gray only; 2 gray + tilt; -1 no band; 3 SDSS extended source (OBJC_TYPE != 6), dm vs PSF recorded, no correction).
     DESI extended sources (MORPHTYPE != PSF): spectrum x MEAN_PSF_TO_FIBER_SPECFLUX vs FIBERFLUX (PCORR_REF F), point
     sources vs the total FLUX (PCORR_REF T).
  2. Milky Way extinction: F99, R_V 3.1, 0.86 x SFD E(B-V), on the corrected spectrum.
  3. rest frame (exact sub-pixel resampling of the v2 walk) and the general basis (v3) fitted redward of Ly-alpha with the
     trough masking of 16_eigenfit.py; AI / BI measured against the model; BAL = BI(C IV) > 0 or BI(Si IV) > 0 or
     AI(Mg II) >= 1000 or AI(Al III) >= 1000 (the v7 rule).  BALs are dropped from the coefficient table (kept in a BAL
     table with the general-basis fit).
  4. non-BALs: the non-BAL basis (v4_nonbal) fitted the same way to F_lambda (1+z) [cgs].  e0 = 1 and e_k = 0 in the
     pivot window, so c0 is the model flux at 3000-3050 A = the eigenspectra normalization: NORM = c0, C_k = c_k / c0.
Output: one parquet per shard in PROD/specphot/<survey>/, --merge -> PROD/specphot/qso_specphot_coeff_<tag>.parquet
(non-BAL) and qso_specphot_bal_<tag>.parquet (BAL).  Inputs: the DESI native-pixel store on exdisk3 (no coadd read), the
SDSS spPlate files, the photometry caches PROD/specphot/phot_*.parquet (scratch phot_cache.py, 2026-11-07).
Usage:
  python3 21_specphot_fit.py --survey desi --shard 0/3 [--nfiles 3]     # also dr17, dr8 (--shard 0/1)
  python3 21_specphot_fit.py --merge
"""
import os, sys, glob, time, argparse, warnings, importlib.util
import numpy as np, pandas as pd, fitsio
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'programs_templates'))
import qsolib as Q
import restframe as RF
_s = importlib.util.spec_from_file_location('eigenfit', os.path.join(HERE, '16_eigenfit.py')); EF = importlib.util.module_from_spec(_s); _s.loader.exec_module(EF)
passB = Q.passB

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{PROD}/specphot'
NATIVE = '/Volumes/exdisk3/sdssdesi_native/desi'
SPECTRO = os.path.expanduser('~/sdssredux'); PLATE_DIRS = {'dr17': [f'{SPECTRO}/v5_13_2'], 'dr8': [f'{SPECTRO}/DR8/26', f'{SPECTRO}/DR8/103', f'{SPECTRO}/DR8/104']}
BASIS_GEN = os.path.join(HERE, '..', 'data_eigenspectra', 'quasar_eigenspectra_v3.fits')
BASIS_NB = os.path.join(HERE, '..', 'data_eigenspectra', 'quasar_eigenspectra_v4_nonbal.fits')
NVEC = 10; NITER = 2; TOL = 0.032; FLOOR = 0.01; MIN_RANGE = 0.95; X0 = 6000.0; TILT_NSIG = 2.0; TILT_MAX = 3.0
STATUS = os.path.expanduser('~/data/sdssdesi_status')
BANDS = {'desi': 'grz', 'dr17': 'gri', 'dr8': 'gri'}


def filters(survey, south=True):
    if survey == 'desi':
        return {'g': 'decam2014-g', 'r': 'decam2014-r', 'z': 'decam2014-z'} if south else {'g': 'BASS-g', 'r': 'BASS-r', 'z': 'MzLS-z'}
    return {b: f'sdss2010-{b}' for b in 'gri'}


def synth(w, flux, ivar, good, names, scale=None):
    """AB mags and errors of the spectrum (cgs) in the named filters; the red tail beyond the spectrum is padded with
    the mean of the last 250 good pixels (as pass B for the z band)."""
    n = len(flux); out = {}
    if scale is not None:
        flux = flux * scale[:, None]; ivar = ivar / scale[:, None]**2
    idx = np.arange(w.size); first = np.where(good.any(axis=1), np.argmax(good, axis=1), w.size); last = np.where(good.any(axis=1), w.size - 1 - np.argmax(good[:, ::-1], axis=1), -1)
    inr = (idx[None, :] >= first[:, None]) & (idx[None, :] <= last[:, None])      # the spectrum's coverage envelope
    hi = max(passB.band_range(passB.filt(nm))[1] for nm in names.values())
    if hi > w[-1]:
        dw = w[-1] - w[-2]; wpad = np.arange(w[-1] + dw, hi + 3 * dw, dw); wz = np.concatenate([w, wpad])
        padval = np.nanmean(np.where(good[:, -250:], flux[:, -250:], np.nan), axis=1)
        fz = np.concatenate([flux, np.repeat(np.nan_to_num(padval)[:, None], len(wpad), axis=1)], axis=1)
        gz = np.concatenate([good, np.ones((n, len(wpad)), bool)], axis=1); inz = np.concatenate([inr, np.ones((n, len(wpad)), bool)], axis=1)
        ivz = np.concatenate([ivar, np.zeros((n, len(wpad)))], axis=1)
    for b, nm in names.items():
        f = passB.filt(nm); lo, bhi = passB.band_range(f)
        if bhi > w[-1]:
            m, e = passB.ab_mags(f, wz, fz, gz, inrange2d=inz, ivar2d=ivz, min_range=MIN_RANGE)
        else:
            m, e = passB.ab_mags(f, w, flux, good, inrange2d=inr, ivar2d=ivar, min_range=MIN_RANGE)
        if lo < w[0]:
            m[:] = np.nan
        out[b] = (m, e, f.effective_wavelength.value)
    return out


def correction(w, syn, pflux, pivar):
    """Smooth photometric correction per spectrum. syn: {band: (m_syn, e_syn, lam_eff)}; pflux / pivar: {band: nanomaggies}.
    Returns factor [n, npix], a0, a1, mode, nband, {band: dm}."""
    bands = list(syn); n = len(next(iter(syn.values()))[0])
    dm = np.full((n, len(bands)), np.nan); wt = np.zeros((n, len(bands))); xb = np.zeros(len(bands))
    for j, b in enumerate(bands):
        m, e, leff = syn[b]; xb[j] = np.log10(leff / X0)
        with np.errstate(all='ignore'):
            ok = np.isfinite(m) & (pflux[b] > 0) & np.isfinite(pflux[b])
            mp = np.where(ok, 22.5 - 2.5 * np.log10(np.where(ok, pflux[b], 1.0)), np.nan)
            ep = np.where(ok & (pivar[b] > 0), 2.5 / np.log(10) / (np.where(ok, pflux[b], 1.0) * np.sqrt(np.where(pivar[b] > 0, pivar[b], 1.0))), 0.05)
            dm[:, j] = mp - m; wt[:, j] = np.where(ok, 1.0 / (ep**2 + np.nan_to_num(e, nan=0.0)**2 + FLOOR**2), 0.0)
    wt[~np.isfinite(dm)] = 0.0; dmz = np.nan_to_num(dm)
    nb = (wt > 0).sum(axis=1)
    S0 = wt.sum(axis=1); Sx = wt @ xb; Sxx = wt @ xb**2; Sy = (wt * dmz).sum(axis=1); Sxy = (wt * dmz) @ xb
    with np.errstate(all='ignore'):
        det = S0 * Sxx - Sx**2
        a1 = np.where(nb >= 2, (S0 * Sxy - Sx * Sy) / det, 0.0); a0 = np.where(nb >= 2, (Sxx * Sy - Sx * Sxy) / det, np.where(nb == 1, Sy / S0, 0.0))
        a1 = np.nan_to_num(a1); a0 = np.nan_to_num(a0)
        sig_a1 = np.sqrt(np.where(nb >= 2, S0 / det, np.inf))              # error of the tilt from the weights
        tilt = (nb >= 2) & (np.abs(a1) > TILT_NSIG * sig_a1) & (np.abs(a1) <= TILT_MAX)   # the tilt is applied only when significant and physical
        a0 = np.where(tilt, a0, np.where(nb > 0, Sy / S0, 0.0)); a1 = np.where(tilt, a1, 0.0)
    within = nb > 0
    within &= (np.abs(np.where(wt > 0, dm, 0.0)) < TOL).all(axis=1)
    mode = np.where(nb == 0, -1, np.where(within, 0, np.where(tilt, 2, 1))).astype(np.int8)
    a0 = np.where(mode > 0, a0, 0.0); a1 = np.where(mode > 0, a1, 0.0)
    xl = np.log10(w / X0)
    xlo = np.where(tilt, np.min(np.where(wt > 0, xb[None, :], np.inf), axis=1), 0.0); xhi = np.where(tilt, np.max(np.where(wt > 0, xb[None, :], -np.inf), axis=1), 0.0)
    xc = np.clip(xl[None, :], xlo[:, None], xhi[:, None])
    fac = 10.0**(-0.4 * (a0[:, None] + a1[:, None] * xc))
    return fac, a0, a1, mode, nb, {b: dm[:, j] for j, b in enumerate(bands)}


def process(w, flux, ivar, mask, meta, survey, E1, ok1, E2, ok2, names_fn, pflux, pivar, sref=None, ext=None):
    """flux / ivar [n, npix] in 1e-17 cgs on the native grid; meta: DataFrame with Z, EBV. Returns a DataFrame of results."""
    n = len(meta); good = (ivar > 0) & (mask == 0)
    fcgs = flux * passB.FLUX_UNIT; ivcgs = np.where(good, ivar / passB.FLUX_UNIT**2, 0.0)
    res = meta.reset_index(drop=True).copy()
    # step 1: photometric correction (two filter systems for DESI)
    fac = np.ones_like(fcgs); a0 = np.zeros(n); a1 = np.zeros(n); mode = np.full(n, -1, np.int8); nb = np.zeros(n, int); dms = {b: np.full(n, np.nan) for b in BANDS[survey]}
    groups = [(np.ones(n, bool), names_fn(True))] if survey != 'desi' else [(meta.PHOTSYS.to_numpy() == 'S', names_fn(True)), (meta.PHOTSYS.to_numpy() != 'S', names_fn(False))]
    for rows, names in groups:
        if not rows.any():
            continue
        syn = synth(w, fcgs[rows], ivcgs[rows], good[rows], names, None if sref is None else sref[rows])
        f_, a0_, a1_, mode_, nb_, dm_ = correction(w, syn, {b: pflux[b][rows] for b in names}, {b: pivar[b][rows] for b in names})
        fac[rows] = f_; a0[rows] = a0_; a1[rows] = a1_; mode[rows] = mode_; nb[rows] = nb_
        for b in names:
            dms[b][rows] = dm_[b]
    if ext is not None:                                       # SDSS extended source: the fiber-flux scale is not calibrated here; record dm vs PSF, no correction
        mode = np.where(ext & (mode >= 0), 3, mode).astype(np.int8); a0 = np.where(mode == 3, 0.0, a0); a1 = np.where(mode == 3, 0.0, a1); fac[mode == 3] = 1.0
    res['PCORR_MODE'] = mode; res['PCORR_NBAND'] = nb.astype(np.int8); res['PCORR_A0'] = a0.astype(np.float32); res['PCORR_A1'] = a1.astype(np.float32)
    for b in BANDS[survey]:
        res[f'PCORR_DM_{b.upper()}'] = dms[b].astype(np.float32)
    fcgs *= fac; ivcgs /= fac**2
    # step 2: extinction
    A1 = passB.extinction.fitzpatrick99(np.asarray(w, float), 1.0, 3.1)
    ebv = np.nan_to_num(meta.EBV.to_numpy(float)); A = (3.1 * 0.86 * ebv)[:, None] * A1[None, :]
    fcgs *= 10**(0.4 * A); ivcgs *= 10**(-0.8 * A)
    # rest frame, F_lambda (1+z) on the quasar lattice
    z = meta.Z.to_numpy(float)
    if survey == 'desi':
        rf, riv = RF.desi_rest(w, fcgs, ivcgs, (~good).astype(np.uint8), z, Q.QRID_MIN, Q.QNREST)
    else:
        rf, riv = RF.sdss_rest(meta.attrs['coeff0'], fcgs, ivcgs, (~good).astype(np.uint8), z, Q.QRID_MIN, Q.QNREST)
    F = rf.astype(np.float64) * (1 + z)[:, None]; IV = riv.astype(np.float64) / (1 + z)[:, None]**2
    F[~np.isfinite(F)] = 0; IV[~np.isfinite(IV)] = 0; IV[~np.isfinite(F)] = 0
    cov = IV > 0; res['REST_LAMMIN'] = np.where(cov.any(axis=1), Q.QWAVE[np.argmax(cov, axis=1)], np.nan).astype(np.float32)
    res['REST_LAMMAX'] = np.where(cov.any(axis=1), Q.QWAVE[Q.QNREST - 1 - np.argmax(cov[:, ::-1], axis=1)], np.nan).astype(np.float32)
    # step 3: general basis, BAL measurement
    c, chi2, npix, nclip, bal = EF.fit_file(E1, ok1, F, IV, NITER)
    for k in range(E1.shape[0]):
        res[f'G{k}'] = c[:, k].astype(np.float32)
    res['G_CHI2'] = chi2.astype(np.float32); res['G_NPIX'] = npix.astype(np.int32); res['G_NCLIP'] = nclip.astype(np.int32)
    for k, v in bal.items():
        if k in EF.BAL_COLS and k not in ('REST_LAMMIN', 'REST_LAMMAX'):
            res[k] = v
    g = lambda col, t=0: np.nan_to_num(res[col].to_numpy(float)) > t
    res['BAL'] = (g('BAL_BI_CIV') | g('BAL_BI_SIIV') | g('BAL_AI_MGII', 999) | g('BAL_AI_ALIII', 999)).astype(np.int8)
    # step 4: non-BAL basis on the non-BALs
    for k in range(1, E2.shape[0]):
        res[f'C{k}'] = np.float32(np.nan)
    res['NORM'] = np.float64(np.nan); res['CHI2'] = np.float32(np.nan); res['NPIX'] = np.int32(0); res['NCLIP'] = np.int32(0)
    nbm = res.BAL.to_numpy() == 0
    if nbm.any():
        c2, chi2, npix, nclip, _ = EF.fit_file(E2, ok2, F[nbm], IV[nbm], NITER)
        with np.errstate(all='ignore'):
            norm = c2[:, 0]; cn = c2[:, 1:] / norm[:, None]
        res.loc[nbm, 'NORM'] = norm
        for k in range(1, E2.shape[0]):
            res.loc[nbm, f'C{k}'] = cn[:, k - 1].astype(np.float32)
        res.loc[nbm, 'CHI2'] = chi2.astype(np.float32); res.loc[nbm, 'NPIX'] = npix.astype(np.int32); res.loc[nbm, 'NCLIP'] = nclip.astype(np.int32)
    return res


def desi_sample():
    cols = ['TARGETID', 'PROGRAM', 'FILE', 'Z', 'ZWARN', 'QCUT', 'EBV', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'Z_SYS', 'Z_SYS_SRC', 'SNR_Z', 'TARGET_RA', 'TARGET_DEC']
    d = pd.read_parquet(f'{PROD}/qso_master_desi.parquet', columns=cols); d = d[d.QCUT].drop(columns='QCUT')
    ph = pd.read_parquet(f'{OUT}/phot_desi_zall.parquet', columns=['TARGETID', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'MORPHTYPE', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z']).drop_duplicates('TARGETID')
    d = d.merge(ph, on='TARGETID', how='left')
    return d


def sdss_sample(survey):
    cols = ['PLATE', 'MJD', 'FIBERID', 'SURVEY', 'Z', 'ZWARN', 'PLATEQUALITY', 'EBV', 'Z_SYS', 'Z_SYS_SRC', 'SNR_ALL', 'TARGET_RA', 'TARGET_DEC', 'THING_ID']
    d = pd.read_parquet(f'{PROD}/qso_master_sdss.parquet', columns=cols); d = d[(d.SURVEY == survey) & (d.ZWARN == 0) & (d.PLATEQUALITY == 'good') & (d.Z > 0) & (d.Z < 5)]
    ph = pd.read_parquet(f'{OUT}/phot_{"dr17_spall" if survey == "dr17" else "dr8_photoplate"}.parquet')
    keep = ['PLATE', 'MJD', 'FIBERID', 'OBJC_TYPE'] + [f'PSFFLUX_{b}' for b in 'GRI'] + [f'PSFFLUX_IVAR_{b}' for b in 'GRI']
    d = d.merge(ph[keep].drop_duplicates(['PLATE', 'MJD', 'FIBERID']), on=['PLATE', 'MJD', 'FIBERID'], how='left')
    return d


def run_desi(a, E1, ok1, E2, ok2):
    d = desi_sample(); files = sorted(d.FILE.unique()); i, m = map(int, a.shard.split('/')); files = files[i::m]
    if a.pattern:
        files = [f for f in files if a.pattern in f]
    if a.nfiles:
        files = files[:a.nfiles]
    out = f'{OUT}/desi/shard{i}of{m}.parquet'; os.makedirs(f'{OUT}/desi', exist_ok=True)
    done = pd.read_parquet(out) if os.path.exists(out) else None; done_files = set(done.FILE) if done is not None else set()
    parts = [done] if done is not None else []; t0 = time.time(); nrow = 0
    byfile = {f: g for f, g in d.groupby('FILE')}
    todo = [f for f in files if f not in done_files]; print(f'desi shard {i}/{m}: {len(files)} files, {len(todo)} to do', flush=True)
    for j, fn in enumerate(todo):
        g = byfile[fn]; p = os.path.join(NATIVE, os.path.basename(fn).replace('.fits', '') + '.npz')
        if not os.path.exists(p):
            print('missing native file', p, flush=True); continue
        z = np.load(p); tid = z['TARGETID']; pos = pd.Series(np.arange(len(tid)), index=tid)
        hit = pos.reindex(g.TARGETID.to_numpy()); okr = hit.notna().to_numpy(); g = g[okr]; rows = hit[okr].to_numpy(int)
        if len(g) == 0:
            continue
        w = z['W0'] + z['DW'] * np.arange(z['FLUX'].shape[1]); flux = z['FLUX'][rows].astype(np.float64); ivar = z['IVAR'][rows].astype(np.float64)
        # point sources (MORPHTYPE PSF): spectrum vs the total flux; extended: spectrum x point-source fiber fraction vs the fiber flux
        pt = g.MORPHTYPE.fillna('').str.strip().isin(['PSF', 'GPSF']).to_numpy(); p2f = np.nan_to_num(g.MEAN_PSF_TO_FIBER_SPECFLUX.to_numpy(float), nan=1.0); p2f = np.where(p2f > 0, p2f, 1.0)
        pflux, pivar = {}, {}
        for b in 'grz':
            tot = g[f'FLUX_{b.upper()}'].to_numpy(float); fib = np.nan_to_num(g[f'FIBERFLUX_{b.upper()}'].to_numpy(float)); iv = np.nan_to_num(g[f'FLUX_IVAR_{b.upper()}'].to_numpy(float))
            with np.errstate(all='ignore'):
                pflux[b] = np.where(pt, tot, fib); pivar[b] = np.where(pt, iv, np.where(fib > 0, iv * (tot / np.where(fib > 0, fib, 1.0))**2, 0.0))
        g = g.assign(PCORR_REF=np.where(pt, 'T', 'F'))
        r = process(w, flux, ivar, np.zeros(flux.shape, np.uint8), g.drop(columns=[c for c in g.columns if c.startswith('FLUX_IVAR') or c.startswith('FIBERFLUX')]), 'desi', E1, ok1, E2, ok2, lambda s: filters('desi', s), pflux, pivar, sref=np.where(pt, 1.0, p2f))
        parts.append(r); nrow += len(r)
        if (j + 1) % a.save == 0 or j == len(todo) - 1:
            pd.concat(parts, ignore_index=True).to_parquet(out + '.tmp', index=False); os.replace(out + '.tmp', out)
            print(f'  {j + 1}/{len(todo)} files, {nrow:,} quasars, {(time.time() - t0) / 60:.1f} min', flush=True)
    print('shard finished', flush=True)


def run_sdss(a, survey, E1, ok1, E2, ok2):
    d = sdss_sample(survey); i, m = map(int, a.shard.split('/'))
    plates = sorted(set(zip(d.PLATE, d.MJD)))[i::m]
    if a.nfiles:
        plates = plates[:a.nfiles]
    out = f'{OUT}/{survey}/shard{i}of{m}.parquet'; os.makedirs(f'{OUT}/{survey}', exist_ok=True)
    done = pd.read_parquet(out) if os.path.exists(out) else None; done_pl = set(zip(done.PLATE, done.MJD)) if done is not None else set()
    parts = [done] if done is not None else []; t0 = time.time(); nrow = 0
    byp = {k: g for k, g in d.groupby(['PLATE', 'MJD'])}
    todo = [p for p in plates if p not in done_pl]; print(f'{survey} shard {i}/{m}: {len(plates)} plates, {len(todo)} to do', flush=True)
    for j, (plate, mjd) in enumerate(todo):
        path = None
        for dd in PLATE_DIRS[survey]:
            cand = glob.glob(f'{dd}/{plate:04d}/spPlate-{plate:04d}-{mjd}.fits') + glob.glob(f'{dd}/spPlate-{plate:04d}-{mjd}.fits')
            if cand:
                path = cand[0]; break
        if path is None:
            print('missing plate', plate, mjd, flush=True); continue
        g = byp[(plate, mjd)]; Fh = fitsio.FITS(path); h = Fh[0].read_header(); coeff0 = float(h['COEFF0'])
        fib = g.FIBERID.to_numpy() - 1
        flux = Fh[0].read()[fib].astype(np.float64); ivar = Fh[1].read()[fib].astype(np.float64); mask = Fh[2].read()[fib]; Fh.close()
        npix = flux.shape[1]; w = 10.0**(coeff0 + np.arange(npix) * Q.DLOG)
        bad = (mask != 0).astype(np.uint8)
        if survey == 'dr17':
            bad[:, (w >= 6322.7) & (w <= 6335.8)] = 1                   # the DR17 sky residual (ms.tex sec:step3)
        meta = g.copy(); meta.attrs['coeff0'] = coeff0
        pflux = {b: g[f'PSFFLUX_{b.upper()}'].to_numpy(float) for b in 'gri'}; pivar = {b: np.nan_to_num(g[f'PSFFLUX_IVAR_{b.upper()}'].to_numpy(float)) for b in 'gri'}
        ext = g.OBJC_TYPE.fillna(6).to_numpy(int) != 6                 # photo OBJC_TYPE 6 = point source, 3 = galaxy
        r = process(w, flux, ivar, bad, meta.drop(columns=[c for c in meta.columns if c.startswith('PSFFLUX')]), survey, E1, ok1, E2, ok2, lambda s: filters(survey), pflux, pivar, ext=ext)
        parts.append(r); nrow += len(r)
        if (j + 1) % a.save == 0 or j == len(todo) - 1:
            pd.concat(parts, ignore_index=True).to_parquet(out + '.tmp', index=False); os.replace(out + '.tmp', out)
            print(f'  {j + 1}/{len(todo)} plates, {nrow:,} quasars, {(time.time() - t0) / 60:.1f} min', flush=True)
    print('shard finished', flush=True)


def merge(a):
    parts = []
    for sv in ['desi', 'dr17', 'dr8']:
        for f in sorted(glob.glob(f'{OUT}/{sv}/shard*.parquet')):
            p = pd.read_parquet(f); p['SURVEY'] = sv
            for c in ['TARGETID', 'PLATE', 'MJD', 'FIBERID', 'THING_ID']:    # nullable ints BEFORE the concat: a float64 pass would round TARGETID > 2^53
                if c in p:
                    p[c] = p[c].astype('Int64')
            parts.append(p)
    d = pd.concat(parts, ignore_index=True)
    nb = d[d.BAL == 0]; bal = d[d.BAL == 1]
    nb.to_parquet(f'{OUT}/qso_specphot_coeff_{a.tag}.parquet', index=False); bal.to_parquet(f'{OUT}/qso_specphot_bal_{a.tag}.parquet', index=False)
    print(f'{len(d):,} quasars: {len(nb):,} non-BAL -> qso_specphot_coeff_{a.tag}.parquet, {len(bal):,} BAL -> qso_specphot_bal_{a.tag}.parquet')
    print(d.groupby('SURVEY').agg(n=('Z', 'size'), bal=('BAL', 'mean'), corrected=('PCORR_MODE', lambda s: (s > 0).mean()), noband=('PCORR_MODE', lambda s: (s < 0).mean()), fitted=('NORM', lambda s: s.notna().mean())).round(4).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', default='desi', choices=['desi', 'dr17', 'dr8']); ap.add_argument('--shard', default='0/1')
    ap.add_argument('--nfiles', type=int, default=0); ap.add_argument('--pattern', default=None, help='test: only DESI coadd files containing this'); ap.add_argument('--save', type=int, default=20); ap.add_argument('--merge', action='store_true'); ap.add_argument('--tag', default='v1')
    ap.add_argument('--basis-gen', default=BASIS_GEN); ap.add_argument('--basis-nonbal', default=BASIS_NB)
    a = ap.parse_args()
    if a.merge:
        merge(a); sys.exit()
    E1, ok1 = EF.load_basis(a.basis_gen, NVEC); E2, ok2 = EF.load_basis(a.basis_nonbal, NVEC)
    if a.survey == 'desi':
        run_desi(a, E1, ok1, E2, ok2)
    else:
        run_sdss(a, a.survey, E1, ok1, E2, ok2)
