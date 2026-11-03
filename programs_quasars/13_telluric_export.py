#!/usr/bin/env python3
"""Phase 0 of the telluric plan (docs/PLAN_telluric.md, Stellar thread's request 2026-11-01): observed-frame residuals of
high-S/N individual quasars for the joint telluric EM-PCA. Per quasar (QCUT, S/N > 20 per SDSS-size pixel):
  model = fit of its rest spectrum (the thread's calibration and exact sub-pixel rest frame, as the stacks) with the v3
          pivoted basis (06_empca.py: e0 + nvec vectors), weighted by ivar, with the observed telluric windows, the Lya forest
          (rest < 1216 A) and the broad lines (+-5000 km/s) at zero weight; one 5-sigma clipping pass (absorbers, cosmics);
  R = f / model - 1 on the OBSERVED native pixels (the model is evaluated at lambda_obs / (1 + z) from the rest fit, so no
      resampling of the data: SDSS spPlates sit exactly on the common lattice, DESI's 0.8 A pixels are interpolated to it);
  W = ivar x model^2, 0 where masked (forest, broad lines, bad pixels) or outside the basis.
Common grid: log10 lambda = 3.5502 + 1e-4 i, i = 0..4446 (3549-9895 A). Metadata: ID, CLASS (slope class SBIN 0-5),
EWCLASS (EBIN), SURVEY, Z (systemic), SNR, PLATE (SDSS plate / DESI dominant TILEID by EXPTIME), MJD (plate MJD / EXPTIME-
weighted mean), AIRMASS (NaN), PROGRAM; DESI also the per-exposure list EXP_OFFSET [n+1], EXP_EXPID, EXP_TILEID, EXP_MJD,
EXP_EXPTIME (expfibermap rows in the coadd, any camera). FLUX_SCALE (Stellar's request): median of the (calibrated,
dereddened) model over 6000-9000 A observed, in 1e-17 erg/s/cm2/A, to separate additive sky from multiplicative telluric
residuals; EBV: SFD E(B-V) of the object (Galactic Ca II K/H and Na D are visible in the mean residual).
Output: /Volumes/exdisk3/sdssdesi_telluric/phase0/quasar_<desi|sdss>.npz (exdisk3 only; nothing else is written).
Usage: python3 13_telluric_export.py --survey sdss            (all 18.6k S/N > 20 quasars; spPlates on the internal disk)
       python3 13_telluric_export.py --survey desi --nmax 20000   (native store on exdisk3; files with the most S/N > 20 quasars first)
"""
import os, sys, time, glob, argparse, importlib
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qsolib as Q
PQS = importlib.import_module('01_passQ_sdss')
import native_hook

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
OUTDIR = '/Volumes/exdisk3/sdssdesi_telluric/phase0'
EXPMAP = os.path.expanduser('~/data/sdssdesi_expfibermap')
LW0, NPIX = 3.5502, 4447
LWAVE = LW0 + 1e-4 * np.arange(NPIX); WAVE = 10.0**LWAVE
TELLURIC = [(6855, 6960), (7150, 7350), (7580, 7720), (8110, 8370), (8920, 9800)]      # observed A
BROAD = [1033.8, 1215.67, 1240.8, 1397.0, 1549.06, 1640.4, 1908.73, 2798.75, 3426.9, 3727.1, 3869.9, 4102.9, 4341.7,
         4862.7, 4960.3, 5008.2, 6564.6]                                                 # rest A, masked +-5000 km/s
VBROAD = 5000.0; LYA = 1215.67; SNR_MIN = 20.0; CLIP = 5.0
LQ = np.log10(Q.QWAVE)


def load_basis(path):
    f = fitsio.FITS(path); E = f['EIGVEC'].read().astype(np.float64); pm = f['PIXMASK'].read().astype(bool); f.close()
    return E * pm[None, :], pm


def rest_weight_mask(z):
    """[n, QNREST] True where a rest pixel may enter the fit: outside the forest, the broad lines and the observed telluric windows."""
    n = len(z); ok = np.ones((n, len(Q.QWAVE)), bool)
    ok[:, Q.QWAVE < LYA] = False
    v = np.log10(1 + VBROAD / 2.998e5)
    for l in BROAD:
        ok[:, np.abs(LQ - np.log10(l)) < v] = False
    for i in range(n):
        lo = LQ + np.log10(1 + z[i])
        for a, b in TELLURIC:
            ok[i, (lo >= np.log10(a)) & (lo <= np.log10(b))] = False
    return ok


def fit_basis(E, rflux, rivar, ok):
    """Weighted least squares per row, one clipping pass. Returns model_rest [n, QNREST] (NaN outside the basis)."""
    n = len(rflux); model = np.full(rflux.shape, np.nan); pm = (E != 0).any(axis=0)
    for i in range(n):
        w = np.where(ok[i] & pm & (rivar[i] > 0), rivar[i], 0.0).astype(np.float64); f = rflux[i].astype(np.float64)
        for it in range(2):
            A = (E * w) @ E.T; b = (E * w) @ f
            try:
                c = np.linalg.solve(A + np.eye(len(E)) * 1e-12 * np.trace(A), b)
            except np.linalg.LinAlgError:
                c = None; break
            m = c @ E
            if it == 0:
                chi = (f - m) * np.sqrt(w); w[np.abs(chi) > CLIP] = 0.0
        if c is not None and (w > 0).sum() >= 500:
            model[i] = np.where(pm, m, np.nan)
    return model


def observed_residual(model_rest, z, w_obs, f_obs, iv_obs, good, mask_obs):
    """R = f/model - 1 and W = ivar model^2 on the object's observed pixels; model from the rest fit at lambda_obs/(1+z)."""
    n = len(z); R = np.full(f_obs.shape, np.nan, np.float32); W = np.zeros(f_obs.shape, np.float32); fscale = np.full(n, np.nan)
    lo = np.log10(w_obs); band = (w_obs >= 6000) & (w_obs <= 9000)
    for i in range(n):
        if not np.isfinite(model_rest[i]).any():
            continue
        lr = lo - np.log10(1 + z[i]); mr = model_rest[i]; fin = np.isfinite(mr)
        m = np.interp(lr, LQ[fin], mr[fin], left=np.nan, right=np.nan) / (1 + z[i])     # rest model is F_lam (1+z)
        # a basis gap (PIXMASK hole) must not be bridged: require both neighbours finite
        j = np.clip(np.searchsorted(LQ, lr), 1, len(LQ) - 1); m = np.where(fin[j] & fin[j - 1], m, np.nan)
        use = good[i] & np.isfinite(m) & (m > 0)
        with np.errstate(all='ignore'):
            R[i] = np.where(use, f_obs[i] / m - 1, np.nan)
        W[i] = np.where(use & mask_obs[i], iv_obs[i] * m * m, 0.0)
        if (use & band).any():
            fscale[i] = np.median(m[use & band]) / Q.passB.FLUX_UNIT          # FLUX_SCALE: model flux 6000-9000 A, 1e-17 cgs
    return R, W, fscale


def obs_mask(z, w_obs):
    """[n, npix] True where the observed pixel is outside the forest and the broad lines (W = 0 elsewhere)."""
    n = len(z); ok = np.ones((n, len(w_obs)), bool); lo = np.log10(w_obs); v = np.log10(1 + VBROAD / 2.998e5)
    for i in range(n):
        lr = lo - np.log10(1 + z[i]); ok[i, lr < np.log10(LYA)] = False
        for l in BROAD:
            ok[i, np.abs(lr - np.log10(l)) < v] = False
    return ok


def to_lattice_sdss(R, W, coeff0):
    k0 = int(round((coeff0 - LW0) / 1e-4)); assert abs((coeff0 - LW0) / 1e-4 - k0) < 1e-3, coeff0
    n, npix = R.shape; Rl = np.full((n, NPIX), np.nan, np.float32); Wl = np.zeros((n, NPIX), np.float32)
    a, b = max(0, -k0), min(npix, NPIX - k0)
    Rl[:, k0 + a:k0 + b] = R[:, a:b]; Wl[:, k0 + a:k0 + b] = W[:, a:b]
    return Rl, Wl


def to_lattice_desi(R, W, w_obs):
    """Linear interpolation of R and W from the 0.8 A pixels to the lattice; NaN / 0 where either neighbour is missing."""
    n = len(R); Rl = np.full((n, NPIX), np.nan, np.float32); Wl = np.zeros((n, NPIX), np.float32)
    j = np.searchsorted(w_obs, WAVE); inside = (j >= 1) & (j < len(w_obs)); jj = np.clip(j, 1, len(w_obs) - 1)
    t = (WAVE - w_obs[jj - 1]) / (w_obs[jj] - w_obs[jj - 1])
    for i in range(n):
        r0, r1 = R[i][jj - 1], R[i][jj]; ok = inside & np.isfinite(r0) & np.isfinite(r1)
        Rl[i] = np.where(ok, (1 - t) * r0 + t * r1, np.nan)
        w0, w1 = W[i][jj - 1], W[i][jj]; okw = ok & (w0 > 0) & (w1 > 0)
        Wl[i] = np.where(okw, (1 - t) * w0 + t * w1, 0.0)
    return Rl, Wl


def classes(d, survey):
    cols = ['SBIN', 'EBIN'] + (['PLATE', 'MJD', 'FIBERID'] if survey == 'sdss' else ['TARGETID', 'PROGRAM'])
    c = pd.read_parquet(f'{PROD}/qso_cells_v3_all.parquet', columns=cols + ['SURVEY'])
    c = c[c.SURVEY == 'desi'] if survey == 'desi' else c[c.SURVEY != 'desi']
    key = ['PLATE', 'MJD', 'FIBERID'] if survey == 'sdss' else ['TARGETID', 'PROGRAM']
    c = c.drop(columns=['SURVEY']).drop_duplicates(key)
    for k in key:
        if k != 'PROGRAM':
            c[k] = c[k].astype(np.int64); d[k] = d[k].astype(np.int64)
    d = d.merge(c, on=key, how='left')
    return d.SBIN.fillna(-1).to_numpy(np.int8), d.EBIN.fillna(-1).to_numpy(np.int8)


def desi_exposures(d):
    """Per quasar: dominant TILEID, EXPTIME-weighted MJD and the exposure list from the expfibermap tables (IN_COADD any camera)."""
    parts = []                                                 # one table per (program, healpix): EXPMAP/<program>/hpx<N>.parquet
    hp = d.FILE.str.extract(r'coadd-main-(\w+)-(\d+)\.fits')
    for (prog, h), g in d.assign(PR=hp[0], HP=hp[1]).groupby(['PR', 'HP']):
        p = f'{EXPMAP}/{prog}/hpx{h}.parquet'
        if not os.path.exists(p):
            print('no expfibermap', p, flush=True); continue
        e = pd.read_parquet(p, columns=['TARGETID', 'PROGRAM', 'TILEID', 'MJD', 'EXPID', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z'])
        e = e[e.TARGETID.isin(g.TARGETID.astype(np.int64)) & (e.IN_COADD_B | e.IN_COADD_R | e.IN_COADD_Z)]
        parts.append(e)
    e = pd.concat(parts, ignore_index=True)
    e = e.merge(d[['TARGETID', 'PROGRAM']].assign(TARGETID=lambda x: x.TARGETID.astype(np.int64), ORD=np.arange(len(d))), on=['TARGETID', 'PROGRAM'])
    e = e.sort_values(['ORD', 'EXPID'])
    n = len(d); plate = np.full(n, -1, np.int64); mjd = np.full(n, np.nan); off = np.zeros(n + 1, np.int64)
    cnt = e.groupby('ORD').size(); off[1:] = np.cumsum(np.array([cnt.get(i, 0) for i in range(n)]))
    for o, g in e.groupby('ORD'):
        tt = g.groupby('TILEID').EXPTIME.sum(); plate[o] = tt.idxmax(); mjd[o] = (g.MJD * g.EXPTIME).sum() / g.EXPTIME.sum()
    return plate, mjd, dict(EXP_OFFSET=off, EXP_EXPID=e.EXPID.to_numpy(np.int64), EXP_TILEID=e.TILEID.to_numpy(np.int64),
                            EXP_MJD=e.MJD.to_numpy(float), EXP_EXPTIME=e.EXPTIME.to_numpy(float))


def run_sdss(E, a):
    m = pd.read_parquet(f'{PROD}/qso_master_sdss.parquet', columns=['PLATE', 'MJD', 'FIBERID', 'SURVEY', 'Z', 'Z_SYS', 'CAL_A', 'EBV', 'SNR_ALL', 'QCUT'])
    m = m[m.QCUT & (m.SNR_ALL > SNR_MIN)].copy()
    if a.nmax:
        m = m.sort_values('SNR_ALL', ascending=False).head(a.nmax)
    m['ZS'] = np.where(np.isfinite(m.Z_SYS), m.Z_SYS, m.Z)
    paths = {}
    for sv in ('dr17', 'dr8'):
        for p in PQS.plate_files(sv):
            paths[os.path.basename(p).replace('.fits', '')] = p
    m['FILE'] = ['spPlate-%04d-%05d' % (p, j) for p, j in zip(m.PLATE, m.MJD)]
    m = m[m.FILE.isin(paths)].reset_index(drop=True)
    print(f'sdss: {len(m):,} quasars (QCUT, S/N > {SNR_MIN:.0f}) in {m.FILE.nunique():,} plates', flush=True)
    Rs, Ws, rows = [], [], []; t0 = time.time()
    for k, (fname, sub) in enumerate(m.groupby('FILE', sort=True)):
        F = fitsio.FITS(paths[fname]); h = F[0].read_header(); coeff0 = float(h['COEFF0'])
        fib = sub.FIBERID.to_numpy().astype(int) - 1
        flux = F[0].read()[fib].astype(float); ivar = F[1].read()[fib].astype(float); mask = F[2].read()[fib]; F.close()
        npix = flux.shape[1]; w_obs = 10.0**(coeff0 + np.arange(npix) * 1e-4); z = sub.ZS.to_numpy(float)
        rid0 = int(round(coeff0 / Q.DLOG)); sv = sub.SURVEY.iloc[0]
        rflux, rivar = Q.rest_from_loggrid(rid0, flux, ivar, mask, z, np.nan_to_num(sub.CAL_A.to_numpy(float), nan=1.0),
                                           sub.EBV.to_numpy(float), coeff0=coeff0, survey=sv)
        model = fit_basis(E, rflux, rivar, rest_weight_mask(z))
        # the observed calibrated spectrum, as rest_from_loggrid builds it
        good = (ivar > 0) & (mask == 0) & Q.RF.sdss_badpix_ok(coeff0, npix, sv)[None, :]
        corr = (np.nan_to_num(sub.CAL_A.to_numpy(float), nan=1.0) * Q.passB.FLUX_UNIT)[:, None]
        A = (3.1 * 0.86 * np.nan_to_num(sub.EBV.to_numpy(float)))[:, None] * Q.passB.extinction.fitzpatrick99(w_obs, 1.0, 3.1)[None, :]
        f_obs = flux * corr * 10**(0.4 * A); iv_obs = np.where(good, ivar / corr**2 * 10**(-0.8 * A), 0.0)
        R, W, fs = observed_residual(model, z, w_obs, f_obs, iv_obs, good, obs_mask(z, w_obs))
        Rl, Wl = to_lattice_sdss(R, W, coeff0); Rs.append(Rl); Ws.append(Wl); rows.append(sub.assign(FLUX_SCALE=fs))
        if k % 200 == 0:
            print(f'{k}/{m.FILE.nunique()} {fname} {time.time()-t0:.0f}s', flush=True)
    d = pd.concat(rows, ignore_index=True); R = np.vstack(Rs); W = np.vstack(Ws)
    sb, eb = classes(d.copy(), 'sdss')
    ids = np.array([f'{p}-{j}-{f}' for p, j, f in zip(d.PLATE, d.MJD, d.FIBERID)], dtype='U24')
    return R, W, dict(ID=ids, CLASS=sb, EWCLASS=eb, SURVEY=np.asarray(d.SURVEY.tolist(), dtype='U8'), Z=d.ZS.to_numpy(float), SNR=d.SNR_ALL.to_numpy(float),
                      PLATE=d.PLATE.to_numpy(np.int64), MJD=d.MJD.to_numpy(float), AIRMASS=np.full(len(d), np.nan), PROGRAM=np.asarray(d.SURVEY.tolist(), dtype='U8'),
                      FLUX_SCALE=d.FLUX_SCALE.to_numpy(float), EBV=d.EBV.to_numpy(float))


def run_desi(E, a):
    cols = ['TARGETID', 'PROGRAM', 'FILE', 'Z', 'Z_SYS', 'CAL_A', 'CAL_B', 'EBV', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'SNR_SDSSPIX', 'QCUT']
    m = pd.read_parquet(f'{PROD}/qso_master_desi.parquet', columns=cols)
    m = m[m.QCUT & (m.SNR_SDSSPIX > SNR_MIN)].copy(); m['ZS'] = np.where(np.isfinite(m.Z_SYS), m.Z_SYS, m.Z)
    cnt = m.groupby('FILE').size().sort_values(ascending=False)              # files with the most high-S/N quasars first
    if a.nmax:
        keep = cnt.index[np.cumsum(cnt.to_numpy()) <= a.nmax]; m = m[m.FILE.isin(keep)]
    m = m.reset_index(drop=True)
    print(f'desi: {len(m):,} quasars (QCUT, S/N > {SNR_MIN:.0f}) in {m.FILE.nunique():,} coadds ({(cnt[cnt.index.isin(m.FILE.unique())] >= 5).sum()} with >= 5)', flush=True)
    Rs, Ws, rows = [], [], []; t0 = time.time(); nf = m.FILE.nunique()
    for k, (fname, sub) in enumerate(m.groupby('FILE', sort=True)):
        npz = native_hook.out(fname)
        if not os.path.exists(npz):
            print('NOT STORED', fname, flush=True); continue
        zf = np.load(npz); tid = zf['TARGETID'].astype(np.int64); flux = zf['FLUX']; ivar = zf['IVAR']
        w_obs = float(zf['W0']) + float(zf['DW']) * np.arange(flux.shape[1]); del zf
        pos = pd.Series(np.arange(len(tid)), index=tid); sub = sub[sub.TARGETID.astype(np.int64).isin(pos.index)]
        if not len(sub):
            continue
        r = pos.loc[sub.TARGETID.astype(np.int64).to_numpy()].to_numpy(); flux = flux[r].astype(float); ivar = ivar[r].astype(float)
        mask = (ivar <= 0).astype(np.int8); df = sub.copy(); df['Z'] = df['ZS']; z = df.ZS.to_numpy(float)
        rflux, rivar = Q.quasar_rest(df, w_obs, flux, ivar, mask, exact=True)
        model = fit_basis(E, rflux, rivar, rest_weight_mask(z))
        # the observed calibrated spectrum, as quasar_rest builds it
        psf = df.MEAN_PSF_TO_FIBER_SPECFLUX.to_numpy(float); psf = np.where(np.isfinite(psf) & (psf > 0), psf, 1.0)
        ca = np.nan_to_num(df.CAL_A.to_numpy(float), nan=1.0); cb = np.nan_to_num(df.CAL_B.to_numpy(float))
        corr = ca[:, None] + cb[:, None] * np.log10(w_obs / 6000.0)[None, :]; corr = np.where(corr > 0, corr, 1.0) * (psf * Q.passB.FLUX_UNIT)[:, None]
        A = (3.1 * 0.86 * np.nan_to_num(df.EBV.to_numpy(float)))[:, None] * Q.passB.extinction.fitzpatrick99(w_obs, 1.0, 3.1)[None, :]
        good = ivar > 0; f_obs = flux * corr * 10**(0.4 * A); iv_obs = np.where(good, ivar / corr**2 * 10**(-0.8 * A), 0.0)
        R, W, fs = observed_residual(model, z, w_obs, f_obs, iv_obs, good, obs_mask(z, w_obs))
        Rl, Wl = to_lattice_desi(R, W, w_obs); Rs.append(Rl); Ws.append(Wl); rows.append(sub.assign(FLUX_SCALE=fs))
        if k % 200 == 0:
            print(f'{k}/{nf} {fname} {time.time()-t0:.0f}s', flush=True)
    d = pd.concat(rows, ignore_index=True); R = np.vstack(Rs); W = np.vstack(Ws)
    sb, eb = classes(d.copy(), 'desi')
    plate, mjd, exp = desi_exposures(d)
    return R, W, dict(ID=np.asarray(d.TARGETID.astype(np.int64).astype(str).tolist(), dtype='U24'), CLASS=sb, EWCLASS=eb, SURVEY=np.array(['desi'] * len(d)),
                      Z=d.ZS.to_numpy(float), SNR=d.SNR_SDSSPIX.to_numpy(float), PLATE=plate, MJD=mjd, AIRMASS=np.full(len(d), np.nan),
                      PROGRAM=np.asarray(d.PROGRAM.tolist(), dtype='U8'),
                      FLUX_SCALE=d.FLUX_SCALE.to_numpy(float), EBV=d.EBV.to_numpy(float), **exp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', required=True, choices=['sdss', 'desi'])
    ap.add_argument('--basis', default=f'{PROD}/qso_empca_v3.fits')
    ap.add_argument('--nmax', type=int, default=0, help='desi: about this many quasars, files with the most S/N > 20 quasars first')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(); out = a.out or f'{OUTDIR}/quasar_{a.survey}.npz'
    E, pm = load_basis(a.basis); print(f'basis {a.basis}: {len(E)} vectors, {pm.sum()} pixels', flush=True)
    R, W, meta = (run_sdss if a.survey == 'sdss' else run_desi)(E, a)
    os.makedirs(os.path.dirname(out), exist_ok=True); tmp = out + '.tmp.npz'
    np.savez(tmp, WAVE=WAVE.astype(np.float64), R=R.astype(np.float32), W=W.astype(np.float32), NVEC=len(E) - 1,
             TELLURIC_WINDOWS=np.array(TELLURIC, float), BASIS=os.path.basename(a.basis), **meta)
    os.replace(tmp, out)
    cov = np.isfinite(R).mean(axis=0)
    print(f'wrote {out}: {len(R):,} quasars; residual rms (W>0) {np.sqrt(np.nanmean(np.where(W > 0, R, np.nan)**2)):.4f}; '
          f'pixel coverage 50 % over {WAVE[cov >= 0.5].min():.0f}-{WAVE[cov >= 0.5].max():.0f} A; plates/tiles {len(np.unique(meta["PLATE"]))}', flush=True)


if __name__ == '__main__':
    main()
