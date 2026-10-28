#!/usr/bin/env python3
"""Step 2b (Nao, 2026-10-24 night): joint spectrophotometric fit on the pass tables (catalog level, no spectra read).

Model, for star i and band b (catalog minus synthetic magnitude of the raw pipeline spectrum, r_ib):
    r_ib = ZP_b + k_b (c_i - 1) + a_i + t_i x_b
  a_i   gray scale of the star's spectrum (mag; the spectrum is multiplied by 10^(-0.4 a_i))
  t_i   tilt (mag per dex), x_b = log10(lambda_eff,b / 5500 A)
  ZP_b, k_b  zero point and color term ((BP-RP) - 1) of band b in its photometric system, solved from all stars
Reference system: Gaia DR3 G, BP, RP have ZP = k = 0 (Gaia G sets the scale, Gaia colors the absolute tilt); the SDSS PSF
u g r i z (SDSS stars) and the Legacy Surveys g r z, separately for the DECam south and BASS/MzLS north systems (DESI stars),
constrain the per-star tilt and get their zero points from the ensemble. Bands with more than 30 % of their weight in the
padded (extrapolated) part of the spectrum are dropped (SDSS u). Pan-STARRS can be added as another system when on disk.
Solver: alternating weighted least squares (per star: a_i, t_i from the normal equations; per band: ZP_b, k_b by weighted
linear regression on the star residuals), 3-sigma clipping per band, 6 iterations. Weights 1/(sigma_syn^2 + floor_b^2),
floors 0.01 (Gaia) and 0.02 (ground). Stars with fewer than 3 bands or a wavelength lever < 0.1 dex get a scale only.
Output: ~/data/sdssdesi_cache/stars/spcal_<survey>.parquet (SID, A_MAG, A_ERR, A_FLUX = 10^(-0.4 A_MAG), TILT, TILT_ERR,
TILT_OK (TILT_ERR < --tilt-maxerr), NBAND, CHI2NU) and zp_<survey>.csv (band, system, ZP, k, scatter, N). step6_stack.py applies
A_FLUX and, where TILT_OK, the tilt.
Usage: python3 step2b_spcal.py --survey sdss|desi [--passdir DIR] [--tag X] [--tilt-maxerr 0.1] [--no-colorterm]
"""
import os, sys, glob, argparse, time
import numpy as np, pandas as pd
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
passB = starlib.passB
CACHE = os.path.expanduser('~/data/sdssdesi_cache/stars'); os.makedirs(CACHE, exist_ok=True)
PASSDIR = {'sdss': os.path.expanduser('~/data/sdssdesi_passS_sdss'), 'desi': os.path.expanduser('~/data/sdssdesi_passS')}
FLOOR = {'gaia': 0.01, 'ground': 0.02}
EXTRAP_MAX = 0.30
ANCHOR = ('gaia:G', 'gaia:BP', 'gaia:RP')


def leff(name):
    return float(passB.filt(name).effective_wavelength.value)


def long_table(df, survey):
    """one row per (star, band): SID, band key 'system:band', x, r, sigma, color."""
    rows = []
    c = (df.GAIA_PHOT_BP_MEAN_MAG - df.GAIA_PHOT_RP_MEAN_MAG).to_numpy(float)
    def add(key, mcat, msyn, msig, x, extrap=0.0, floor=FLOOR['ground']):
        if extrap > EXTRAP_MAX:
            return
        r = np.asarray(mcat, float) - np.asarray(msyn, float)
        s = np.sqrt(np.nan_to_num(np.broadcast_to(np.asarray(msig, float), r.shape), nan=0.0)**2 + floor**2)
        ok = np.isfinite(r) & np.isfinite(c)
        rows.append(pd.DataFrame({'SID': df.SID[ok].to_numpy(), 'KEY': key, 'X': x, 'R': r[ok], 'SIG': s[ok], 'C': c[ok]}))
    for b, name in starlib.GAIA_FILTERS.items():
        add(f'gaia:{b}', df[f'GAIA_PHOT_{b}_MEAN_MAG'].where((df[f'GAIA_PHOT_{b}_MEAN_MAG'] > 0) & (df[f'GAIA_PHOT_{b}_MEAN_MAG'] < 30)),
            df[f'SYN_GAIA_{b}'], 0.0, np.log10(leff(name) / 5500), floor=FLOOR['gaia'])
    if survey == 'sdss':
        for b, name in starlib.SDSS_FILTERS.items():
            mcat = df[f'PSFMAG_{b}'].where((df[f'PSFMAG_{b}'] > 10) & (df[f'PSFMAG_{b}'] < 30))
            add(f'sdss:{b}', mcat, df[f'SYN_SDSS_{b}'], df[f'SYN_SDSS_{b}_ERR'], np.log10(leff(name) / 5500), extrap=float(df[f'SYN_SDSS_EXTRAP_{b}'].iloc[0]))
    else:
        psf2fib = df['MEAN_PSF_TO_FIBER_SPECFLUX'].to_numpy(float)
        for south, sysname in ((True, 'ls-s'), (False, 'ls-n')):
            m = (df.PHOTSYS.to_numpy() == 'S') == south
            if not m.any():
                continue
            sub = df[m]
            for b, name in starlib.LS_FILTERS[south].items():
                with np.errstate(all='ignore'):
                    mcat = 22.5 - 2.5 * np.log10(sub[f'FLUX_{b}'].to_numpy(float))
                    msyn = sub[f'SYN_FIBER_{b}_RAW'].to_numpy(float) + 2.5 * np.log10(psf2fib[m])     # point-source (total) synthetic mag
                err = sub[f'SYN_FIBER_{b}_ERR'].to_numpy(float) if f'SYN_FIBER_{b}_ERR' in sub else np.full(m.sum(), 0.02)
                r = mcat - msyn; ok = np.isfinite(r) & np.isfinite(c[m])
                rows.append(pd.DataFrame({'SID': sub.SID[ok].to_numpy(), 'KEY': f'{sysname}:{b}', 'X': np.log10(leff(name) / 5500),
                                          'R': r[ok], 'SIG': np.sqrt(np.nan_to_num(err[ok], nan=0.0)**2 + FLOOR['ground']**2), 'C': c[m][ok]}))
    return pd.concat(rows, ignore_index=True)


def solve(L, niter=6, colorterm=True, clip=3.0):
    sid, sidx = np.unique(L.SID.to_numpy(), return_inverse=True)
    key, kidx = np.unique(L.KEY.to_numpy(), return_inverse=True)
    anchor = np.array([k in ANCHOR for k in key])
    x = L.X.to_numpy(float); r = L.R.to_numpy(float); w0 = 1.0 / L.SIG.to_numpy(float)**2; cc = L.C.to_numpy(float) - 1.0
    ns, nk = len(sid), len(key)
    zp = np.zeros(nk); kc = np.zeros(nk); w = w0.copy()
    # per-star lever arm: tilt only where the bands span >= 0.1 dex and n >= 3
    xmin = np.full(ns, np.inf); xmax = np.full(ns, -np.inf)
    np.minimum.at(xmin, sidx, x); np.maximum.at(xmax, sidx, x)
    nb = np.bincount(sidx, minlength=ns)
    tilt_ok = (xmax - xmin >= 0.1) & (nb >= 3)
    a = np.zeros(ns); t = np.zeros(ns)
    for it in range(niter):
        y = r - zp[kidx] - kc[kidx] * cc
        S0 = np.bincount(sidx, w, ns); Sx = np.bincount(sidx, w * x, ns); Sxx = np.bincount(sidx, w * x * x, ns)
        Sy = np.bincount(sidx, w * y, ns); Sxy = np.bincount(sidx, w * x * y, ns)
        det = S0 * Sxx - Sx**2
        with np.errstate(all='ignore'):
            t_new = np.where(tilt_ok & (det > 0), (S0 * Sxy - Sx * Sy) / det, 0.0)
            a_new = np.where(tilt_ok & (det > 0), (Sxx * Sy - Sx * Sxy) / det, np.where(S0 > 0, Sy / S0, np.nan))
        a, t = a_new, t_new
        e = y - a[sidx] - t[sidx] * x
        # per band: zero point and color term from the stars (anchors fixed), robust scatter, clipping
        for k in range(nk):
            m = (kidx == k) & (w > 0) & np.isfinite(e)
            if not anchor[k] and m.sum() > 50:
                ww = w[m]; ek = e[m] + zp[k] + kc[k] * cc[m]           # residual before this band's own terms
                if colorterm:
                    A = np.stack([np.ones(m.sum()), cc[m]], axis=1); N = (A * ww[:, None]).T @ A; bvec = (A * ww[:, None]).T @ ek
                    zp[k], kc[k] = np.linalg.solve(N, bvec)
                else:
                    zp[k] = np.sum(ww * ek) / np.sum(ww)
        y = r - zp[kidx] - kc[kidx] * cc; e = y - a[sidx] - t[sidx] * x
        sig_b = np.array([1.4826 * np.median(np.abs(e[(kidx == k) & np.isfinite(e)])) if ((kidx == k) & np.isfinite(e)).any() else np.nan for k in range(nk)])
        thresh = clip * np.maximum(1.0 / np.sqrt(w0), sig_b[kidx])
        w = np.where(np.abs(e) > thresh, 0.0, w0)
    # errors from the normal equations, scaled by sqrt(chi2nu) when > 1
    S0 = np.bincount(sidx, w, ns); Sx = np.bincount(sidx, w * x, ns); Sxx = np.bincount(sidx, w * x * x, ns); det = S0 * Sxx - Sx**2
    chi2 = np.bincount(sidx, w * e**2, ns); nused = np.bincount(sidx, (w > 0).astype(float), ns)
    dof = np.where(tilt_ok, nused - 2, nused - 1); chi2nu = np.where(dof > 0, chi2 / np.where(dof > 0, dof, 1), np.nan)
    infl = np.sqrt(np.where(np.isfinite(chi2nu) & (chi2nu > 1), chi2nu, 1.0))
    with np.errstate(all='ignore'):
        a_err = np.where(tilt_ok & (det > 0), np.sqrt(Sxx / det), np.sqrt(1.0 / S0)) * infl
        t_err = np.where(tilt_ok & (det > 0), np.sqrt(S0 / det), np.nan) * infl
    stars = pd.DataFrame({'SID': sid, 'A_MAG': a, 'A_ERR': a_err, 'TILT': t, 'TILT_ERR': t_err, 'NBAND': nused.astype(int), 'CHI2NU': chi2nu})
    bands = pd.DataFrame({'KEY': key, 'SYSTEM': [k.split(':')[0] for k in key], 'BAND': [k.split(':')[1] for k in key], 'ANCHOR': anchor,
                          'ZP': zp, 'COLOR_TERM': kc, 'SCATTER': sig_b, 'N': [int(((kidx == k) & (w > 0)).sum()) for k in range(nk)],
                          'N_CLIPPED': [int(((kidx == k) & (w == 0)).sum()) for k in range(nk)], 'X': [x[kidx == k][0] for k in range(nk)]})
    return stars, bands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='sdss', choices=['sdss', 'desi']); ap.add_argument('--passdir'); ap.add_argument('--tag')
    ap.add_argument('--tilt-maxerr', type=float, default=0.1); ap.add_argument('--no-colorterm', action='store_true'); ap.add_argument('--snmin', type=float, default=3.0)
    a = ap.parse_args(); t0 = time.time()
    files = sorted(glob.glob(os.path.join(a.passdir or PASSDIR[a.survey], '*', '*.parquet')))
    files = [f for f in files if 'FILE' in pq.ParquetFile(f).schema_arrow.names]
    need = ['ZWARN', 'SNR_ALL', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'SYN_GAIA_G', 'SYN_GAIA_BP', 'SYN_GAIA_RP', 'GAIA_DG', 'CAL_B_EST']
    need += (['PLATE', 'MJD', 'FIBERID'] + [f'{p}_{b}' for b in 'UGRIZ' for p in ('PSFMAG', 'SYN_SDSS', 'SYN_SDSS_EXTRAP')] + [f'SYN_SDSS_{b}_ERR' for b in 'UGRIZ']) if a.survey == 'sdss' \
            else ['TARGETID', 'PROGRAM', 'PHOTSYS', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'SYN_FIBER_G_RAW', 'SYN_FIBER_R_RAW', 'SYN_FIBER_Z_RAW', 'SYN_FIBER_G_ERR', 'SYN_FIBER_R_ERR']
    df = pd.concat([pd.read_parquet(f, columns=[c for c in need if c in pq.ParquetFile(f).schema_arrow.names]) for f in files], ignore_index=True)
    df['SID'] = (df.PLATE.astype(str) + '-' + df.MJD.astype(str) + '-' + df.FIBERID.astype(str)) if a.survey == 'sdss' else starlib.desi_sid(df)
    df = df[(df.ZWARN == 0) & (df.SNR_ALL > a.snmin)].reset_index(drop=True)
    L = long_table(df, a.survey)
    print(f'{a.survey}: {len(files)} files, {len(df):,} stars (ZWARN=0, S/N>{a.snmin}), {len(L):,} star-band pairs, bands {sorted(L.KEY.unique())}; {time.time()-t0:.0f}s', flush=True)
    stars, bands = solve(L, colorterm=not a.no_colorterm)
    stars['A_FLUX'] = 10**(-0.4 * stars.A_MAG); stars['TILT_OK'] = np.isfinite(stars.TILT_ERR) & (stars.TILT_ERR < a.tilt_maxerr)
    stars = stars.merge(df[['SID', 'SNR_ALL', 'GAIA_DG', 'CAL_B_EST']], on='SID', how='left')
    tag = a.tag or a.survey
    stars.to_parquet(os.path.join(CACHE, f'spcal_{tag}.parquet')); bands.to_csv(os.path.join(CACHE, f'zp_{tag}.csv'), index=False)
    print('\n== zero points relative to Gaia (ZP = catalog - synthetic after the per-star scale and tilt; color term per mag of (BP-RP) - 1)')
    print(bands.sort_values('X')[['KEY', 'X', 'ZP', 'COLOR_TERM', 'SCATTER', 'N', 'N_CLIPPED', 'ANCHOR']].to_string(index=False, float_format=lambda v: f'{v:+.4f}'))
    def q(v):
        v = np.asarray(v, float); v = v[np.isfinite(v)]
        return (np.median(v), np.diff(np.percentile(v, [16, 84]))[0] / 2) if len(v) else (np.nan, np.nan)
    print('\n== per-star solution')
    print(f'  scale A_MAG: median {q(stars.A_MAG)[0]:+.3f} sigma68 {q(stars.A_MAG)[1]:.3f}   (Gaia G alone, GAIA_DG: {q(stars.GAIA_DG)[0]:+.3f} / {q(stars.GAIA_DG)[1]:.3f})   median A_ERR {np.nanmedian(stars.A_ERR):.3f}')
    for lo, hi in ((3, 10), (10, 30), (30, 1e9)):
        s = stars[(stars.SNR_ALL >= lo) & (stars.SNR_ALL < hi)]
        print(f'  S/N {lo:>3}-{hi if hi < 1e9 else "":<4} n={len(s):>7,}: tilt median {q(s.TILT)[0]:+.3f} sigma68 {q(s.TILT)[1]:.3f} mag/dex, median TILT_ERR {np.nanmedian(s.TILT_ERR):.3f}, TILT_OK {s.TILT_OK.mean()*100:4.1f} %, chi2nu median {np.nanmedian(s.CHI2NU):.2f}')
    m = np.isfinite(stars.TILT) & np.isfinite(stars.CAL_B_EST) & stars.TILT_OK
    if m.sum() > 30:
        print(f'  tilt vs pass-B two-band tilt estimate CAL_B_EST (flux ratio per dex): corr {np.corrcoef(stars.TILT[m], stars.CAL_B_EST[m])[0,1]:+.2f} (n={m.sum()})')
    print(f'  tilt scatter beyond the errors (sqrt(sigma68^2 - median err^2)): {np.sqrt(max(q(stars.TILT[stars.TILT_OK])[1]**2 - np.nanmedian(stars.TILT_ERR[stars.TILT_OK])**2, 0)):.3f} mag/dex  -> the real per-spectrum tilt to solve per exposure/petal/fiber')
    print(f'written {CACHE}/spcal_{tag}.parquet, zp_{tag}.csv  {time.time()-t0:.0f}s')
    try:
        import json; S = os.path.expanduser('~/data/sdssdesi_status')
        json.dump({'done': 1 if a.survey == 'sdss' else 2, 'total': 2, 'status': 'waiting' if a.survey == 'sdss' else 'done', 'note': f'{a.survey}: {len(stars):,} stars, tilt scatter beyond errors {np.sqrt(max(q(stars.TILT[stars.TILT_OK])[1]**2 - np.nanmedian(stars.TILT_ERR[stars.TILT_OK])**2, 0)):.2f} mag/dex' + ('; DESI after pass B' if a.survey == 'sdss' else '')}, open(f'{S}/step_star_spcal.json', 'w'))
    except Exception:
        pass


if __name__ == '__main__':
    main()
