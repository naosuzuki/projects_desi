#!/usr/bin/env python3
"""Telluric Phase 0, galaxies (docs/PLAN_telluric.md; Nao 2026-11-01): observed-frame residuals of high-S/N individual
galaxies for the joint telluric basis, in the common format of programs_stars/telluric_phase0.py. No new DESI / SDSS read:
the spectra come from the v2 rest-frame stores on exdisk3 (float16 normalized flux + ivar on the 1e-4 log grid), the
selection (S/N > 20, SNR_SDSSPIX) and the type (DB_EFF) from the v3 per-galaxy tables. Whole healpix files / plates are taken,
so each tile or plate has many galaxies for the coherence test.
Per galaxy: weighted LSQ fit with the galaxy basis (--basis, default the 28 Oct release basis, rest 1400-9900 A) on the basis
pixels OUTSIDE the observed telluric windows and the rest-frame emission lines (+-1000 km/s), so the fit is not pulled by
telluric residuals; R = f / model - 1 and W = ivar x model^2 (0 where masked, NaN / 0 outside the data), moved to the observed
frame on the stellar grid SWAVE (log10 lambda = 3.5502 + 1e-4 i, i < 4447) by linear interpolation of the log grid.
Metadata: ID (DESI TARGETID-program, SDSS plate-mjd-fiber), CLASS DB<DB_EFF>, Z, SNR, PLATE (SDSS plate; DESI the TILEID with
the most EXPTIME in the coadd), MJD (SDSS plate MJD; DESI EXPTIME-weighted mean), AIRMASS NaN; and for DESI the exposure rows
EXP_ID / EXP_EXPID / EXP_TILEID / EXP_MJD / EXP_EXPTIME (flat arrays, from ~/data/sdssdesi_expfibermap, IN_COADD any camera).
FLUX_SCALE (fibre r photometry as F_lambda at 6400 A, 1e-17 erg/s/cm2/A) and EBV_SFD are added for the sky / telluric
separation. MAPPED_FROM = 'rest_store_v2' records that R was mapped back from the sub-pixel rest-frame store (one extra resample).
Output: /Volumes/exdisk3/sdssdesi_telluric/phase0/galaxy_<desi|sdss>.npz.
Usage: python3 26_telluric_phase0.py --survey desi|sdss [--target 20000] [--basis FITS]
"""
import os, sys, glob, json, time, argparse
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')

OUT = '/Volumes/exdisk3/sdssdesi_telluric/phase0'
STORE = {'desi': '/Volumes/exdisk3/sdssdesi_restframe_v2', 'sdss': '/Volumes/exdisk3/sdssdesi_restframe_sdss_v2'}
TABLES = {'desi': os.path.expanduser('~/data/sdssdesi_passC_v3/tables'), 'sdss': os.path.expanduser('~/data/sdssdesi_passG_sdss_v3/tables')}
EXPF = os.path.expanduser('~/data/sdssdesi_expfibermap')
BASIS = os.path.join(HERE, '..', 'data_eigenspectra', 'galaxy_empca_desisdss_v3.fits')
TELLURIC = [(6855, 6960), (7150, 7350), (7580, 7720), (8110, 8370), (8920, 9800)]     # observed A (the plan's windows)
LINES = [3727.09, 3729.88, 4862.68, 4960.30, 5008.24, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67]   # rest vacuum A
VLINE = 1000.0 / 299792.458
SWAVE = 10 ** (3.5502 + 1e-4 * np.arange(4447)); LWS = np.log10(SWAVE)
SNR_MIN = 20.0
STATUS = os.path.expanduser('~/data/sdssdesi_status/step_gal_telluric_phase0.json')


def select(survey, target, stride):
    """S/N > 20 galaxies from every stride-th v3 table file (whole files) until target is reached."""
    cols = ['SNR_SDSSPIX', 'Z', 'DB_EFF', 'EBV'] + (['TARGETID', 'FIBERFLUX_R'] if survey == 'desi' else ['PLATE', 'MJD', 'FIBERID', 'FIB_r'])
    files = sorted(glob.glob(os.path.join(TABLES[survey], '*', '*.parquet')))[::stride]
    parts = []; n = 0
    for f in files:
        try:
            t = pd.read_parquet(f, columns=cols)
        except Exception:
            continue
        t = t[(t.SNR_SDSSPIX > SNR_MIN) & np.isfinite(t.Z) & (t.Z > 0)]
        if not len(t):
            continue
        sub = os.path.basename(os.path.dirname(f)); t['SUB'] = sub
        t['FILE'] = os.path.basename(f).replace('.parquet', '')
        parts.append(t); n += len(t)
        if n >= target:
            break
    return pd.concat(parts, ignore_index=True)


def read_store(survey, sub, fname, keys):
    """store spectra of the wanted keys (TARGETID / FIBERID) in one file: dict key -> (rid0, flux, ivar)."""
    p = os.path.join(STORE[survey], sub, (fname if survey == 'desi' else fname) + '.npz')
    if not os.path.exists(p):
        return {}
    z = np.load(p); ids = z['TARGETID'] if survey == 'desi' else z['FIBERID']
    st, off, fl, iv, rid0 = z['START'], z['OFFSET'], z['FLUX'], z['IVAR'], int(z['RID_MIN'])
    want = set(keys); out = {}
    for i, k in enumerate(ids):
        if int(k) in want:
            out[int(k)] = (rid0 + int(st[i]), fl[off[i]:off[i + 1]].astype(np.float32), iv[off[i]:off[i + 1]].astype(np.float32))
    return out


def fit_one(rid0, fl, iv, z, B, pm, bw, lw):
    """fit on the basis grid, residual and weight moved to the observed SWAVE grid; None if too few pixels."""
    f = np.zeros(len(bw), np.float32); v = np.zeros(len(bw), np.float32)
    i0 = rid0 - pb.RID_MIN; i1 = i0 + len(fl)
    j0, j1 = max(i0, 0), min(i1, len(bw))
    if j1 <= j0:
        return None
    f[j0:j1] = fl[j0 - i0:j1 - i0]; v[j0:j1] = iv[j0 - i0:j1 - i0]
    wo = bw * (1 + z); tel = np.zeros(len(bw), bool)
    for lo, hi in TELLURIC:
        tel |= (wo >= lo) & (wo <= hi)
    em = np.zeros(len(bw), bool)
    for l in LINES:
        em |= (bw >= l * (1 - VLINE)) & (bw <= l * (1 + VLINE))
    ok = pm & (v > 0) & np.isfinite(f); fit = ok & ~tel & ~em
    if fit.sum() < 500:
        return None
    sw = np.sqrt(v[fit]); coef = np.linalg.lstsq((B[:, fit] * sw).T, f[fit] * sw, rcond=None)[0]
    model = coef @ B
    good = ok & (model > 0) & ~em
    r = np.where(good, f / np.where(good, model, 1) - 1, np.nan); wr = np.where(good, v * model**2, 0.0)
    sh = np.log10(1 + z)
    R = np.interp(LWS, lw + sh, r, left=np.nan, right=np.nan)
    Wi = np.interp(LWS, lw + sh, wr, left=0.0, right=0.0)
    bad = np.interp(LWS, lw + sh, (~good).astype(float), left=1.0, right=1.0) > 0
    R[bad] = np.nan; Wi[bad] = 0.0
    return R.astype(np.float32), Wi.astype(np.float32)


def desi_exposures(sel):
    """per galaxy: dominant TILEID, EXPTIME-weighted MJD; and the flat exposure rows."""
    parts = []
    for (sub, fname), g in sel.groupby(['SUB', 'FILE']):
        p = os.path.join(EXPF, sub, fname + '.parquet')
        if not os.path.exists(p):
            continue
        e = pd.read_parquet(p, columns=['TARGETID', 'EXPID', 'TILEID', 'MJD', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z'])
        e = e[e.TARGETID.isin(g.TARGETID) & (e.IN_COADD_B | e.IN_COADD_R | e.IN_COADD_Z) & (e.EXPTIME > 0)].copy()
        e['ID'] = e.TARGETID.astype(str) + '-' + sub; parts.append(e)
    e = pd.concat(parts, ignore_index=True)
    e['WM'] = e.MJD * e.EXPTIME
    mjd = e.groupby('ID').agg(WM=('WM', 'sum'), EXP=('EXPTIME', 'sum')); mjd['MJD'] = mjd.WM / mjd.EXP
    tile = e.groupby(['ID', 'TILEID']).EXPTIME.sum().reset_index().sort_values('EXPTIME').drop_duplicates('ID', keep='last').set_index('ID')
    return tile.TILEID.reindex(sel.ID).fillna(-1).to_numpy(np.int64), mjd.MJD.reindex(sel.ID).to_numpy(float), e


def status(done, total, state, note):
    json.dump(dict(done=done, total=total, status=state, started=T0, note=note), open(STATUS, 'w'))


def main():
    global T0
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', required=True, choices=['desi', 'sdss'])
    ap.add_argument('--target', type=int, default=20000); ap.add_argument('--stride', type=int, default=0)
    ap.add_argument('--basis', default=BASIS); ap.add_argument('--nvec', type=int, default=11)
    a = ap.parse_args(); T0 = time.time()
    stride = a.stride or (26 if a.survey == 'desi' else 4)
    F = fitsio.FITS(a.basis); B = F['EIGVEC'].read()[:a.nvec].astype(float); pm = F['PIXMASK'].read().astype(bool); bw = F['WAVE'].read(); F.close()
    lw = np.log10(bw)
    if len(bw) != pb.NREST:
        raise SystemExit(f'basis grid {len(bw)} != pass-C grid {pb.NREST}')
    sel = select(a.survey, a.target, stride)
    if a.survey == 'desi':
        sel['ID'] = sel.TARGETID.astype(str) + '-' + sel.SUB; sel['KEY'] = sel.TARGETID.astype(np.int64)
    else:
        sel['ID'] = sel.PLATE.astype(str) + '-' + sel.MJD.astype(str) + '-' + sel.FIBERID.astype(str); sel['KEY'] = sel.FIBERID.astype(np.int64)
    sel = sel.drop_duplicates('ID').reset_index(drop=True); n = len(sel)
    groups = list(sel.groupby(['SUB', 'FILE']))
    print(f'{a.survey}: {n} galaxies with S/N > {SNR_MIN:.0f} in {len(groups)} files (stride {stride}); basis {os.path.basename(a.basis)} {a.nvec} vectors', flush=True)
    status(0, len(groups), 'running', f'{a.survey}: {n} galaxies in {len(groups)} files')
    R = np.full((n, len(SWAVE)), np.nan, np.float32); W = np.zeros((n, len(SWAVE)), np.float32); fitted = np.zeros(n, bool)
    for k, ((sub, fname), g) in enumerate(groups):
        spec = read_store(a.survey, sub, fname, g.KEY.tolist())
        for i, key, z in zip(g.index, g.KEY, g.Z):
            if key not in spec:
                continue
            res = fit_one(*spec[key], float(z), B, pm, bw, lw)
            if res is None:
                continue
            R[i], W[i] = res; fitted[i] = True
        if k % 50 == 0 or k == len(groups) - 1:
            status(k + 1, len(groups), 'running', f'{a.survey}: {int(fitted.sum())} of {n} galaxies fitted, {time.time() - T0:.0f}s')
    keep = fitted; sel = sel[keep].reset_index(drop=True); R, W = R[keep], W[keep]; n = len(sel)
    cls = ('DB' + sel.DB_EFF.astype(int).astype(str)).to_numpy(str)
    extra = {}
    if a.survey == 'sdss':
        plate = sel.PLATE.to_numpy(np.int64); mjd = sel.MJD.to_numpy(float)
    else:
        plate, mjd, e = desi_exposures(sel)
        extra = dict(EXP_ID=e.ID.to_numpy(str), EXP_EXPID=e.EXPID.to_numpy(np.int64), EXP_TILEID=e.TILEID.to_numpy(np.int64),
                     EXP_MJD=e.MJD.to_numpy(float), EXP_EXPTIME=e.EXPTIME.to_numpy(float))
    # flux scale (Stars' request): the fibre r flux from the photometry, as F_lambda at 6400 A in 1e-17 erg/s/cm2/A
    # (1 nanomaggy = 3.631e-29 erg/s/cm2/Hz; F_lambda = F_nu c / lambda^2 = 2.66e-18 per nanomaggy)
    nmgy = sel.FIBERFLUX_R.to_numpy(float) if a.survey == 'desi' else 10 ** (-0.4 * (sel.FIB_r.to_numpy(float) - 22.5))
    flux_scale = 0.266 * nmgy
    os.makedirs(OUT, exist_ok=True); out = os.path.join(OUT, f'galaxy_{a.survey}.npz'); tmp = out + '.tmp.npz'
    np.savez(tmp, WAVE=SWAVE.astype(np.float64), R=R, W=W, ID=sel.ID.to_numpy(str), CLASS=cls, SURVEY=np.array([a.survey] * n),
             Z=sel.Z.to_numpy(float), SNR=sel.SNR_SDSSPIX.to_numpy(float), PLATE=plate, MJD=mjd, AIRMASS=np.full(n, np.nan),
             NVEC=a.nvec, TELLURIC_WINDOWS=np.array(TELLURIC, float), MAPPED_FROM='rest_store_v2', BASIS=os.path.basename(a.basis),
             LINE_MASK_KMS=1000.0, FLUX_SCALE=flux_scale, FLUX_SCALE_NOTE='fibre r-band photometry as F_lambda(6400 A), 1e-17 erg/s/cm2/A',
             EBV_SFD=sel.EBV.to_numpy(float), **extra)
    os.replace(tmp, out)
    npl = pd.Series(plate[plate >= 0]).value_counts()
    msg = (f'{out}: {n} galaxies fitted; classes ' + ', '.join(f'{c} {np.sum(cls == c)}' for c in np.unique(cls)) +
           f'; {len(npl)} plates/tiles, {int((npl >= 5).sum())} with >= 5 galaxies ({int(npl[npl >= 5].sum())} galaxies); {time.time() - T0:.0f}s')
    print(msg, flush=True); status(len(groups), len(groups), 'done', msg.split(': ', 1)[1][:200])


if __name__ == '__main__':
    main()
