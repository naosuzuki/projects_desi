#!/usr/bin/env python3
"""Telluric Phase 0, stars (docs/PLAN_telluric.md; Nao 2026-11-01): observed-frame residuals of the cached v2 individual stars
for the joint telluric basis. No new DESI / SDSS read: the S/N > 20 individuals of step 7 (individuals_<tag>_5_20.npz,
dereddened absolute F_lambda at 10 pc on the rest grid SWAVE, with ivar).
Per star: weighted LSQ fit with its class's v2 basis (MEAN + EIGVEC, step7_empca.star_class of the star's stack) on the basis
pixels (PIXMASK) OUTSIDE the observed telluric windows, so the fit is not pulled by telluric residuals and they stay in
R = f / model - 1;  W = ivar x model^2 (0 outside PIXMASK). R and W are then moved to the observed frame
(log lambda_obs = log lambda_rest + log10(1 + z), linear interpolation on the log grid; the common grid is SWAVE itself).
Metadata: SDSS PLATE / MJD from the SID; DESI the TILEID with the most EXPTIME among the coadd's exposures and the
EXPTIME-weighted mean MJD (expfibermap tables, IN_COADD in any camera). AIRMASS is NaN until Phase 1.
Output (common Phase 0 format): /Volumes/exdisk3/sdssdesi_telluric/phase0/stars_<survey>.npz with WAVE, R, W [n, 4447] float32,
ID, CLASS, SURVEY, Z, SNR, PLATE, MJD, AIRMASS, and NVEC / TELLURIC_WINDOWS for the record.
Usage: python3 telluric_phase0.py --survey desi|sdss [--nvec 10]
"""
import os, sys, argparse
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
from step7_empca import star_class

V2 = os.path.expanduser('~/data/sdssdesi_stars/v2'); EXPF = os.path.expanduser('~/data/sdssdesi_expfibermap')
OUT = '/Volumes/exdisk3/sdssdesi_telluric/phase0'
TELLURIC = [(6855, 6960), (7150, 7350), (7580, 7720), (8110, 8370), (8920, 9800)]     # observed A (the plan's windows)
SWAVE = starlib.SWAVE; LW = np.log10(SWAVE); DLOG = starlib.DLOG


def redshifts(survey, ids):
    if survey == 'desi':
        t = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_products/stars/step6_selection_desi_0.5x0.1_gspphot.parquet'), columns=['SID', 'Z'])
    else:
        t = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_products/sdss_stars_gaia.parquet'), columns=['PLATE', 'MJD', 'FIBERID', 'Z'])
        t['SID'] = t.PLATE.astype(str) + '-' + t.MJD.astype(str) + '-' + t.FIBERID.astype(str)
    return pd.DataFrame({'SID': ids}).merge(t[['SID', 'Z']].drop_duplicates('SID'), on='SID', how='left').Z.to_numpy(float)


def desi_tiles(ids, files):
    """dominant TILEID (most EXPTIME) and EXPTIME-weighted mean MJD per star, from the expfibermap tables of its healpix."""
    d = pd.DataFrame({'SID': ids, 'FILE': files})
    d['TARGETID'] = d.SID.str.split('-').str[0].astype(np.int64); d['PROGRAM'] = d.SID.str.split('-').str[1]
    d['HPX'] = d.FILE.str.extract(r'-(\d+)\.fits$')[0].astype(int)
    parts = []
    for (prog, hpx), g in d.groupby(['PROGRAM', 'HPX']):
        p = os.path.join(EXPF, prog, f'hpx{hpx}.parquet')
        if not os.path.exists(p):
            continue
        e = pd.read_parquet(p, columns=['TARGETID', 'TILEID', 'MJD', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z'])
        e = e[e.TARGETID.isin(g.TARGETID) & (e.IN_COADD_B | e.IN_COADD_R | e.IN_COADD_Z) & (e.EXPTIME > 0)]
        e['PROGRAM'] = prog; parts.append(e)
    e = pd.concat(parts, ignore_index=True)
    e['WM'] = e.MJD * e.EXPTIME
    mjd = e.groupby(['TARGETID', 'PROGRAM']).agg(WM=('WM', 'sum'), EXP=('EXPTIME', 'sum')); mjd['MJD'] = mjd.WM / mjd.EXP
    tile = e.groupby(['TARGETID', 'PROGRAM', 'TILEID']).EXPTIME.sum().reset_index().sort_values('EXPTIME').drop_duplicates(['TARGETID', 'PROGRAM'], keep='last')
    d = d.merge(tile[['TARGETID', 'PROGRAM', 'TILEID']], on=['TARGETID', 'PROGRAM'], how='left').merge(mjd[['MJD']].reset_index(), on=['TARGETID', 'PROGRAM'], how='left')
    return d.TILEID.fillna(-1).to_numpy(np.int64), d.MJD.to_numpy(float)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', default='desi', choices=['desi', 'sdss']); ap.add_argument('--nvec', type=int, default=10)
    a = ap.parse_args(); tag = f'{a.survey}_0.5x0.1_gspphot'
    z0 = np.load(f'{V2}/individuals_{tag}_5_20.npz')
    F, IV, ids, stk = z0['FLUX'].astype(float), z0['W'].astype(float), z0['IDS'].astype(str), np.char.strip(z0['STACKS'].astype(str))
    cells = fitsio.read(f'{V2}/stacks_{tag}.fits', ext='CELLS')
    cls_of = {s.strip(): star_class(float(m), float(c)) for s, m, c in zip(cells['STACK'].astype(str), cells['MG0_MED'], cells['COL0_MED'])}
    cls = np.array([cls_of.get(s, '?') for s in stk])
    z = redshifts(a.survey, ids)
    mem = pd.read_csv(f'{V2}/members_{tag}.csv', usecols=['SID', 'FILE', 'SNR_ALL']).drop_duplicates('SID').set_index('SID')
    snr = mem.SNR_ALL.reindex(ids).to_numpy(float)
    n = len(ids); R = np.full((n, len(SWAVE)), np.nan, np.float32); W = np.zeros((n, len(SWAVE)), np.float32)
    nfit = 0
    for c in np.unique(cls):
        p = f'{V2}/empca_{tag}_{c}.fits'
        if not os.path.exists(p):
            print(f'  no basis for class {c}: {np.sum(cls == c)} stars skipped'); continue
        B = np.vstack([fitsio.read(p, ext='MEAN')[None], fitsio.read(p, ext='EIGVEC')[:a.nvec]]).astype(float); pm = fitsio.read(p, ext='PIXMASK').astype(bool)
        for i in np.where(cls == c)[0]:
            if not np.isfinite(z[i]):
                continue
            wo = SWAVE * (1 + z[i]); tel = np.zeros(len(SWAVE), bool)
            for lo, hi in TELLURIC: tel |= (wo >= lo) & (wo <= hi)
            ok = pm & (IV[i] > 0) & np.isfinite(F[i]); fit = ok & ~tel
            if fit.sum() < 200:
                continue
            sw = np.sqrt(IV[i, fit]); coef = np.linalg.lstsq((B[:, fit] * sw).T, F[i, fit] * sw, rcond=None)[0]
            model = coef @ B
            good = ok & (model > 0)
            r = np.where(good, F[i] / np.where(good, model, 1) - 1, np.nan); wr = np.where(good, IV[i] * model**2, 0.0)
            sh = np.log10(1 + z[i])                         # to the observed frame on the same log grid
            R[i] = np.interp(LW, LW + sh, r, left=np.nan, right=np.nan)
            wi = np.interp(LW, LW + sh, wr, left=0.0, right=0.0); bad = np.interp(LW, LW + sh, (~good).astype(float), left=1.0, right=1.0) > 0
            W[i] = np.where(bad, 0.0, wi); R[i, bad] = np.nan; nfit += 1
    if a.survey == 'sdss':
        plate = np.array([int(s.split('-')[0]) for s in ids], np.int64); mjd = np.array([float(s.split('-')[1]) for s in ids])
    else:
        plate, mjd = desi_tiles(ids, mem.FILE.reindex(ids).to_numpy(str))
    os.makedirs(OUT, exist_ok=True); out = os.path.join(OUT, f'stars_{a.survey}.npz'); tmp = out + '.tmp.npz'
    np.savez(tmp, WAVE=SWAVE.astype(np.float64), R=R, W=W, ID=ids, CLASS=cls, SURVEY=np.array([a.survey] * n), Z=z, SNR=snr,
             PLATE=plate, MJD=mjd, AIRMASS=np.full(n, np.nan), NVEC=a.nvec, TELLURIC_WINDOWS=np.array(TELLURIC, float))
    os.replace(tmp, out)
    npl = pd.Series(plate[plate >= 0]).value_counts()
    print(f'{out}: {n} stars, {nfit} fitted; classes ' + ', '.join(f'{c} {np.sum(cls == c)}' for c in np.unique(cls))
          + f'; {len(npl)} plates/tiles, {int((npl >= 5).sum())} with >= 5 stars ({int(npl[npl >= 5].sum())} stars)')


if __name__ == '__main__':
    main()
