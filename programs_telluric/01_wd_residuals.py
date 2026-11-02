#!/usr/bin/env python3
"""Telluric step 1 (Nao 2026-11-01): observed-frame residuals of white dwarfs per survey, for the telluric-only basis (02_basis.py).
Candidates: SDSS subclass WD / WDcooler / WDhotter or the Gaia WD box (M_G0 > 9, (BP-RP)_0 < 0.8), S/N >= 5; DESI Redrock SUBTYPE WD
with G < --gmax, BP-RP < 0.8 and M_G > 8 where the parallax allows it. No stellar basis, which may carry telluric structure itself:
per object a polynomial in log lambda (--deg) over 5600 A to the red end, fitted outside the bands of the external atmosphere
(DECam standard, tau > TAU_MIN) and outside the WD lines, with sigma clipping. R = f / continuum - 1, W = ivar x continuum^2, on SWAVE.
Spectra: spPlate rows (SDSS; already on the SWAVE lattice) or the native store (DESI; rows by seek, one flux-conserving rebin). No DESI coadd read.
FEATURELESS: no H or He I line (|EW| below LIMIT and EW error below ERRMAX at H-alpha, H-beta, He I 5877, 6680) and a clean continuum
(robust rms of R / sigma outside the bands below 1.3). Only 2 surveys of SDSS: dr8 = SDSS-I/II spPlates, dr17 = BOSS / eBOSS spPlates.
Output: /Volumes/exdisk3/sdssdesi_telluric/wd/wd_<survey>.npz (Phase 0 format: WAVE, R, W, ID, CLASS = DC | WD, SURVEY, Z, SNR, PLATE, MJD,
FLUX_SCALE, EW_*, FEATURELESS; DESI also EXP_ID, EXP_EXPID, EXP_TILEID, EXP_MJD, EXP_EXPTIME). PLATE on DESI: dominant TILEID.
Usage: python3 01_wd_residuals.py --survey dr8|dr17|desi [--workers 3] [--max N] [--gmax 18.8]
"""
import os, sys, zipfile, struct, argparse, time
from multiprocessing import Pool
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, 'programs_stars'), os.path.join(ROOT, 'programs_templates')]
import starlib, restframe
from telluric_guide import band_mask
passB = starlib.passB
SWAVE = starlib.SWAVE; LW = np.log10(SWAVE); DLOG = starlib.DLOG; SRID0 = int(round(LW[0] / DLOG))
OUT = '/Volumes/exdisk3/sdssdesi_telluric/wd'; NATIVE = '/Volumes/exdisk3/sdssdesi_native/desi'; EXPF = os.path.expanduser('~/data/sdssdesi_expfibermap')
LINES = {'HB': 4862.7, 'HA': 6564.6, 'HE5877': 5877.3, 'HE6680': 6680.0}; LIMIT = {'HB': 2.0, 'HA': 1.5, 'HE5877': 1.5, 'HE6680': 1.5}; ERRMAX = 1.5
MASKLINES = [6564.6, 5877.3, 6680.0, 7067.2, 8204.0, 8440.0, 8469.0, 8504.0, 8547.0, 8600.0, 8667.0, 8752.0, 8865.0, 9017.0, 9232.0, 9548.0]
WLO = 5600.0; BAND = band_mask(SWAVE)


def line_ew(f, iv, c):
    """EW [A] and its error of each line against the sidebands 25-60 A away (absorption > 0)."""
    out = {}
    for nm, l0 in LINES.items():
        d = np.abs(SWAVE - l0); core = d < 20; side = (d > 25) & (d < 60); dl = SWAVE * DLOG * np.log(10)
        if (iv[core] > 0).sum() < 0.7 * core.sum() or (iv[side] > 0).sum() < 0.7 * side.sum():
            out[nm] = (np.nan, np.nan); continue
        s = side & (iv > 0); cont = np.sum(f[s] * iv[s]) / np.sum(iv[s]); g = core & (iv > 0)
        out[nm] = (float(np.sum((1 - f[g] / cont) * dl[g])), float(np.sqrt(np.sum(dl[g]**2 / iv[g])) / abs(cont)) if cont != 0 else np.nan)
    return out


def one(f, iv, deg):
    """f, iv on SWAVE -> R, W, SNR, FLUX_SCALE, EW, FEATURELESS, or None."""
    ok = (iv > 0) & np.isfinite(f); red = ok & (SWAVE > WLO)
    if red.sum() < 500:
        return None
    x = (LW - LW[red].mean()) / 0.1; fit = red & ~BAND
    for l0 in MASKLINES: fit &= np.abs(SWAVE - l0) > 30
    if fit.sum() < 300:
        return None
    for _ in range(4):
        sw = np.sqrt(iv[fit]); co = np.linalg.lstsq(np.vander(x[fit], deg + 1) * sw[:, None], f[fit] * sw, rcond=None)[0]
        c = np.polyval(co, x); chi = (f - c) * np.sqrt(iv); new = fit & (np.abs(chi) < 3)
        if new.sum() < 300 or new.sum() == fit.sum():
            break
        fit = new
    if not (c[red] > 0).all():
        return None
    R = np.where(red, f / np.where(red, c, 1) - 1, np.nan); W = np.where(red, iv * c**2, 0.0)
    clean = red & ~BAND; rms = 1.4826 * np.median(np.abs(R[clean] * np.sqrt(W[clean])))
    ew = line_ew(f, iv, c); fl = rms < 1.3
    for nm, (e, de) in ew.items():
        fl = fl and np.isfinite(e) and abs(e) < LIMIT[nm] and de < ERRMAX
    snr = float(np.median(np.sqrt(W[clean])))
    return R.astype(np.float32), W.astype(np.float32), snr, float(np.median(c[red])), ew, bool(fl)


def plate_rows(args):
    path, g, survey, deg = args; res = []
    F = fitsio.FITS(path); h = F[0].read_header(); coeff0 = float(h['COEFF0']); rid0 = int(round(coeff0 / DLOG)); npix = int(h['NAXIS1'])
    good = restframe.sdss_badpix_ok(coeff0, npix, survey)
    for r in g.itertuples():
        i = int(r.FIBERID) - 1; fl = F[0][i:i + 1, :][0].astype(float); iv = F[1][i:i + 1, :][0].astype(float); mk = F[2][i:i + 1, :][0]
        iv = np.where((mk == 0) & good, iv, 0.0); f2 = np.zeros(len(SWAVE)); iv2 = np.zeros(len(SWAVE))
        s = np.arange(npix) + rid0 - SRID0; k = (s >= 0) & (s < len(SWAVE)); f2[s[k]] = fl[k]; iv2[s[k]] = iv[k]
        o = one(f2, iv2, deg)
        if o is not None:
            res.append((f'{r.PLATE}-{r.MJD}-{r.FIBERID}', float(r.Z), int(r.PLATE), float(r.MJD)) + o)
    F.close(); return res


def npy_offset(fh, info):
    fh.seek(info.header_offset); hdr = fh.read(30); nl, xl = struct.unpack('<HH', hdr[26:30]); fh.seek(info.header_offset + 30 + nl + xl)
    v = np.lib.format.read_magic(fh); shape, _, dtype = np.lib.format._read_array_header(fh, v); return fh.tell(), shape, dtype


def native_rows(args):
    prog, hpx, g, deg = args; p = f'{NATIVE}/coadd-main-{prog}-{hpx}.npz'; res = []
    if not os.path.exists(p):
        return res
    z = zipfile.ZipFile(p); tid = np.load(p)['TARGETID']; pos = {t: i for i, t in enumerate(tid)}
    with open(p, 'rb') as fh:
        of, shape, dt = npy_offset(fh, z.getinfo('FLUX.npy')); oi, _, _ = npy_offset(fh, z.getinfo('IVAR.npy')); npix = shape[1]; nb = npix * dt.itemsize
        w = 3600.0 + 0.8 * np.arange(npix)
        for r in g.itertuples():
            i = pos.get(r.TARGETID)
            if i is None:
                continue
            fh.seek(of + i * nb); f = np.frombuffer(fh.read(nb), dt).astype(float); fh.seek(oi + i * nb); iv = np.frombuffer(fh.read(nb), dt).astype(float)
            rid, rf, ri = passB.rebin(w, f[None], iv[None], (iv <= 0)[None].astype(np.int8)); f2 = np.zeros(len(SWAVE)); iv2 = np.zeros(len(SWAVE))
            s = rid - SRID0; k = (s >= 0) & (s < len(SWAVE)); f2[s[k]] = rf[0][k]; iv2[s[k]] = ri[0][k]
            o = one(f2, iv2, deg)
            if o is not None:
                res.append((f'{r.TARGETID}-{prog}', float(r.Z), -1, float(r.MEAN_MJD)) + o + (int(hpx),))
    return res


def select(survey, gmax):
    if survey == 'desi':
        z = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_cache/zall_stars_main.parquet'), columns=['TARGETID', 'SUBTYPE', 'PROGRAM', 'HEALPIX', 'Z', 'MEAN_MJD', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'PARALLAX'])
        z = z[(z.SUBTYPE.astype(str).str.strip() == 'WD') & z.PROGRAM.isin(['bright', 'dark'])]
        G = z.GAIA_PHOT_G_MEAN_MAG; col = z.GAIA_PHOT_BP_MEAN_MAG - z.GAIA_PHOT_RP_MEAN_MAG
        with np.errstate(all='ignore'):
            mg = G + 5 * np.log10(np.where(z.PARALLAX > 0, z.PARALLAX, np.nan) / 100)
        return z[(G > 0) & (G < gmax) & (col < 0.8) & ~((z.PARALLAX > 1) & (mg < 8))].sort_values(['PROGRAM', 'HEALPIX'])
    t = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_products/sdss_stars_gaia.parquet'), columns=['SURVEY', 'PLATE', 'MJD', 'FIBERID', 'SUBCLASS', 'Z', 'SNR_CAT', 'M_G0', 'BPRP0'])
    t = t[t.SURVEY == survey]; sub = t.SUBCLASS.astype(str).str.strip()
    return t[(sub.isin(['WD', 'WDcooler', 'WDhotter']) | ((t.M_G0 > 9) & (t.BPRP0 < 0.8))) & ~sub.isin(['WDmagnetic', 'CarbonWD', 'CalciumWD']) & (t.SNR_CAT >= 5)]


def desi_exposures(ids, hpx):
    """per object the IN_COADD exposures (flat lists as the galaxy thread's) and the dominant TILEID; expfibermap tables, local."""
    d = pd.DataFrame({'ID': ids, 'HPX': hpx}); d['TARGETID'] = d.ID.str.split('-').str[0].astype(np.int64); d['PROGRAM'] = d.ID.str.split('-').str[1]; parts = []
    for (prog, h), g in d.groupby(['PROGRAM', 'HPX']):
        p = os.path.join(EXPF, prog, f'hpx{h}.parquet')
        if not os.path.exists(p):
            continue
        e = pd.read_parquet(p, columns=['TARGETID', 'EXPID', 'TILEID', 'MJD', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z'])
        e = e[e.TARGETID.isin(g.TARGETID) & (e.IN_COADD_B | e.IN_COADD_R | e.IN_COADD_Z) & (e.EXPTIME > 0)].drop_duplicates(['TARGETID', 'EXPID'])
        parts.append(e.merge(g[['ID', 'TARGETID']], on='TARGETID'))
    e = pd.concat(parts, ignore_index=True)
    top = e.groupby(['ID', 'TILEID']).EXPTIME.sum().reset_index().sort_values('EXPTIME').drop_duplicates('ID', keep='last').set_index('ID').TILEID
    return top.reindex(ids).fillna(-1).to_numpy(np.int64), e


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', required=True, choices=['dr8', 'dr17', 'desi']); ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--max', type=int, default=0); ap.add_argument('--gmax', type=float, default=18.8); ap.add_argument('--deg', type=int, default=5)
    a = ap.parse_args(); t0 = time.time(); s = select(a.survey, a.gmax)
    if a.max: s = s.sample(a.max, random_state=1)
    print(f'{a.survey}: {len(s):,} candidates', flush=True)
    if a.survey == 'desi':
        jobs = [(p, h, g, a.deg) for (p, h), g in s.groupby(['PROGRAM', 'HEALPIX'])]; fn = native_rows
    else:
        jobs = []
        for (pl, mj), g in s.groupby(['PLATE', 'MJD']):
            base = f'spPlate-{int(pl):04d}-{int(mj)}.fits'; dirs = [f'~/sdssredux/v5_13_2'] if a.survey == 'dr17' else ['~/sdssredux/DR8/26', '~/sdssredux/DR8/103', '~/sdssredux/DR8/104']
            for d in dirs:
                p = os.path.join(os.path.expanduser(d), f'{int(pl):04d}', base)
                if os.path.exists(p):
                    jobs.append((p, g, a.survey, a.deg)); break
        fn = plate_rows
    print(f'{len(jobs):,} files', flush=True); rows = []
    with Pool(a.workers) as P:
        for k, r in enumerate(P.imap_unordered(fn, jobs, chunksize=4)):
            rows += r
            if (k + 1) % 500 == 0: print(f'  {k + 1:,}/{len(jobs):,} files, {len(rows):,} objects, {time.time() - t0:.0f} s', flush=True)
    n = len(rows); R = np.empty((n, len(SWAVE)), np.float32); W = np.empty((n, len(SWAVE)), np.float32)
    for i, r in enumerate(rows): R[i] = r[4]; W[i] = r[5]
    ids = np.array([r[0] for r in rows]); fl = np.array([r[9] for r in rows]); ex = {}
    plate = np.array([r[2] for r in rows], np.int64)
    if a.survey == 'desi':
        plate, e = desi_exposures(ids, np.array([r[10] for r in rows])); ex = dict(EXP_ID=e.ID.to_numpy(str), EXP_EXPID=e.EXPID.to_numpy(np.int64), EXP_TILEID=e.TILEID.to_numpy(np.int64), EXP_MJD=e.MJD.to_numpy(float), EXP_EXPTIME=e.EXPTIME.to_numpy(float))
    ew = {f'EW_{nm}': np.array([r[8][nm][0] for r in rows]) for nm in LINES}; ew.update({f'EWERR_{nm}': np.array([r[8][nm][1] for r in rows]) for nm in LINES})
    os.makedirs(OUT, exist_ok=True); out = f'{OUT}/wd_{a.survey}{"_test" if a.max else ""}.npz'; tmp = out + '.tmp.npz'
    np.savez(tmp, WAVE=SWAVE.astype(np.float64), R=R, W=W, ID=ids, CLASS=np.where(fl, 'DC', 'WD'), SURVEY=np.array([a.survey] * n), Z=np.array([r[1] for r in rows]),
             SNR=np.array([r[6] for r in rows]), PLATE=plate, MJD=np.array([r[3] for r in rows]), AIRMASS=np.full(n, np.nan), FLUX_SCALE=np.array([r[7] for r in rows]), FEATURELESS=fl, DEG=a.deg, **ew, **ex)
    os.replace(tmp, out)
    print(f'{out}: {n:,} objects, {int(fl.sum()):,} featureless; S/N of featureless {np.percentile(np.array([r[6] for r in rows])[fl], [10, 50, 90]).round(1) if fl.any() else "-"}; {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main()
