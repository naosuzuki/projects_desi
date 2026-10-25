#!/usr/bin/env python3
"""Pass Q, SDSS driver: quasar measurements over the SDSS spPlate files (DR17 v5_13_2 = BOSS/eBOSS, and DR8 = SDSS-I/II).
Runs first (internal disk, so it does not compete with pass B on the external disks).

Per plate: rows with CLASS == 'QSO' and SPECPRIMARY == 1 in spAll (DR17) / specObj (DR8); ZWARNING is kept as a column.
  Step 1  spectrophotometric scale: synthetic SDSS g, r of the spectrum vs the PSF magnitudes (the SDSS pipeline
          calibrates to PSF fluxes of the standards, so a point source's spectrum matches its PSF flux, not its 2"/3"
          fiber flux); ivar-weighted mean of the two flux ratios, scale only (CAL_A; tilt stored as CAL_B_EST).
          DR17: PSFMAG from spAll. DR8: no PSF mags in specObj -> no scale (CAL_A = 1, CAL_MODE 0), fiber offsets stored.
          Quasars vary between the imaging and spectroscopic epochs, so the scatter is larger than for stars.
  Step 2  Milky Way extinction: F99, R_V = 3.1, E(B-V) = 0.86 x SFD_EBV (plugmap)
  Step 3  rest frame: the spPlate grid is already log10 with dlog = 1e-4 -> integer-pixel shift (no rebin), F x (1+z)
  Step 4  absolute magnitude M_I (i-band, z=2 normalization): m_i = PSFMAG[3] (DR17) or plugmap fiber MAG[3] - 0.3 (DR8;
          M_I_SRC = 2), dereddened with EXTINCTION or 1.698 x E(B-V); DM from flat LCDM h=0.7 Om=0.3; K from Richards+06
  Step 5  measurements: see qsolib.measure_rest (S/N, galaxy lines, Dn4000, slope, broad + oxygen lines, Z_SYS, forest, BAL)
Output: one Parquet per plate in OUTDIR/<dr17|dr8>/spPlate-PPPP-MMMMM.parquet (resumable; empty when no QSO rows).

Usage:
  python3 01_passQ_sdss.py --survey dr17 --nfiles 3 --seed 1
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 01_passQ_sdss.py --survey all --all > ~/data/sdssdesi_status/passQ_sdss.log 2>&1 &
"""
import os, sys, time, glob, argparse, warnings, gc
import numpy as np
import pandas as pd
import fitsio
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qsolib as Q

OUTDIR = os.path.expanduser('~/data/sdssdesi_passQ_sdss')
SPECTRO = os.path.expanduser('~/sdssredux')
SPALL = f'{SPECTRO}/spAll-v5_13_2.fits'
SPECOBJ = f'{SPECTRO}/specObj-dr8.fits'
KCORR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'quasar', 'richards06_table4_kcorr.dat')
PLATE_DIRS = {'dr17': [f'{SPECTRO}/v5_13_2'], 'dr8': [f'{SPECTRO}/DR8/26', f'{SPECTRO}/DR8/103', f'{SPECTRO}/DR8/104']}
SDSS_EXT = {'u': 4.239, 'g': 3.303, 'r': 2.285, 'i': 1.698, 'z': 1.263}     # A/E(B-V), Schlafly & Finkbeiner 2011 (R_V=3.1)
FIBER_TO_PSF_DR8 = -0.3     # mag: SDSS-I/II 3" fiber mag -> PSF mag for a point source (approximate; flagged M_I_SRC=2)
CHUNK = Q.passB.CHUNK
F32 = np.float32
_K = np.loadtxt(KCORR)


# ---------------------------------------------------------------- catalogs (QSO rows only, read once)
def load_catalog(survey):
    if survey == 'dr17':
        cols = ['PLATE', 'MJD', 'FIBERID', 'CLASS', 'Z', 'Z_ERR', 'ZWARNING', 'SPECPRIMARY', 'PLATEQUALITY', 'PROGRAMNAME',
                'SN_MEDIAN_ALL', 'PSFMAG', 'FIBER2MAG', 'EXTINCTION', 'THING_ID', 'RA', 'DEC']
        t = fitsio.read(SPALL, columns=cols)
    else:
        cols = ['PLATE', 'MJD', 'FIBERID', 'CLASS', 'Z', 'Z_ERR', 'ZWARNING', 'SPECPRIMARY', 'PLATEQUALITY', 'PROGRAMNAME',
                'SN_MEDIAN', 'PLUG_RA', 'PLUG_DEC']
        t = fitsio.read(SPECOBJ, columns=cols)
    cls = np.char.strip(t['CLASS'].astype(str))
    sel = (cls == 'QSO') & (t['SPECPRIMARY'] == 1)
    t = t[sel]
    df = pd.DataFrame({'PLATE': t['PLATE'].astype(np.int32), 'MJD': t['MJD'].astype(np.int32), 'FIBERID': t['FIBERID'].astype(np.int16),
                       'Z': t['Z'].astype(float), 'ZERR': t['Z_ERR'].astype(float), 'ZWARN': t['ZWARNING'].astype(np.int32),
                       'PLATEQUALITY': np.char.strip(t['PLATEQUALITY'].astype(str)), 'PROGRAMNAME': np.char.strip(t['PROGRAMNAME'].astype(str))})
    df['SPECTYPE'] = 'QSO'; df['ZCAT_PRIMARY'] = True; df['DELTACHI2'] = np.nan
    if survey == 'dr17':
        df['SNR_CAT'] = t['SN_MEDIAN_ALL'].astype(float)
        df['PSFMAG_I'] = t['PSFMAG'][:, 3].astype(float); df['EXT_I'] = t['EXTINCTION'][:, 3].astype(float)
        df['PSFMAG_G'] = t['PSFMAG'][:, 1].astype(float); df['PSFMAG_R'] = t['PSFMAG'][:, 2].astype(float)
        df['FIBER2MAG_G'] = t['FIBER2MAG'][:, 1].astype(float); df['FIBER2MAG_R'] = t['FIBER2MAG'][:, 2].astype(float)
        df['THING_ID'] = t['THING_ID'].astype(np.int64); df['TARGET_RA'] = t['RA'].astype(float); df['TARGET_DEC'] = t['DEC'].astype(float)
    else:
        sn = t['SN_MEDIAN']; df['SNR_CAT'] = np.nanmedian(np.where(sn > 0, sn, np.nan)[:, 1:4], axis=1) if sn.ndim == 2 else sn.astype(float)
        df['PSFMAG_I'] = np.nan; df['EXT_I'] = np.nan; df['PSFMAG_G'] = np.nan; df['PSFMAG_R'] = np.nan
        df['FIBER2MAG_G'] = np.nan; df['FIBER2MAG_R'] = np.nan
        df['THING_ID'] = -1; df['TARGET_RA'] = t['PLUG_RA'].astype(float); df['TARGET_DEC'] = t['PLUG_DEC'].astype(float)
    df['SURVEY'] = survey
    return df


def plate_files(survey):
    out = []
    for d in PLATE_DIRS[survey]:
        out += glob.glob(f'{d}/*/spPlate-*.fits')
    return sorted(out)


# ---------------------------------------------------------------- Step 1 for SDSS: scale to the fiber magnitudes
def sdss_scale(w, flux, ivar, good, fibmag_g, fibmag_r):
    """ivar-weighted mean of the catalog/synthetic flux ratios in g and r (sdss2010 filters); returns a, b_est, dg, dr."""
    n = flux.shape[0]
    fl = flux * Q.passB.FLUX_UNIT; iv = ivar / Q.passB.FLUX_UNIT**2
    mg, eg = Q.passB.ab_mags(Q.passB.filt('sdss2010-g'), w, fl, good, ivar2d=iv)
    mr, er = Q.passB.ab_mags(Q.passB.filt('sdss2010-r'), w, fl, good, ivar2d=iv)
    with np.errstate(all='ignore'):
        dg = fibmag_g - mg; dr = fibmag_r - mr
        Rg, Rr = 10**(-0.4 * dg), 10**(-0.4 * dr)
        wg = np.where(np.isfinite(Rg), 1.0 / np.where(np.isfinite(eg) & (eg > 0), eg, 0.3)**2, 0.0)
        wr = np.where(np.isfinite(Rr), 1.0 / np.where(np.isfinite(er) & (er > 0), er, 0.3)**2, 0.0)
        a = (np.nan_to_num(Rg) * wg + np.nan_to_num(Rr) * wr) / np.where(wg + wr > 0, wg + wr, np.nan)
        lg = Q.passB.filt('sdss2010-g').effective_wavelength.value; lr = Q.passB.filt('sdss2010-r').effective_wavelength.value
        xg, xr = np.log10(lg / 6000), np.log10(lr / 6000)
        b_est = np.where(np.isfinite(eg) & np.isfinite(er) & (eg < Q.passB.TILT_MAX_ERR) & (er < Q.passB.TILT_MAX_ERR), (Rr - Rg) / (xr - xg), np.nan)
    ok = np.isfinite(a) & (a > 0) & (np.abs(np.nan_to_num(dg)) < 0.75) & (np.abs(np.nan_to_num(dr)) < 0.75) & (np.isfinite(dg) | np.isfinite(dr))
    return np.where(ok, a, np.nan), b_est, dg, dr, mg, mr, ok


# ---------------------------------------------------------------- per plate
def process_plate(path, cat, survey):
    t0 = time.time()
    base = os.path.basename(path).replace('.fits', '')
    _, plate, mjd = base.split('-'); plate, mjd = int(plate), int(mjd)
    rows = cat[(cat.PLATE == plate) & (cat.MJD == mjd)]
    if len(rows) == 0:
        return pd.DataFrame(), 0.0
    F = fitsio.FITS(path)
    h = F[0].read_header(); coeff0, coeff1 = float(h['COEFF0']), float(h['COEFF1'])
    assert abs(coeff1 - Q.DLOG) < 1e-9, f'{path}: COEFF1={coeff1}'
    rid0 = int(round(coeff0 / Q.DLOG))
    fib = rows.FIBERID.to_numpy() - 1                     # row index in the plate arrays
    flux = F[0].read()[fib].astype(float); ivar = F[1].read()[fib].astype(float); mask = F[2].read()[fib]
    plug = F[5].read(columns=['FIBERID', 'MAG', 'SFD_EBV', 'OBJTYPE'])
    F.close()
    order = np.argsort(plug['FIBERID']); plug = plug[order][fib]
    assert np.all(plug['FIBERID'] == rows.FIBERID.to_numpy())
    t_read = time.time() - t0
    df = rows.reset_index(drop=True).copy()
    df['EBV'] = plug['SFD_EBV'].astype(float)
    mag = plug['MAG'].astype(float)
    df['FIBERMAG_G'], df['FIBERMAG_R'], df['FIBERMAG_I'] = mag[:, 1], mag[:, 2], mag[:, 3]
    df['OBJTYPE'] = np.char.strip(plug['OBJTYPE'].astype(str))
    npix = flux.shape[1]; w = 10.0**((rid0 + np.arange(npix)) * Q.DLOG)
    good = (ivar > 0) & (mask == 0)
    with np.errstate(all='ignore'):
        snr = np.where(good, flux * np.sqrt(ivar), np.nan)
        df['SNR_ALL'] = np.nanmedian(snr, axis=1); df['SNR_SDSSPIX'] = df['SNR_ALL']; df['NGOOD'] = good.sum(axis=1)
        for cam, (lo, hi) in {'B': (3600, 6000), 'R': (6000, 10400)}.items():
            s = (w >= lo) & (w < hi); df[f'SNR_{cam}'] = np.nanmedian(snr[:, s], axis=1)
    df['SNR_Z'] = np.nan
    del snr
    # Step 1 reference. The SDSS/BOSS pipeline calibrates spectra to the PSF magnitudes of the standard stars, so a
    # point source's spectrum reproduces its PSF flux, not the flux in the 2" (BOSS) / 3" (SDSS-I/II) fiber aperture
    # (test: synthetic - fiber mag = -0.70 for BOSS, -0.30 for SDSS-I/II). DR17: scale to PSFMAG g, r (spAll).
    # DR8: specObj carries no PSF magnitudes -> no scale applied (CAL_A = 1, CAL_MODE 0); the offsets to the plugmap
    # fiber magnitudes are stored (CAL_DG/DR_FIBER) for the record.
    fibg = np.where((df.FIBERMAG_G > 10) & (df.FIBERMAG_G < 30), df.FIBERMAG_G, np.nan).astype(float)
    fibr = np.where((df.FIBERMAG_R > 10) & (df.FIBERMAG_R < 30), df.FIBERMAG_R, np.nan).astype(float)
    pg = np.where((df.PSFMAG_G > 10) & (df.PSFMAG_G < 30), df.PSFMAG_G, np.nan).astype(float)
    pr = np.where((df.PSFMAG_R > 10) & (df.PSFMAG_R < 30), df.PSFMAG_R, np.nan).astype(float)
    a, b_est, dg, dr, mg, mr, calok = sdss_scale(w, flux, ivar, good, pg, pr)
    df['SYN_G_RAW'], df['SYN_R_RAW'] = mg, mr
    df['CAL_DG_FIBER'], df['CAL_DR_FIBER'] = fibg - mg, fibr - mr
    if survey != 'dr17':
        a = np.where(np.isfinite(mg) | np.isfinite(mr), 1.0, np.nan); calok = np.isfinite(a)
    df['CAL_A'], df['CAL_B'], df['CAL_B_EST'] = a, 0.0, b_est
    df['CAL_DG'], df['CAL_DR'], df['CAL_OK'] = dg, dr, calok
    df['CAL_MODE'] = np.where(calok & (survey == 'dr17'), 1, 0).astype(np.int8)
    # Step 4: M_I (i band, z=2 normalization à la Richards+06)
    ebv = df.EBV.to_numpy(dtype=float)
    m_i = np.where(np.isfinite(df.PSFMAG_I) & (df.PSFMAG_I > 0), df.PSFMAG_I - np.nan_to_num(df.EXT_I, nan=0.0), np.nan)
    src = np.where(np.isfinite(m_i), 1, 0).astype(np.int8)
    alt = df.FIBERMAG_I.to_numpy(dtype=float) + FIBER_TO_PSF_DR8 - SDSS_EXT['i'] * np.nan_to_num(ebv)
    m_i = np.where(np.isfinite(m_i), m_i, np.where((alt > 10) & (alt < 30), alt, np.nan)); src = np.where((src == 0) & np.isfinite(m_i), 2, src)
    z = df.Z.to_numpy(dtype=float); zok = np.isfinite(z) & (z > 0)
    DM = np.where(zok, Q.passB.COSMO.distmod(np.clip(z, 1e-4, None)).value, np.nan)
    df['DM'] = DM; df['M_I_OBS_DERED'] = m_i; df['M_I_SRC'] = src
    df['MABS_I_Z2'] = m_i - DM - np.interp(np.clip(z, 0, 5.49), _K[:, 0], _K[:, 1])
    # Steps 2-3 + 5 in chunks
    parts = []
    for i0 in range(0, len(df), CHUNK):
        sl = slice(i0, min(i0 + CHUNK, len(df)))
        rflux, rivar = Q.rest_from_loggrid(rid0, flux[sl], ivar[sl], mask[sl], z[sl], np.nan_to_num(a[sl], nan=1.0), ebv[sl])
        parts.append(Q.measure_rest(df.iloc[sl].copy(), rflux, rivar))
        del rflux, rivar; gc.collect()
    out = pd.concat(parts, ignore_index=True)
    out['FILE'] = base; out['T_READ'] = t_read; out['PROC_SEC'] = time.time() - t0
    return out, t_read


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='all', choices=['dr17', 'dr8', 'all'])
    ap.add_argument('--nfiles', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--files', nargs='*')
    ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--shard', default='0/1')
    a = ap.parse_args()
    surveys = ['dr17', 'dr8'] if a.survey == 'all' else [a.survey]
    t0 = time.time(); done = 0; nq = 0
    for survey in surveys:
        cat = load_catalog(survey)
        print(f'{survey}: {len(cat):,} primary QSO rows in the catalog, {time.time()-t0:.0f}s', flush=True)
        paths = a.files if a.files else plate_files(survey)
        if not a.files and not a.all:
            rng = np.random.default_rng(a.seed); paths = [paths[i] for i in rng.choice(len(paths), a.nfiles, replace=False)]
        si, sn = (int(v) for v in a.shard.split('/')); paths = paths[si::sn]
        for k, p in enumerate(paths):
            out = os.path.join(a.outdir, survey, os.path.basename(p).replace('.fits', '.parquet'))
            if os.path.exists(out) and not a.overwrite:
                continue
            os.makedirs(os.path.dirname(out), exist_ok=True)
            try:
                df, t_read = process_plate(p, cat, survey)
                df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
                done += 1; nq += len(df)
                print(f'{survey} {k+1}/{len(paths)} {os.path.basename(p)} nqso={len(df)} read {t_read:.1f}s total {time.time()-t0:.0f}s', flush=True)
            except Exception as e:
                import traceback; traceback.print_exc()
                print('ERROR', p, repr(e), flush=True)
    print(f'done: {done} plates, {nq} quasars, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
