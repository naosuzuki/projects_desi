#!/usr/bin/env python3
"""Steps 2-3 (SDSS): spectrophotometric calibration, RV frame and per-star measurements over the SDSS spPlate files
(DR17 v5_13_2 = BOSS/eBOSS/SEGUE-on-BOSS, DR8 = SDSS-I/II + SEGUE), for the stars of step 1 (sdss_stars_gaia.parquet,
Gaia parallax selection; --all-stars includes the rest). Reuses the quasar thread's spPlate reader pieces
(programs_quasars/01_passQ_sdss.py: sdss_scale, plate_files) and starlib.

Per plate (rows of the catalog on that plate/MJD):
  Step 2  scale: synthetic SDSS g, r of the pipeline spectrum vs PSF g, r (the pipelines calibrate to PSF fluxes of
          point sources; DR8 PSF mags come from the Gaia crossmatch) -> CAL_A (scale), CAL_B_EST (tilt, stored),
          CAL_DG/CAL_DR, CAL_OK; offsets to the plugmap fiber mags stored (CAL_DG_FIBER/CAL_DR_FIBER);
          Gaia check from the raw spectrum (starlib.gaia_mags): SYN_GAIA_G/BP/RP, GAIA_DG/DBP/DRP/DBPRP, and GAIA_DG_CAL =
          the G residual after the PSF-mag scale (independent test of that scale).
  Step 3  RV frame: integer shift by the pipeline z on the log grid (spPlates are already on dlog = 1e-4: no rebin),
          RV_KMS, RV_RESID_KMS; indices, band ratios, colors (starlib); S/N per half of the spectrum.
  3b/4    distance, DM, ABS_SCALE, EBV_SFD (plugmap), EBV_3D (from the catalog when step3b_ebv3d.py has run), EBV_USED,
          A_G_USED, M_G_SYN0 (starlib.distance_columns).
Output: OUTDIR/<dr17|dr8>/spPlate-PPPP-MMMMM.parquet (resumable; empty when the plate has no selected star).

Usage:
  python3 step2_3_pass_sdss.py --survey dr17 --nfiles 3 --seed 1 --outdir ~/data/sdssdesi_passS_sdss_test
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 step2_3_pass_sdss.py --survey all --all > ~/data/sdssdesi_status/passS_sdss.log 2>&1 &
"""
import os, sys, time, argparse, warnings, gc, importlib.util
import numpy as np, pandas as pd, fitsio

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
QDIR = os.path.join(HERE, '..', 'programs_quasars'); sys.path.insert(0, QDIR)
_spec = importlib.util.spec_from_file_location('passQ_sdss', os.path.join(QDIR, '01_passQ_sdss.py'))
Qs = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(Qs)

OUTDIR = os.path.expanduser('~/data/sdssdesi_passS_sdss')
CATALOG = os.path.expanduser('~/data/sdssdesi_cache/sdss_stars_gaia.parquet')
DLOG = starlib.DLOG
CAT_COLS = ['SURVEY', 'PLATE', 'MJD', 'FIBERID', 'SUBCLASS', 'SUBCLASS1', 'Z', 'ZERR', 'ZWARN', 'PLATEQUALITY', 'PROGRAMNAME', 'SNR_CAT',
            'PSFMAG_U', 'PSFMAG_G', 'PSFMAG_R', 'PSFMAG_I', 'PSFMAG_Z', 'THING_ID', 'RA', 'DEC',
            'ELODIE_TEFF', 'ELODIE_LOGG', 'ELODIE_FEH', 'ELODIE_Z', 'ELODIE_SPTYPE',
            'GAIA_SOURCE_ID', 'GAIA_PARALLAX', 'GAIA_PARALLAX_OVER_ERROR', 'GAIA_RUWE', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG',
            'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_BP_RP', 'GAIA_RADIAL_VELOCITY', 'GAIA_HAS_XP_CONTINUOUS', 'GAIA_NON_SINGLE_STAR',
            'GAIA_PHOT_VARIABLE_FLAG', 'GAIA_TEFF_GSPPHOT', 'GAIA_LOGG_GSPPHOT', 'GAIA_MH_GSPPHOT', 'GAIA_AG_GSPPHOT', 'GAIA_EBPMINRP_GSPPHOT',
            'GAIA_PLX', 'PLX_CORR', 'QUALITY', 'M_G', 'M_G0', 'BPRP0', 'RV_PIPE_KMS']
OPT_COLS = ['EBV_B19', 'EBV_E24', 'EBV_3D', 'EBV_3D_SRC', 'EBV_3D_ERR']


def load_catalog(all_stars):
    import pyarrow.parquet as pq
    have = pq.ParquetFile(CATALOG).schema_arrow.names
    cols = [c for c in CAT_COLS + OPT_COLS if c in have]
    cat = pd.read_parquet(CATALOG, columns=cols)
    if not all_stars:
        cat = cat[cat.GAIA_PLX]
    cat['SPECTYPE'] = 'STAR'
    return cat.reset_index(drop=True)


def process_plate(path, cat, survey):
    t0 = time.time()
    base = os.path.basename(path).replace('.fits', '')
    _, plate, mjd = base.split('-'); plate, mjd = int(plate), int(mjd)
    rows = cat[(cat.PLATE == plate) & (cat.MJD == mjd)]
    if len(rows) == 0:
        return pd.DataFrame(), 0.0
    F = fitsio.FITS(path)
    h = F[0].read_header(); coeff0, coeff1 = float(h['COEFF0']), float(h['COEFF1'])
    assert abs(coeff1 - DLOG) < 1e-9, f'{path}: COEFF1={coeff1}'
    rid0 = int(round(coeff0 / DLOG))
    fib = rows.FIBERID.to_numpy() - 1
    flux = F[0].read()[fib].astype(float); ivar = F[1].read()[fib].astype(float); mask = F[2].read()[fib]
    plug = F[5].read(columns=['FIBERID', 'MAG', 'SFD_EBV', 'OBJTYPE'])
    F.close()
    order = np.argsort(plug['FIBERID']); plug = plug[order][fib]
    assert np.all(plug['FIBERID'] == rows.FIBERID.to_numpy())
    t_read = time.time() - t0
    df = rows.reset_index(drop=True).copy()
    df['EBV'] = plug['SFD_EBV'].astype(float)
    mag = plug['MAG'].astype(float)
    df['FIBERMAG_G'], df['FIBERMAG_R'] = mag[:, 1], mag[:, 2]
    df['OBJTYPE'] = np.char.strip(plug['OBJTYPE'].astype(str))
    npix = flux.shape[1]; w = 10.0**((rid0 + np.arange(npix)) * DLOG)
    good = (ivar > 0) & (mask == 0)
    with np.errstate(all='ignore'):
        snr = np.where(good, flux * np.sqrt(ivar), np.nan)
        df['SNR_ALL'] = np.nanmedian(snr, axis=1); df['NGOOD'] = good.sum(axis=1)
        for cam, (lo, hi) in {'B': (3600, 6000), 'R': (6000, 10400)}.items():
            s = (w >= lo) & (w < hi); df[f'SNR_{cam}'] = np.nanmedian(snr[:, s], axis=1)
    del snr
    # --- step 2: scale to the PSF magnitudes
    pg = np.where((df.PSFMAG_G > 10) & (df.PSFMAG_G < 30), df.PSFMAG_G, np.nan).astype(float)
    pr = np.where((df.PSFMAG_R > 10) & (df.PSFMAG_R < 30), df.PSFMAG_R, np.nan).astype(float)
    a, b_est, dg, dr, mg, mr, calok = Qs.sdss_scale(w, flux, ivar, good, pg, pr)
    df['SYN_G_RAW'], df['SYN_R_RAW'] = mg, mr
    fibg = np.where((df.FIBERMAG_G > 10) & (df.FIBERMAG_G < 30), df.FIBERMAG_G, np.nan).astype(float)
    fibr = np.where((df.FIBERMAG_R > 10) & (df.FIBERMAG_R < 30), df.FIBERMAG_R, np.nan).astype(float)
    df['CAL_DG_FIBER'], df['CAL_DR_FIBER'] = fibg - mg, fibr - mr
    df['CAL_A'], df['CAL_B'], df['CAL_B_EST'] = a, 0.0, b_est
    df['CAL_DG'], df['CAL_DR'], df['CAL_OK'] = dg, dr, calok
    df['CAL_MODE'] = np.where(calok, 1, 0).astype(np.int8)
    # synthetic SDSS u g r i z of the raw spectrum for the joint fit of step 2b (u is almost entirely extrapolated: SDSS spectra
    # start at 3800 A; its SYN_SDSS_EXTRAP_U says so and step2b drops bands with > 30 % padding)
    for k, v in starlib.synth_mags(w, flux, ivar, mask, starlib.SDSS_FILTERS, 'SYN_SDSS').items():
        df[k] = v
    # --- steps 2 (Gaia check), 3, 3b/4 (starlib)
    df = starlib.measure_stars(df, w, flux, ivar, mask, rid0=rid0)
    with np.errstate(all='ignore'):
        df['GAIA_DG_CAL'] = df['GAIA_DG'] + 2.5 * np.log10(df['CAL_A'])      # G residual after the PSF-mag scale
    df['FILE'] = base; df['T_READ'] = t_read; df['PROC_SEC'] = time.time() - t0
    return df, t_read


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
    ap.add_argument('--all-stars', action='store_true', help='include the stars without a Gaia parallax')
    a = ap.parse_args()
    surveys = ['dr17', 'dr8'] if a.survey == 'all' else [a.survey]
    t0 = time.time(); done = 0; nst = 0
    cat_all = load_catalog(a.all_stars)
    for survey in surveys:
        cat = cat_all[cat_all.SURVEY == survey]
        print(f'{survey}: {len(cat):,} stars in the catalog, {time.time()-t0:.0f}s', flush=True)
        paths = a.files if a.files else Qs.plate_files(survey)
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
                done += 1; nst += len(df)
                print(f'{survey} {k+1}/{len(paths)} {os.path.basename(p)} nstar={len(df)} read {t_read:.1f}s total {time.time()-t0:.0f}s', flush=True)
            except Exception as e:
                import traceback; traceback.print_exc()
                print('ERROR', p, repr(e), flush=True)
            gc.collect()
    print(f'done: {done} plates, {nst} stars, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
