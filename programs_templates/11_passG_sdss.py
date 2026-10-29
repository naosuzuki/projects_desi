#!/usr/bin/env python3
"""SDSS galaxy pass: Steps 1-6 for the SDSS DR8 (SDSS-I/II) and DR17 (BOSS/eBOSS) galaxy spectra, on the internal disk.

Same recipe and grid as the DESI pass C (08_passC_stack.py): vacuum line table, normalization windows (normwin.py with
the survey's observed range), 5-axis cells (gridcells.py), per-cell accumulators. SDSS galaxies are stacked into
SEPARATE accumulators per survey, so the combination with DESI is decided later (DR8 main sample = colour-unbiased
low-z anchor; BOSS = validation).

Step 1 (SDSS): the SDSS pipeline calibrates to PSF magnitudes, so a galaxy's spectrum carries its fiber light divided
by the point-source fiber fraction. Reference = fiber magnitude (DR8: FIBERMAG 3", photoPlate-dr8; DR17: FIBER2MAG 2",
spAll) + the plate's point-source offset DELTA_PLATE = median(PSFMAG - FIBERMAG) of the plate's stars (r band; g from
its own median), the analogue of DESI's MEAN_PSF_TO_FIBER_SPECFLUX. Scale-only, ivar-weighted over g and r; spectra
failing the +-0.5 mag check keep scale 1 (flag CALC=0) and are still stacked.
Step 2: F99, R_V=3.1, 0.86 x SFD (plugmap SFD_EBV). Step 3: spPlate grid is already log10 (dlog=1e-4) -> v2 (Nao 2026-10-29): sub-pixel Lanczos-3
shift at the exact z and COEFF0 (restframe.sdss_rest); v1 shifted by round(log10(1+z)/dlog) whole pixels (~/data/sdssdesi_passG_sdss).
Step 4: K from the spectrum (SDSS r observed vs rest DECam r/g, SDSS u, NUV tophat); M from CMODELMAG r dereddened
(A_r = 2.285 E(B-V), SF11). Colour axis: CMODELMAG g-r (z<0.5) / r-z (z>=0.5), dereddened with SF11 coefficients.
Step 5: lines (vacuum), Dn4000. Step 6: normalize, accumulate.

Usage:
  python3 11_passG_sdss.py --prep                                  # build the galaxy catalog + plate offsets (once, ~3 min)
  OPENBLAS_NUM_THREADS=1 python3 11_passG_sdss.py --survey dr8 --shard 0/2 [--store /Volumes/exdisk3/sdssdesi_restframe_sdss_v2]
Outputs: ~/data/sdssdesi_passG_sdss_v2/{tables/<survey>/<plate>.parquet, acc/acc_<survey>_shard<i>of<n>.npz}
"""
import os, sys, glob, time, json, argparse, signal, warnings
import numpy as np
import pandas as pd
import fitsio

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')
pc = importlib.import_module('08_passC_stack')
import normwin, restframe, gridcells as gc

assert pb.LINES_FRAME == 'vacuum'
SPECTRO = os.path.expanduser('~/sdssredux')
CAT = os.path.expanduser('~/data/sdssdesi_products/sdss_galaxy_catalog.parquet')
OFFS = os.path.expanduser('~/data/sdssdesi_products/sdss_plate_psf_offsets.parquet')
OUTDIR = os.path.expanduser('~/data/sdssdesi_passG_sdss_v2')   # v1 (whole-pixel shift): ~/data/sdssdesi_passG_sdss
PLATE_DIRS = {'dr17': [f'{SPECTRO}/v5_13_2'], 'dr8': [f'{SPECTRO}/DR8/26', f'{SPECTRO}/DR8/103', f'{SPECTRO}/DR8/104']}
OBS_RANGE = {'dr8': (3800.0, 9200.0), 'dr17': (3600.0, 10300.0)}
SF11 = {'u': 4.239, 'g': 3.303, 'r': 2.285, 'i': 1.698, 'z': 1.263}
NREST, RWAVE, DLOG = pb.NREST, pb.RWAVE, pb.DLOG
F32 = np.float32
STOP = False
signal.signal(signal.SIGTERM, lambda s, f: globals().__setitem__('STOP', True))
signal.signal(signal.SIGINT, lambda s, f: globals().__setitem__('STOP', True))


# ---------------------------------------------------------------- one-time catalog preparation
def prep():
    t0 = time.time()
    a = fitsio.read(f'{SPECTRO}/spAll-v5_13_2.fits', columns=['PLATE', 'MJD', 'FIBERID', 'CLASS', 'Z', 'ZWARNING', 'SPECPRIMARY',
                                                               'PLATEQUALITY', 'PSFMAG', 'FIBER2MAG', 'CMODELMAG', 'SN_MEDIAN_ALL'])
    d17 = pd.DataFrame({'PLATE': a['PLATE'], 'MJD': a['MJD'], 'FIBERID': a['FIBERID'], 'CLASS': np.char.strip(a['CLASS'].astype(str)),
                        'Z': a['Z'].astype(float), 'ZWARN': a['ZWARNING'], 'PRIMARY': a['SPECPRIMARY'] == 1,
                        'PQ': np.char.strip(a['PLATEQUALITY'].astype(str)), 'SNR_CAT': a['SN_MEDIAN_ALL'].astype(float)})
    for k, b in enumerate('ugriz'):
        d17[f'PSF_{b}'] = a['PSFMAG'][:, k]; d17[f'FIB_{b}'] = a['FIBER2MAG'][:, k]; d17[f'CMOD_{b}'] = a['CMODELMAG'][:, k]
    d17['SURVEY'] = 'dr17'; del a
    print(f'spAll read {time.time()-t0:.0f}s', flush=True)
    s = fitsio.read(f'{SPECTRO}/specObj-dr8.fits', columns=['PLATE', 'MJD', 'FIBERID', 'CLASS', 'Z', 'ZWARNING', 'SPECPRIMARY', 'PLATEQUALITY'])
    p = fitsio.read(f'{SPECTRO}/photoPlate-dr8.fits', columns=['PSFMAG', 'FIBERMAG', 'CMODELMAG'])     # row-aligned with specObj
    d8 = pd.DataFrame({'PLATE': s['PLATE'], 'MJD': s['MJD'], 'FIBERID': s['FIBERID'], 'CLASS': np.char.strip(s['CLASS'].astype(str)),
                       'Z': s['Z'].astype(float), 'ZWARN': s['ZWARNING'], 'PRIMARY': s['SPECPRIMARY'] == 1,
                       'PQ': np.char.strip(s['PLATEQUALITY'].astype(str)), 'SNR_CAT': np.nan})
    for k, b in enumerate('ugriz'):
        d8[f'PSF_{b}'] = p['PSFMAG'][:, k]; d8[f'FIB_{b}'] = p['FIBERMAG'][:, k]; d8[f'CMOD_{b}'] = p['CMODELMAG'][:, k]
    d8['SURVEY'] = 'dr8'; del s, p
    d = pd.concat([d17, d8], ignore_index=True); del d17, d8
    for c in [c for c in d.columns if c[:4] in ('PSF_', 'FIB_', 'CMOD')]:
        d[c] = d[c].where((d[c] > 10) & (d[c] < 30))
    # per-plate point-source offset from the plate's stars (PSF - fiber), r and g
    st = d[(d.CLASS == 'STAR') & (d.PSF_r > 14) & (d.PSF_r < 20.5)]
    o = st.groupby(['SURVEY', 'PLATE', 'MJD']).agg(N_STAR=('PSF_r', 'size'),
                                                   DELTA_R=('PSF_r', lambda x: np.nan), DELTA_G=('PSF_g', lambda x: np.nan)).reset_index()
    dr = (st.PSF_r - st.FIB_r).groupby([st.SURVEY, st.PLATE, st.MJD]).median().rename('DELTA_R')
    dg = (st.PSF_g - st.FIB_g).groupby([st.SURVEY, st.PLATE, st.MJD]).median().rename('DELTA_G')
    o = o.drop(columns=['DELTA_R', 'DELTA_G']).merge(dr.reset_index(), on=['SURVEY', 'PLATE', 'MJD']).merge(dg.reset_index(), on=['SURVEY', 'PLATE', 'MJD'])
    glob_med = o.groupby('SURVEY')[['DELTA_R', 'DELTA_G']].median()
    print('global median point-source offset PSF-fiber:\n', glob_med.round(3), flush=True)
    o.to_parquet(OFFS)
    g = d[(d.CLASS == 'GALAXY') & d.PRIMARY & (d.ZWARN == 0) & (d.PQ == 'good') & (d.Z > 0.001) & (d.Z < gc.ZMAX)]
    g = g.drop(columns=['CLASS', 'PRIMARY', 'PQ']).reset_index(drop=True)
    g.to_parquet(CAT)
    print(f'galaxies: {len(g):,} ({g.SURVEY.value_counts().to_dict()}); plates with offsets {len(o):,}; {time.time()-t0:.0f}s')


# ---------------------------------------------------------------- per plate
def process_plate(path, survey, cat, offs, cellmap, acc, store_dir):
    t0 = time.time()
    base = os.path.basename(path).replace('.fits', ''); _, plate, mjd = base.split('-'); plate, mjd = int(plate), int(mjd)
    rows = cat[(cat.PLATE == plate) & (cat.MJD == mjd)]
    if len(rows) == 0:
        return pd.DataFrame(), 0, time.time() - t0
    F = fitsio.FITS(path); h = F[0].read_header()
    assert abs(float(h['COEFF1']) - DLOG) < 1e-9
    coeff0 = float(h['COEFF0']); rid0 = int(round(coeff0 / DLOG))
    fib = rows.FIBERID.to_numpy() - 1
    flux = F[0].read()[fib].astype(float); ivar = F[1].read()[fib].astype(float); mask = F[2].read()[fib]
    plug = F[5].read(columns=['FIBERID', 'SFD_EBV']); F.close()
    plug = plug[np.argsort(plug['FIBERID'])][fib]
    df = rows.reset_index(drop=True).copy(); df['EBV'] = plug['SFD_EBV'].astype(float)
    n, npix = flux.shape; w = 10.0**(coeff0 + np.arange(npix) * DLOG)
    good = (ivar > 0) & (mask == 0)
    with np.errstate(all='ignore'):
        df['SNR_SDSSPIX'] = np.nanmedian(np.where(good, flux * np.sqrt(ivar), np.nan), axis=1)
    # Step 1: fiber magnitude + plate point-source offset
    o = offs.get((survey, plate, mjd))
    dR, dG = (o if o is not None else GLOBAL_DELTA[survey])
    fl = flux * pb.FLUX_UNIT; iv = ivar / pb.FLUX_UNIT**2
    mg, eg = pb.ab_mags(pb.filt('sdss2010-g'), w, fl, good, ivar2d=iv)
    mr, er = pb.ab_mags(pb.filt('sdss2010-r'), w, fl, good, ivar2d=iv)
    with np.errstate(all='ignore'):
        dg = (df.FIB_g.to_numpy() + dG) - mg; dr = (df.FIB_r.to_numpy() + dR) - mr
        Rg, Rr = 10**(-0.4 * dg), 10**(-0.4 * dr)
        wg = np.where(np.isfinite(Rg), 1 / np.where(np.isfinite(eg) & (eg > 0), eg, 0.3)**2, 0)
        wr = np.where(np.isfinite(Rr), 1 / np.where(np.isfinite(er) & (er > 0), er, 0.3)**2, 0)
        a = (np.nan_to_num(Rg) * wg + np.nan_to_num(Rr) * wr) / np.where(wg + wr > 0, wg + wr, np.nan)
    calok = np.isfinite(a) & (a > 0) & (np.abs(np.nan_to_num(dg)) < 0.5) & (np.abs(np.nan_to_num(dr)) < 0.5) & (np.isfinite(dg) | np.isfinite(dr))
    a = np.where(calok, a, 1.0)
    df['CAL_DG'], df['CAL_DR'], df['CAL_A'], df['CALC'], df['DELTA_PLATE_R'] = dg, dr, a, calok, dR
    # Step 2
    ebv = np.nan_to_num(df.EBV.to_numpy(float))
    A = (3.1 * 0.86 * ebv)[:, None] * pb.extinction.fitzpatrick99(w, 1.0, 3.1)[None, :]
    fder = fl * a[:, None] * 10**(0.4 * A); ider = np.where(good, iv / a[:, None]**2 * 10**(-0.8 * A), 0.0)
    m_obs_r, _ = pb.ab_mags(pb.filt('sdss2010-r'), w, fder, good)
    # Step 3: sub-pixel shift on the native log grid (Lanczos-3)
    z = df.Z.to_numpy(float)
    rflux, rivar = restframe.sdss_rest(coeff0, fder, ider, (~good).astype(np.int8), z, pb.RID_MIN, NREST)
    rflux *= (1 + z)[:, None].astype(F32); rivar /= ((1 + z)**2)[:, None].astype(F32)
    lobs = RWAVE[None, :] * (1 + z)[:, None]
    rin = (lobs >= w[0]) & (lobs <= w[-1])
    del fder, ider, fl, iv, lobs
    # Step 4
    DM = pb.COSMO.distmod(np.clip(z, 1e-4, None)).value
    m_tot = df.CMOD_r.to_numpy(float) - SF11['r'] * ebv
    for b, nm in {'R': 'decam2014-r', 'G': 'decam2014-g', 'U': 'sdss2010-u', 'NUV': 'nuv'}.items():
        mrest, _ = pb.ab_mags(pb.filt(nm), RWAVE, rflux, rivar > 0, inrange2d=rin)
        df[f'MABS_{b}'] = m_tot - DM - (m_obs_r - mrest)
    # Step 5
    for k, v in pb.fit_lines(RWAVE, rflux, rivar).items():
        df[k] = v
    df['DN4000'], df['DN4000_ERR'] = pb.dn4000(RWAVE, rflux, rivar)
    df['MABS'] = gc.pick_mabs(z, df.MABS_R, df.MABS_G, df.MABS_U, df.MABS_NUV); df['MBAND'] = gc.mband(z)
    cg = df.CMOD_g - SF11['g'] * ebv; cr = df.CMOD_r - SF11['r'] * ebv; cz = df.CMOD_z - SF11['z'] * ebv
    df['COLOR'] = np.where(z < 0.5, cg - cr, cr - cz)
    df['EW'] = gc.line_ew(z, df.HALPHA_EW, df.OII_SUM_EW)
    ix = gc.indices(z, df.MABS, df.COLOR, df.EW, df.DN4000)
    for k in gc.KEY:
        df[k] = ix[k]
    df['VALID'] = ix['valid']
    code = gc.cell_code(ix)
    df['ACC_INDEX'] = np.array([cellmap.get(int(c), -1) if v else -1 for c, v in zip(code, ix['valid'])], dtype=np.int32)
    # Step 6
    s = np.full(n, np.nan); nok = np.zeros(n, bool)
    for zb in np.unique(ix['ZB']):
        r = ix['ZB'] == zb
        win = normwin.slice_window(zb * gc.DZ, (zb + 1) * gc.DZ, 'GALAXY', OBS_RANGE[survey])
        if win is not None:
            s[r], nok[r] = normwin.normalize(RWAVE, rflux[r], rivar[r], win)
    df['NORM_SCALE'], df['NORM_OK'] = s, nok
    use = nok & (df.ACC_INDEX.to_numpy() >= 0); nacc = int(use.sum())
    if use.any():
        fn = rflux[use] * s[use, None]; ivn = rivar[use] / s[use, None]**2; ci = df.ACC_INDEX.to_numpy()[use]
        for c in np.unique(ci):
            m = ci == c
            acc['wf'][c] += (ivn[m] * fn[m]).sum(0); acc['w'][c] += ivn[m].sum(0)
            acc['wf2'][c] += (ivn[m] * fn[m]**2).sum(0); acc['n'][c] += (ivn[m] > 0).sum(0); acc['ncell'][c] += m.sum()
    if store_dir is not None and nok.any():
        sd = os.path.join(store_dir, survey); os.makedirs(sd, exist_ok=True)
        keep = np.where(nok)[0]; starts, flx, ivv, lens = [], [], [], []
        for i in keep:
            cov = np.where(rivar[i] > 0)[0]
            if len(cov) == 0:
                starts.append(0); lens.append(0); continue
            a0, a1 = cov[0], cov[-1] + 1
            starts.append(a0); lens.append(a1 - a0)
            flx.append((rflux[i, a0:a1] * s[i]).astype(np.float16)); ivv.append(np.clip(rivar[i, a0:a1] / s[i]**2, 0, 65000).astype(np.float16))
        tmp = os.path.join(sd, base + '.tmp.npz')
        np.savez(tmp, FIBERID=df.FIBERID.to_numpy()[keep], ACC_INDEX=df.ACC_INDEX.to_numpy()[keep], START=np.array(starts, np.int32),
                 OFFSET=np.concatenate([[0], np.cumsum(lens)]).astype(np.int64), FLUX=np.concatenate(flx) if flx else np.zeros(0, np.float16),
                 IVAR=np.concatenate(ivv) if ivv else np.zeros(0, np.float16), RID_MIN=pb.RID_MIN)
        os.replace(tmp, os.path.join(sd, base + '.npz'))
    df['FILE'] = base
    return df, nacc, time.time() - t0


GLOBAL_DELTA = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prep', action='store_true'); ap.add_argument('--survey', choices=['dr8', 'dr17'])
    ap.add_argument('--shard', default='0/1'); ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--store', default=None); ap.add_argument('--nfiles', type=int, default=0)
    a = ap.parse_args()
    if a.prep:
        prep(); return
    t0 = time.time()
    cat = pd.read_parquet(CAT); cat = cat[cat.SURVEY == a.survey]
    o = pd.read_parquet(OFFS); o = o[(o.SURVEY == a.survey) & (o.N_STAR >= 3)]
    GLOBAL_DELTA[a.survey] = (float(o.DELTA_R.median()), float(o.DELTA_G.median()))
    offs = {(a.survey, int(p), int(m)): (r, g if np.isfinite(g) else GLOBAL_DELTA[a.survey][1]) for p, m, r, g in zip(o.PLATE, o.MJD, o.DELTA_R, o.DELTA_G)}
    cells, cellmap = pc.load_cells(); nc = len(cells)
    paths = sorted(sum([glob.glob(f'{d}/*/spPlate-*.fits') for d in PLATE_DIRS[a.survey]], []))
    if a.nfiles:
        paths = paths[::max(1, len(paths) // a.nfiles)][:a.nfiles]
    si, sn = (int(v) for v in a.shard.split('/')); paths = paths[si::sn]
    accp = os.path.join(a.outdir, 'acc', f'acc_{a.survey}_shard{si}of{sn}.npz'); os.makedirs(os.path.dirname(accp), exist_ok=True)
    if os.path.exists(accp):
        zz = np.load(accp); acc = {k: zz[k] for k in ('wf', 'w', 'wf2', 'n', 'ncell')}; done = set(zz['done'].tolist())
    else:
        acc = dict(wf=np.zeros((nc, NREST), F32), w=np.zeros((nc, NREST), F32), wf2=np.zeros((nc, NREST), F32),
                   n=np.zeros((nc, NREST), np.int32), ncell=np.zeros(nc, np.int64)); done = set()
    meta = dict(cells=pc.CELLTABLE, acc_min_n=pc.ACC_MIN_N, nrest=NREST, rid_min=pb.RID_MIN, dlog=DLOG, lines='vacuum', ncells=nc,
                survey=a.survey, delta_global=GLOBAL_DELTA[a.survey])
    print(f'{a.survey}: {len(cat):,} galaxies, {len(paths)} plates in this shard, global PSF-fiber offset r/g {GLOBAL_DELTA[a.survey]}, '
          f'{len(done)} done, setup {time.time()-t0:.0f}s', flush=True)
    status = os.path.expanduser(f'~/data/sdssdesi_status/step_gal_sdss_{a.survey}_shard{si}of{sn}.json')
    since = 0; ngal = 0; nacc = 0
    for k, p in enumerate(paths):
        if STOP:
            break
        key = os.path.basename(p)
        if key in done:
            continue
        try:
            tab, na, dt = process_plate(p, a.survey, cat, offs, cellmap, acc, a.store)
        except Exception as e:
            import traceback; traceback.print_exc(); print('ERROR', p, repr(e), flush=True); continue
        td = os.path.join(a.outdir, 'tables', a.survey); os.makedirs(td, exist_ok=True)
        if len(tab):
            out = os.path.join(td, key.replace('.fits', '.parquet')); tab.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
        done.add(key); since += 1; ngal += len(tab); nacc += na
        print(f'{a.survey} {k+1}/{len(paths)} {key} gal={len(tab)} stacked={na} {dt:.1f}s elapsed {time.time()-t0:.0f}s', flush=True)
        if since >= 50:
            pc.save_acc(accp, acc, done, meta); since = 0
        json.dump(dict(done=len(done), total=len(paths), status='running', note=f'{a.survey} galaxies {ngal:,}, stacked {nacc:,}',
                       started=t0), open(status, 'w'))
    pc.save_acc(accp, acc, done, meta)
    json.dump(dict(done=len(done), total=len(paths), status=('stopped' if STOP else 'done'), note=f'{a.survey} galaxies {ngal:,}, stacked {nacc:,}',
                   started=t0), open(status, 'w'))
    print(f'{"stopped" if STOP else "done"}: plates {len(done)}, galaxies {ngal}, stacked {nacc}, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
