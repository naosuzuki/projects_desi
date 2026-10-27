#!/usr/bin/env python3
"""Pass C: re-process every DESI galaxy spectrum and accumulate the cell stacks (ms.tex Sec. 3.9 Step 6).

Per galaxy (SPECTYPE=GALAXY, ZWARN=0, COADD_FIBERSTATUS=0, ZCAT_PRIMARY, DELTACHI2>25, 0<z<1.6):
  Step 1  fiber-photometry scale (as pass B; if the +-0.5 mag check fails the scale is 1 and CALC=0 - such
          spectra are kept, because every spectrum is renormalized below and K is scale-free)
  Step 2  Milky Way extinction (F99, R_V=3.1, 0.86 x SFD)
  Step 3  rest frame on the log10 grid (dlog=1e-4, 1400-9900 A), F_lambda*(1+z)
  Step 4  K-corrections from the galaxy's own spectrum; band-switched absolute magnitude from the whole-galaxy flux
  Step 5  lines re-measured with the VACUUM table; Dn4000; cell indices (gridcells.py)
  Step 6  normalization to unit mean F_lambda in the line-free window of the z slice (normwin.py);
          accumulation per cell: sum(w f), sum(w), sum(w f^2), count, with w = normalized ivar
Optional (--store DIR): the normalized rest-frame spectra (float16 flux and ivar, cropped to coverage) are
written per coadd, for median stacks and the EM-PCA.

Outputs
  ~/data/sdssdesi_passC/tables/<program>/hpx<N>.parquet   per-galaxy table (atomic; also the record of done files)
  ~/data/sdssdesi_passC/acc/acc_shard<i>of<n>.npz         accumulator checkpoint (atomic) + list of files included
Resumable: a file is skipped when it is in the shard's accumulator checkpoint; files whose table exists but which
are not in the checkpoint are processed again (their table is rewritten identically).

Usage: OPENBLAS_NUM_THREADS=1 python3 08_passC_stack.py --all --shard i/4 [--store /Volumes/exdisk3/sdssdesi_restframe]
       python3 08_passC_stack.py --files <coadd paths> --outdir <test dir>
"""
import os, sys, time, json, argparse, signal, warnings
import numpy as np
import pandas as pd
import fitsio

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')
import normwin
import gridcells as gc

assert getattr(pb, 'LINES_FRAME', 'air') == 'vacuum', 'pass C requires the vacuum LINES table'
OUTDIR = os.path.expanduser('~/data/sdssdesi_passC')
CELLTABLE = os.path.expanduser('~/data/sdssdesi_products/galaxy_cell_table.csv')
ACC_MIN_N = 100                 # accumulate cells that had >= 100 members in the pass-B assignment
CHECKPOINT_EVERY = 100          # files
NREST, RWAVE, DLOG = pb.NREST, pb.RWAVE, pb.DLOG
F32 = np.float32
STOP = False


def _term(sig, frm):
    global STOP
    STOP = True                 # finish the current file, checkpoint, exit


signal.signal(signal.SIGTERM, _term)
signal.signal(signal.SIGINT, _term)


def load_cells():
    t = pd.read_csv(CELLTABLE)
    t = t[t.N >= ACC_MIN_N].reset_index(drop=True)
    ix = {k: t[k].to_numpy() for k in gc.KEY}
    codes = gc.cell_code(ix)
    return t, {int(c): i for i, c in enumerate(codes)}


def galaxies_of_file(path, zl, offsets, prog, hpx):
    s, e = offsets.get((0 if prog == 'dark' else 1) * 10**6 + hpx, (0, 0))
    z = pd.DataFrame({'TARGETID': zl['TARGETID'][s:e], 'Z': zl['Z'][s:e].astype(float), 'ZWARN': zl['ZWARN'][s:e],
                      'SCODE': zl['SPECTYPE_CODE'][s:e], 'DELTACHI2': zl['DELTACHI2'][s:e].astype(float),
                      'ZCAT_PRIMARY': zl['ZCAT_PRIMARY'][s:e]})
    return z[(z.SCODE == 0) & (z.ZWARN == 0) & z.ZCAT_PRIMARY & (z.DELTACHI2 > 25) & (z.Z > 0.001) & (z.Z < gc.ZMAX)]


def process_chunk(df, w, flux, ivar, mask):
    """Steps 1-5 for galaxy rows; returns (df with measurements, rflux, rivar, rinrange) on the rest grid."""
    n = len(df)
    good = (ivar > 0) & (mask == 0)
    with np.errstate(all='ignore'):
        df['SNR_SDSSPIX'] = np.nanmedian(np.where(good, flux * np.sqrt(ivar), np.nan) * np.sqrt(w * np.log(10) * DLOG / 0.8).astype(F32), axis=1)
    # Step 1
    psf2fib = df['MEAN_PSF_TO_FIBER_SPECFLUX'].to_numpy(dtype=float)
    psf2fib = np.where(np.isfinite(psf2fib) & (psf2fib > 0), psf2fib, np.nan)
    scale0 = (np.nan_to_num(psf2fib, nan=1.0) * pb.FLUX_UNIT)[:, None]
    ffib = flux * scale0; ivfib = ivar / scale0**2
    south = df['PHOTSYS'].to_numpy() == 'S'
    names = {True: ('decam2014-g', 'decam2014-r'), False: ('BASS-g', 'BASS-r')}
    syn = {b: np.full(n, np.nan) for b in 'gr'}; err = {b: np.full(n, np.nan) for b in 'gr'}
    for sf in (True, False):
        rows = south == sf
        if rows.any():
            for j, b in enumerate('gr'):
                syn[b][rows], err[b][rows] = pb.ab_mags(pb.filt(names[sf][j]), w, ffib[rows], good[rows], ivar2d=ivfib[rows])
    with np.errstate(all='ignore'):
        dg = 22.5 - 2.5 * np.log10(df.FIBERFLUX_G.to_numpy(float)) - syn['g']
        dr = 22.5 - 2.5 * np.log10(df.FIBERFLUX_R.to_numpy(float)) - syn['r']
        Rg, Rr = 10**(-0.4 * dg), 10**(-0.4 * dr)
        wg = np.where(np.isfinite(Rg), 1 / np.where(np.isfinite(err['g']) & (err['g'] > 0), err['g'], 0.3)**2, 0)
        wr = np.where(np.isfinite(Rr), 1 / np.where(np.isfinite(err['r']) & (err['r'] > 0), err['r'], 0.3)**2, 0)
        a = (np.nan_to_num(Rg) * wg + np.nan_to_num(Rr) * wr) / np.where(wg + wr > 0, wg + wr, np.nan)
    calok = np.isfinite(a) & (a > 0) & (np.abs(np.nan_to_num(dg)) < 0.5) & (np.abs(np.nan_to_num(dr)) < 0.5) & (np.isfinite(dg) | np.isfinite(dr))
    a = np.where(calok, a, 1.0)
    df['CAL_A'], df['CALC'] = a, calok
    fcal = ffib * a[:, None]; icar = np.where(good, ivfib / a[:, None]**2, 0.0)
    del ffib, ivfib
    # Step 2
    ebv = np.nan_to_num(df.EBV.to_numpy(float), nan=0.0)
    if pb.A1_F99 is None or len(pb.A1_F99) != len(w):
        pb.A1_F99 = pb.extinction.fitzpatrick99(w, 1.0, 3.1)
    A = (3.1 * 0.86 * ebv)[:, None] * pb.A1_F99[None, :]
    fder = fcal * 10**(0.4 * A); ider = icar * 10**(-0.8 * A)
    del fcal, icar, A
    m_obs_r = np.full(n, np.nan)
    for sf in (True, False):
        rows = south == sf
        if rows.any():
            m_obs_r[rows], _ = pb.ab_mags(pb.filt(names[sf][1]), w, fder[rows], good[rows])
    # Step 3
    rid_obs, rflux_o, rivar_o = pb.rebin(w, fder, ider, mask)
    del fder, ider
    z = df.Z.to_numpy(float)
    rflux = np.zeros((n, NREST), F32); rivar = np.zeros((n, NREST), F32); rin = np.zeros((n, NREST), bool)
    shift = np.round(np.log10(1 + z) / DLOG).astype(int)
    for i in range(n):
        r0 = rid_obs - shift[i] - pb.RID_MIN
        ok = (r0 >= 0) & (r0 < NREST)
        rflux[i, r0[ok]] = rflux_o[i, ok] * (1 + z[i]); rivar[i, r0[ok]] = rivar_o[i, ok] / (1 + z[i])**2; rin[i, r0[ok]] = True
    del rflux_o, rivar_o
    # Step 4
    DM = pb.COSMO.distmod(np.clip(z, 1e-4, None)).value
    with np.errstate(all='ignore'):
        m_tot = 22.5 - 2.5 * np.log10(df.FLUX_R.to_numpy(float)) - pb.LS_EXT['R'] * ebv
        for b, nm in {'R': 'decam2014-r', 'G': 'decam2014-g', 'U': 'sdss2010-u', 'NUV': 'nuv'}.items():
            mr, _ = pb.ab_mags(pb.filt(nm), RWAVE, rflux, rivar > 0, inrange2d=rin)
            df[f'MABS_{b}'] = m_tot - DM - (m_obs_r - mr)
    # Step 5
    for k, v in pb.fit_lines(RWAVE, rflux, rivar).items():
        df[k] = v
    df['DN4000'], df['DN4000_ERR'] = pb.dn4000(RWAVE, rflux, rivar)
    return df, rflux, rivar, rin


FM_COLS = ['TARGETID', 'EBV', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R',
           'FIBERFLUX_Z', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_FIBERSTATUS']


def galaxy_rows(fm, path, prog, hpx, zl, offsets):
    """Galaxy rows of a coadd: DataFrame with the fibermap columns, the redshift columns and ROW (index into the arrays)."""
    gal = galaxies_of_file(path, zl, offsets, prog, hpx)
    fmd = pd.DataFrame({c: (np.char.strip(fm[c].astype(str)) if fm[c].dtype.kind in 'SU' else fm[c].astype(fm[c].dtype.newbyteorder('='))) for c in fm.dtype.names})
    fmd['ROW'] = np.arange(len(fmd))
    df = fmd.merge(gal, on='TARGETID', how='inner')
    return df[df.COADD_FIBERSTATUS == 0].sort_values('ROW').reset_index(drop=True)


def process_file(path, prog, hpx, zl, offsets, cellmap, acc, store_dir):
    """Standalone pass C: read the coadd, then process_arrays."""
    t0 = time.time()
    f = fitsio.FITS(path)
    fm = f['FIBERMAP'].read(columns=FM_COLS)
    df = galaxy_rows(fm, path, prog, hpx, zl, offsets)
    if len(df) == 0:
        f.close(); return pd.DataFrame(), 0, 0, time.time() - t0
    w, flux, ivar, mask, _ = pb.combine_cameras(f)
    f.close()
    out, ng, na = process_arrays(df, w, flux, ivar, mask, prog, hpx, cellmap, acc, store_dir)
    return out, ng, na, time.time() - t0


def process_arrays(df, w, flux, ivar, mask, prog, hpx, cellmap, acc, store_dir):
    """Steps 1-6 for the galaxy rows `df` (from galaxy_rows) given the combined-camera arrays of the whole coadd.
    Used by the standalone pass C and by the shared DESI walk (galaxy_hook.py)."""
    rows = df.ROW.to_numpy()
    flux, ivar, mask = flux[rows], ivar[rows], mask[rows]
    parts = []; n_acc = 0; store = []
    for i0 in range(0, len(df), pb.CHUNK):
        sl = slice(i0, min(i0 + pb.CHUNK, len(df)))
        d, rflux, rivar, rin = process_chunk(df.iloc[sl].copy(), w, flux[sl], ivar[sl], mask[sl])
        z = d.Z.to_numpy()
        d['MABS'] = gc.pick_mabs(z, d.MABS_R, d.MABS_G, d.MABS_U, d.MABS_NUV)
        d['MBAND'] = gc.mband(z)
        d['COLOR'] = gc.color(z, d.FLUX_G.to_numpy(float), d.FLUX_R.to_numpy(float), d.FLUX_Z.to_numpy(float), d.EBV)
        d['EW'] = gc.line_ew(z, d.HALPHA_EW, d.OII_SUM_EW)
        ix = gc.indices(z, d.MABS, d.COLOR, d.EW, d.DN4000)
        for k in gc.KEY:
            d[k] = ix[k]
        d['VALID'] = ix['valid']
        code = gc.cell_code(ix)
        d['ACC_INDEX'] = np.array([cellmap.get(int(c), -1) if v else -1 for c, v in zip(code, ix['valid'])], dtype=np.int32)
        # Step 6: normalization in the window of the z slice
        s = np.full(len(d), np.nan); nok = np.zeros(len(d), bool)
        for zb in np.unique(ix['ZB']):
            r = ix['ZB'] == zb
            win = normwin.slice_window(zb * gc.DZ, (zb + 1) * gc.DZ, 'GALAXY')
            if win is None:
                continue
            s[r], nok[r] = normwin.normalize(RWAVE, rflux[r], rivar[r], win)
        d['NORM_SCALE'], d['NORM_OK'] = s, nok
        use = nok & (d.ACC_INDEX.to_numpy() >= 0)
        if use.any():
            fn = rflux[use] * s[use, None]; iv = rivar[use] / s[use, None]**2
            ci = d.ACC_INDEX.to_numpy()[use]
            for c in np.unique(ci):
                m = ci == c
                acc['wf'][c] += (iv[m] * fn[m]).sum(axis=0); acc['w'][c] += iv[m].sum(axis=0)
                acc['wf2'][c] += (iv[m] * fn[m]**2).sum(axis=0); acc['n'][c] += (iv[m] > 0).sum(axis=0)
                acc['ncell'][c] += m.sum()
            n_acc += int(use.sum())
        if store_dir is not None:
            keep = nok
            for i in np.where(keep)[0]:
                cov = np.where(rivar[i] > 0)[0]
                if len(cov) == 0:
                    continue
                a0, a1 = cov[0], cov[-1] + 1
                store.append((d.TARGETID.iloc[i], int(d.ACC_INDEX.iloc[i]), a0, (rflux[i, a0:a1] * s[i]).astype(np.float16),
                              np.clip(rivar[i, a0:a1] / s[i]**2, 0, 65000).astype(np.float16)))
        parts.append(d.drop(columns=['ROW']))
        del rflux, rivar, rin
    out = pd.concat(parts, ignore_index=True)
    if store_dir is not None and store:
        sd = os.path.join(store_dir, prog); os.makedirs(sd, exist_ok=True)
        lens = np.array([len(x[3]) for x in store]); offs = np.concatenate([[0], np.cumsum(lens)])
        tmp = os.path.join(sd, f'hpx{hpx}.npz.tmp.npz')
        np.savez(tmp, TARGETID=np.array([x[0] for x in store]), ACC_INDEX=np.array([x[1] for x in store], np.int32),
                 START=np.array([x[2] for x in store], np.int32), OFFSET=offs.astype(np.int64),
                 FLUX=np.concatenate([x[3] for x in store]), IVAR=np.concatenate([x[4] for x in store]), RID_MIN=pb.RID_MIN)
        os.replace(tmp, os.path.join(sd, f'hpx{hpx}.npz'))
    return out, len(out), n_acc


def save_acc(path, acc, done, meta):
    tmp = path + '.tmp.npz'
    np.savez(tmp, wf=acc['wf'], w=acc['w'], wf2=acc['wf2'], n=acc['n'], ncell=acc['ncell'],
             done=np.array(sorted(done)), meta=json.dumps(meta))
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true'); ap.add_argument('--files', nargs='*')
    ap.add_argument('--shard', default='0/1'); ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--store', default=None, help='directory for the normalized rest-frame spectra (optional)')
    a = ap.parse_args()
    if a.files:
        paths = a.files
    else:
        paths = [l.split()[1] for l in open(pb.MANIFEST)]
    si, sn = (int(v) for v in a.shard.split('/'))
    paths = paths[si::sn]
    t0 = time.time()
    zl = dict(np.load(pb.ZLOOKUP_COMPACT))
    offsets = {int(k): (int(s), int(e)) for k, s, e in zip(zl['OFF_KEY'], zl['OFF_START'], zl['OFF_STOP'])}
    cells, cellmap = load_cells()
    nc = len(cells)
    accp = os.path.join(a.outdir, 'acc', f'acc_shard{si}of{sn}.npz')
    os.makedirs(os.path.dirname(accp), exist_ok=True)
    if os.path.exists(accp):
        z = np.load(accp, allow_pickle=False)
        acc = {k: z[k] for k in ('wf', 'w', 'wf2', 'n', 'ncell')}; done = set(z['done'].tolist())
        print(f'resumed accumulator: {len(done)} files', flush=True)
    else:
        acc = dict(wf=np.zeros((nc, NREST), F32), w=np.zeros((nc, NREST), F32), wf2=np.zeros((nc, NREST), F32),
                   n=np.zeros((nc, NREST), np.int32), ncell=np.zeros(nc, np.int64)); done = set()
    meta = dict(cells=CELLTABLE, acc_min_n=ACC_MIN_N, nrest=NREST, rid_min=pb.RID_MIN, dlog=DLOG, lines='vacuum', ncells=nc)
    print(f'cells accumulated: {nc}, accumulator {sum(v.nbytes for v in acc.values())/1e9:.2f} GB, setup {time.time()-t0:.0f}s', flush=True)
    since = 0; ngal = 0; nacc = 0
    for k, p in enumerate(paths):
        if STOP:
            break
        prog = 'dark' if '/dark/' in p else 'bright'
        hpx = int(os.path.basename(p).split('-')[-1].replace('.fits', ''))
        key = f'{prog}/{hpx}'
        if key in done:
            continue
        waited = 0
        while p.startswith('/Volumes/exdisk2/') and not os.path.ismount('/Volumes/exdisk2'):
            if waited == 0:
                print(f'DISK UNMOUNTED /Volumes/exdisk2 at {time.strftime("%H:%M:%S")} - waiting', flush=True)
            if STOP:
                break
            time.sleep(30); waited += 30
        lp = pb.local_path(p)
        if lp is None:
            print('MISSING', p, flush=True); continue
        try:
            tab, ng, na, dt = process_file(lp, prog, hpx, zl, offsets, cellmap, acc, a.store)
        except Exception as e:
            import traceback; traceback.print_exc(); print('ERROR', p, repr(e), flush=True); continue
        td = os.path.join(a.outdir, 'tables', prog); os.makedirs(td, exist_ok=True)
        if len(tab):
            tab.to_parquet(os.path.join(td, f'hpx{hpx}.parquet.tmp')); os.replace(os.path.join(td, f'hpx{hpx}.parquet.tmp'), os.path.join(td, f'hpx{hpx}.parquet'))
        done.add(key); since += 1; ngal += ng; nacc += na
        print(f'{k+1}/{len(paths)} {os.path.basename(lp)} gal={ng} stacked={na} {dt:.1f}s elapsed {time.time()-t0:.0f}s', flush=True)
        if since >= CHECKPOINT_EVERY:
            save_acc(accp, acc, done, meta); since = 0
    save_acc(accp, acc, done, meta)
    print(f'{"stopped" if STOP else "done"}: files {len(done)}, galaxies this run {ngal}, stacked {nacc}, {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
