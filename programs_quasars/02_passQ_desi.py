#!/usr/bin/env python3
"""Pass Q, DESI driver: quasar measurements over the DESI DR2 loa healpix coadds (QSO rows only).
Measurements: see qsolib.py (and docs/HANDOVER_quasars.md sec. 8.3/8.12/8.15). Output: one Parquet per coadd in
OUTDIR/<program>/hpx<N>.parquet (resumable; empty Parquet when a file has no QSO rows).
With --with-stars the same read also feeds the Stellar thread's STAR-row measurements through its hook
programs_stars/step2_3_pass_desi.write_star_file (output --star-outdir, default ~/data/sdssdesi_passS/<program>/hpx<N>.parquet):
one read per coadd for both tracks. With --expfibermap the EXP_FIBERMAP of every coadd (all rows) is dumped to
--expfm-outdir (default ~/data/sdssdesi_expfibermap/<program>/hpx<N>.parquet) for calibration 2 (all tracks).
Run after pass B has finished (the external disks are the bottleneck):
  python3 02_passQ_desi.py --nfiles 5 --seed 1
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --with-stars --expfibermap --with-galaxies --shard 0/4 > ~/data/sdssdesi_status/passQ_desi0.log 2>&1 &
"""
import os, sys, time, argparse, warnings, gc, signal
import numpy as np
import pandas as pd
import fitsio
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qsolib as Q

OUTDIR = os.path.expanduser('~/data/sdssdesi_passQ')

FM_COLS = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'EBV', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R',
           'FIBERFLUX_Z', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_FIBERSTATUS']
STARS = None          # programs_stars/step2_3_pass_desi module when --with-stars (Stellar thread's hook write_star_file)
STAR_OUTDIR = os.path.expanduser('~/data/sdssdesi_passS')
EXPFM_OUTDIR = os.path.expanduser('~/data/sdssdesi_expfibermap')
EXPFM_COLS = ['TARGETID', 'NIGHT', 'EXPID', 'MJD', 'TILEID', 'EXPTIME', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'FIBERSTATUS',
              'FIBER_X', 'FIBER_Y', 'DELTA_X', 'DELTA_Y', 'FIBER_RA', 'FIBER_DEC', 'PSF_TO_FIBER_SPECFLUX', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z']


def rows_df(fm, rows, program, zsub, percam):
    fm = fm[rows]
    df = pd.DataFrame({c: (np.char.strip(fm[c].astype(str)) if fm[c].dtype.kind in 'SU' else fm[c].astype(fm[c].dtype.newbyteorder('='))) for c in fm.dtype.names})
    df['PROGRAM'] = program
    df = df.merge(zsub, on=['TARGETID', 'PROGRAM'], how='left')
    assert len(df) == len(rows)
    df['SNR_B'], df['SNR_R'], df['SNR_Z'] = percam['B'][rows], percam['R'][rows], percam['Z'][rows]
    return df


def run_chunks(df, w, flux, ivar, mask, rows, fn):
    parts = []
    for i0 in range(0, len(df), Q.passB.CHUNK):
        sl = slice(i0, min(i0 + Q.passB.CHUNK, len(df))); r = rows[sl]
        parts.append(fn(df.iloc[sl].copy(), w, flux[r], ivar[r], mask[r]))
        gc.collect()
    return pd.concat(parts, ignore_index=True)


def dump_expfibermap(f, path, program, outdir, overwrite=False):
    """EXP_FIBERMAP of the coadd (all rows, all classes) -> outdir/<program>/hpx<N>.parquet: EXPID/PETAL/FIBER per
    exposure for the Stellar thread's calibration 2 (shared by all tracks). Returns rows written, -1 if skipped."""
    hpx = os.path.basename(path).split('-')[-1].replace('.fits', '')
    out = os.path.join(outdir, program, f'hpx{hpx}.parquet')
    if os.path.exists(out) and not overwrite:
        return -1
    os.makedirs(os.path.dirname(out), exist_ok=True)
    e = f['EXP_FIBERMAP'].read(columns=[c for c in EXPFM_COLS if c in f['EXP_FIBERMAP'].get_colnames()])
    df = pd.DataFrame({c: e[c].astype(e[c].dtype.newbyteorder('=')) for c in e.dtype.names})
    df['PROGRAM'] = program; df['HEALPIX'] = int(hpx)
    df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
    return len(df)


def process_file(path, zl, program, star_lookup=None, want_stars=False, want_q=True, expfm_dir=None, overwrite=False,
                 G=None, need_g=False, hpx=None, NH=None, need_n=False):
    """One read of the coadd: quasar measurements for the QSO rows (returned), the Stellar thread's hook for the STAR
    rows (written by S.write_star_file) and the EXP_FIBERMAP dump. Returns (qso_df, n_star, n_expfm, t_read)."""
    t0 = time.time()
    f = fitsio.FITS(path)
    n_expfm = dump_expfibermap(f, path, program, expfm_dir, overwrite) if expfm_dir else 0
    fm = f['FIBERMAP'].read(columns=FM_COLS)
    gfm = f['FIBERMAP'].read(columns=G_FM_COLS) if need_g else None      # Galaxy thread's hook (pass C on the same read)
    tid = fm['TARGETID']
    qso = zl[(zl.SPECTYPE == 'QSO') & zl.ZCAT_PRIMARY]
    qrows = np.where(np.isin(tid, qso.TARGETID.to_numpy()))[0] if want_q else np.array([], int)
    nstar_rows = int(np.isin(tid, zl[(zl.SPECTYPE == 'STAR') & zl.ZCAT_PRIMARY].TARGETID.to_numpy()).sum()) if want_stars else 0
    if len(qrows) == 0 and nstar_rows == 0 and not need_g and not need_n:
        f.close()
        return pd.DataFrame(), 0, n_expfm, time.time() - t0
    w, flux, ivar, mask, percam = Q.passB.combine_cameras(Q.NoMaskFITS(f))   # MASK HDUs skipped: masked pixels have IVAR = 0
    f.close()
    t_read = time.time() - t0
    if need_g:
        G.process(path, gfm, w, flux, ivar, mask, program, hpx)
    if need_n:
        NH.process(path, tid, w, flux, ivar)                              # v2 native-pixel store (quasars + stars)
    n_star = 0
    if want_stars:
        n_star = STARS.write_star_file(path, fm, w, flux, ivar, mask, percam, zl, star_lookup, program, outdir=STAR_OUTDIR_RUN, overwrite=overwrite)
    qdf = pd.DataFrame()
    if len(qrows):
        df = rows_df(fm, qrows, program, qso, percam)
        qdf = run_chunks(df, w, flux, ivar, mask, qrows, Q.measure_qso_chunk_desi)
        qdf['FILE'] = os.path.basename(path); qdf['T_READ'] = t_read; qdf['PROC_SEC'] = time.time() - t0
    return qdf, n_star, n_expfm, t_read


STAR_OUTDIR_RUN = STAR_OUTDIR
G_FM_COLS = None      # galaxy_hook.FM_COLS when --with-galaxies


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default=Q.passB.MANIFEST)
    ap.add_argument('--nfiles', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--files', nargs='*')
    ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--with-stars', action='store_true', help='also run the Stellar thread\'s measurements on the STAR rows (one read per coadd)')
    ap.add_argument('--star-outdir', default=STAR_OUTDIR)
    ap.add_argument('--expfibermap', action='store_true', help='dump EXP_FIBERMAP of every coadd (all rows) to --expfm-outdir')
    ap.add_argument('--expfm-outdir', default=EXPFM_OUTDIR)
    ap.add_argument('--with-galaxies', action='store_true', help="also run the Galaxy thread's pass-C hook (programs_templates/galaxy_hook.py) on the same read")
    ap.add_argument('--with-native', action='store_true', help='v2: store the native-pixel quasar + star spectra on exdisk3 (native_hook.py)')
    ap.add_argument('--gal-store', default='/Volumes/exdisk3/sdssdesi_restframe_v2', help='rest-frame spectrum store of the galaxy hook')
    ap.add_argument('--gal-outdir', default=None, help='override the galaxy hook state dir (testing only; default = pass C state)')
    a = ap.parse_args()
    global STARS, STAR_OUTDIR_RUN
    star_lookup = None; STAR_OUTDIR_RUN = os.path.expanduser(a.star_outdir)
    if a.with_stars:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'programs_stars'))
        import step2_3_pass_desi as STARS_mod
        STARS = STARS_mod
        star_lookup = STARS.load_star_lookup()
    if a.files:
        paths = a.files
    else:
        paths = [l.split()[1] for l in open(a.manifest)]
        if not a.all:
            rng = np.random.default_rng(a.seed); paths = [paths[i] for i in rng.choice(len(paths), a.nfiles, replace=False)]
    si, sn = (int(v) for v in a.shard.split('/'))
    paths = paths[si::sn]
    t0 = time.time()
    zl = dict(np.load(Q.passB.ZLOOKUP_COMPACT))
    offsets = {int(k): (int(s), int(e)) for k, s, e in zip(zl['OFF_KEY'], zl['OFF_START'], zl['OFF_STOP'])}
    print(f'zlookup loaded: {len(zl["TARGETID"]):,} rows, {time.time()-t0:.0f}s', flush=True)
    STYPE = np.array(['GALAXY', 'QSO', 'STAR', 'OTHER'])
    global G_FM_COLS
    G = None
    if a.with_galaxies:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'programs_templates'))
        import galaxy_hook as GH
        G_FM_COLS = GH.FM_COLS
        kw = dict(store=a.gal_store, zl=zl, offsets=offsets)
        if a.gal_outdir:
            kw['outdir'] = os.path.expanduser(a.gal_outdir)
        G = GH.GalaxyHook(si, sn, **kw)
        print(f'galaxy hook: shard {si}/{sn}, {len(G.done):,} coadds already in the accumulator', flush=True)
    NH = None
    if a.with_native:
        import native_hook; NH = native_hook.NativeHook()
        print(f'native hook: {len(NH.by_file):,} coadds with quasar/star rows', flush=True)
    # graceful stop: SIGTERM only sets a flag; the loop finishes the current coadd, then checkpoints and exits, so the
    # galaxy accumulator is never saved half-way through a coadd (a redone coadd would otherwise be counted twice)
    stop = {'flag': False}
    signal.signal(signal.SIGTERM, lambda sig, frm: stop.__setitem__('flag', True))
    done = 0; missing = []; nq = 0
    for k, p in enumerate(paths):
        if stop['flag']:
            print('SIGTERM: stopping after the last complete coadd', flush=True); break
        prog = 'dark' if '/dark/' in p else 'bright'
        hpx = os.path.basename(p).split('-')[-1].replace('.fits', '')
        out = os.path.join(a.outdir, prog, f'hpx{hpx}.parquet')
        sout = os.path.join(STAR_OUTDIR_RUN, prog, f'hpx{hpx}.parquet')
        eout = os.path.join(a.expfm_outdir, prog, f'hpx{hpx}.parquet')
        need_q = a.overwrite or not os.path.exists(out)
        need_s = a.with_stars and (a.overwrite or not os.path.exists(sout))
        need_e = a.expfibermap and (a.overwrite or not os.path.exists(eout))
        need_g = G is not None and G.needs(prog, int(hpx))
        need_n = NH is not None and os.path.basename(p) in NH.by_file and not os.path.exists(native_hook.out(p))
        if not (need_q or need_s or need_e or need_g or need_n):
            continue
        lp = Q.passB.local_path(p)
        if lp is None:
            print('MISSING', p, flush=True); missing.append(p); continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        try:
            s, e = offsets.get((0 if prog == 'dark' else 1) * 10**6 + int(hpx), (0, 0))
            zsub = pd.DataFrame({'TARGETID': zl['TARGETID'][s:e], 'Z': zl['Z'][s:e].astype(float), 'ZERR': zl['ZERR'][s:e].astype(float),
                                 'ZWARN': zl['ZWARN'][s:e], 'SPECTYPE': STYPE[zl['SPECTYPE_CODE'][s:e]], 'DELTACHI2': zl['DELTACHI2'][s:e].astype(float),
                                 'ZCAT_PRIMARY': zl['ZCAT_PRIMARY'][s:e]})
            zsub['PROGRAM'] = prog
            qdf, n_star, n_expfm, t_read = process_file(lp, zsub, prog, star_lookup=star_lookup, want_stars=need_s, want_q=need_q,
                                                        expfm_dir=(a.expfm_outdir if need_e else None), overwrite=a.overwrite,
                                                        G=G, need_g=need_g, hpx=int(hpx), NH=NH, need_n=need_n)
            if need_q:
                qdf.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
            done += 1; nq += len(qdf)
            print(f'{k+1}/{len(paths)} {os.path.basename(lp)} nqso={len(qdf)}' + (f' nstar={n_star}' if need_s else '') + (f' nexp={n_expfm}' if need_e else '')
                  + (' gal' if need_g else '') + (' nat' if need_n else '')
                  + f' read {t_read:.1f}s total {time.time()-t0:.0f}s', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print('ERROR', p, repr(e), flush=True)
    if G is not None:
        G.close(); print(f'galaxy hook checkpoint written: {len(G.done):,} coadds, {G.ngal:,} galaxies this run', flush=True)
    print(f'{"stopped" if stop["flag"] else "done"}: {done} files, {nq} quasars, missing {len(missing)}, {time.time()-t0:.0f}s')
    if missing:
        with open(os.path.join(a.outdir, f'missing_shard{si}of{sn}.txt'), 'w') as fh:
            fh.write('\n'.join(missing) + '\n')


if __name__ == '__main__':
    main()
