#!/usr/bin/env python3
"""Pass S: stellar per-spectrum measurements over the DESI DR2 loa healpix coadds (STAR rows only).

For the STAR rows of each coadd (SPECTYPE == STAR and ZCAT_PRIMARY from the pass-B compact lookup; ZWARN kept as a
column) this script
  1. runs the pass-B chain (02_passB_measure.measure_chunk: Step 1 fiber-flux scale CAL_*, S/N, synthetic LS g/r/z,
     F99 extinction, rest frame, Dn4000, galaxy lines) -> the same columns as pass B (stars are the Step-1 validators);
  2. attaches the stellar catalog columns from ~/data/sdssdesi_cache/zall_stars_main.parquet (00_build_star_lookup.py):
     SUBTYPE, REF_ID (Gaia DR3 source_id), GAIA_PHOT_G/BP/RP_MEAN_MAG, PARALLAX, PMRA, PMDEC, MWS_TARGET, SCND_TARGET;
  3. runs starlib.measure_stars: synthetic Gaia G/BP/RP vs catalog, RV frame, spectral indices, band ratios, colors.
Output: one Parquet per coadd in OUTDIR/<program>/hpx<N>.parquet (resumable; empty Parquet when the file has no stars).
Design: docs/HANDOVER_stars.md. Production: call from the pass-Q DESI walk so each coadd is read once (or run standalone
after pass B with the same recipe as pass B: two workers, OPENBLAS_NUM_THREADS=1, nice 10, detached).

Usage (same CLI as pass B / pass Q):
  python3 02_passS_measure.py --files /Volumes/exdisk2/.../coadd-main-bright-0.fits --outdir ~/data/sdssdesi_passS_test
  python3 02_passS_measure.py --nfiles 5 --seed 1
  python3 02_passS_measure.py --all --shard 0/2
"""
import os, sys, time, argparse, warnings, gc
import numpy as np, pandas as pd, fitsio

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
passB = starlib.passB

OUTDIR = os.path.expanduser('~/data/sdssdesi_passS')
STAR_LOOKUP = os.path.expanduser('~/data/sdssdesi_cache/zall_stars_main.parquet')
STAR_COLS = ['TARGETID', 'PROGRAM', 'HEALPIX', 'SUBTYPE', 'REF_ID', 'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG',
             'GAIA_PHOT_RP_MEAN_MAG', 'PARALLAX', 'PMRA', 'PMDEC', 'MWS_TARGET', 'SCND_TARGET', 'COADD_EXPTIME', 'COADD_NUMEXP', 'MEAN_MJD']
FM_COLS = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'EBV', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R',
           'FIBERFLUX_Z', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_FIBERSTATUS']


def load_star_lookup():
    t0 = time.time()
    st = pd.read_parquet(STAR_LOOKUP, columns=STAR_COLS)
    st = st.set_index(['PROGRAM', 'HEALPIX']).sort_index()
    print(f'star lookup: {len(st):,} rows, {time.time()-t0:.0f}s', flush=True)
    return st


def measure_star_chunk(df, w, flux, ivar, mask):
    df = passB.measure_chunk(df, w, flux, ivar, mask)
    return starlib.measure_stars(df, w, flux, ivar, mask)


def process_file(path, zl, stars, program):
    t0 = time.time()
    f = fitsio.FITS(path)
    fm = f['FIBERMAP'].read(columns=FM_COLS)
    tid = fm['TARGETID']
    st = zl[(zl.SPECTYPE == 'STAR') & zl.ZCAT_PRIMARY]
    rows = np.where(np.isin(tid, st.TARGETID.to_numpy()))[0]
    if len(rows) == 0:
        f.close()
        return pd.DataFrame(), time.time() - t0
    w, flux, ivar, mask, percam = passB.combine_cameras(f)
    f.close()
    flux, ivar, mask = flux[rows], ivar[rows], mask[rows]
    percam = {c: v[rows] for c, v in percam.items()}
    fm = fm[rows]
    df = pd.DataFrame({c: (np.char.strip(fm[c].astype(str)) if fm[c].dtype.kind in 'SU' else fm[c].astype(fm[c].dtype.newbyteorder('='))) for c in fm.dtype.names})
    df['PROGRAM'] = program
    df = df.merge(st, on=['TARGETID', 'PROGRAM'], how='left')
    if stars is not None and len(stars):
        df = df.merge(stars.drop(columns=['PROGRAM', 'HEALPIX'], errors='ignore'), on='TARGETID', how='left')
    assert len(df) == len(rows)
    df['SNR_B'], df['SNR_R'], df['SNR_Z'] = percam['B'], percam['R'], percam['Z']
    t_read = time.time() - t0
    parts = []
    for i0 in range(0, len(df), passB.CHUNK):
        sl = slice(i0, min(i0 + passB.CHUNK, len(df)))
        parts.append(measure_star_chunk(df.iloc[sl].copy(), w, flux[sl], ivar[sl], mask[sl]))
        gc.collect()
    out = pd.concat(parts, ignore_index=True)
    out['FILE'] = os.path.basename(path)
    out['T_READ'] = t_read
    out['PROC_SEC'] = time.time() - t0
    return out, t_read


def stars_for(lookup, prog, hpx):
    try:
        s = lookup.loc[(prog, int(hpx))]
    except KeyError:
        return None
    s = s.reset_index()
    return s.drop_duplicates('TARGETID')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default=passB.MANIFEST)
    ap.add_argument('--nfiles', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--files', nargs='*')
    ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--shard', default='0/1')
    a = ap.parse_args()
    if a.files:
        paths = a.files
    else:
        paths = [l.split()[1] for l in open(a.manifest)]
        if not a.all:
            rng = np.random.default_rng(a.seed); paths = [paths[i] for i in rng.choice(len(paths), a.nfiles, replace=False)]
    si, sn = (int(v) for v in a.shard.split('/'))
    paths = paths[si::sn]
    t0 = time.time()
    zl = dict(np.load(passB.ZLOOKUP_COMPACT))
    offsets = {int(k): (int(s), int(e)) for k, s, e in zip(zl['OFF_KEY'], zl['OFF_START'], zl['OFF_STOP'])}
    print(f'zlookup loaded: {len(zl["TARGETID"]):,} rows, {time.time()-t0:.0f}s', flush=True)
    lookup = load_star_lookup()
    STYPE = np.array(['GALAXY', 'QSO', 'STAR', 'OTHER'])
    done = 0; missing = []; nst = 0
    for k, p in enumerate(paths):
        prog = 'dark' if '/dark/' in p else 'bright'
        hpx = os.path.basename(p).split('-')[-1].replace('.fits', '')
        out = os.path.join(a.outdir, prog, f'hpx{hpx}.parquet')
        if os.path.exists(out) and not a.overwrite:
            continue
        lp = passB.local_path(p)
        if lp is None:
            print('MISSING', p, flush=True); missing.append(p); continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        try:
            s, e = offsets.get((0 if prog == 'dark' else 1) * 10**6 + int(hpx), (0, 0))
            zsub = pd.DataFrame({'TARGETID': zl['TARGETID'][s:e], 'Z': zl['Z'][s:e].astype(float), 'ZERR': zl['ZERR'][s:e].astype(float),
                                 'ZWARN': zl['ZWARN'][s:e], 'SPECTYPE': STYPE[zl['SPECTYPE_CODE'][s:e]], 'DELTACHI2': zl['DELTACHI2'][s:e].astype(float),
                                 'ZCAT_PRIMARY': zl['ZCAT_PRIMARY'][s:e]})
            zsub['PROGRAM'] = prog
            df, t_read = process_file(lp, zsub, stars_for(lookup, prog, hpx), prog)
            df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
            done += 1; nst += len(df)
            print(f'{k+1}/{len(paths)} {os.path.basename(lp)} nstar={len(df)} read {t_read:.1f}s total {time.time()-t0:.0f}s', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print('ERROR', p, repr(e), flush=True)
    print(f'done: {done} files, {nst} stars, missing {len(missing)}, {time.time()-t0:.0f}s')
    if missing:
        with open(os.path.join(a.outdir, f'missing_shard{si}of{sn}.txt'), 'w') as fh:
            fh.write('\n'.join(missing) + '\n')


if __name__ == '__main__':
    main()
