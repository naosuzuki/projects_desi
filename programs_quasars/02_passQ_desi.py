#!/usr/bin/env python3
"""Pass Q, DESI driver: quasar measurements over the DESI DR2 loa healpix coadds (QSO rows only).
Measurements: see qsolib.py (and docs/HANDOVER_quasars.md sec. 8.3/8.12/8.15). Output: one Parquet per coadd in
OUTDIR/<program>/hpx<N>.parquet (resumable; empty Parquet when a file has no QSO rows).
Run after pass B has finished (the external disks are the bottleneck):
  python3 02_passQ_desi.py --nfiles 5 --seed 1
  OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --shard 0/2 > ~/data/sdssdesi_status/passQ_desi0.log 2>&1 &
"""
import os, sys, time, argparse, warnings, gc
import numpy as np
import pandas as pd
import fitsio
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qsolib as Q

OUTDIR = os.path.expanduser('~/data/sdssdesi_passQ')

def process_file(path, zl, program):
    t0 = time.time()
    f = fitsio.FITS(path)
    fm = f['FIBERMAP'].read(columns=['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'EBV', 'PHOTSYS',
                                      'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
                                      'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_FIBERSTATUS'])
    tid = fm['TARGETID']
    qso = zl[(zl.SPECTYPE == 'QSO') & zl.ZCAT_PRIMARY]
    rows = np.where(np.isin(tid, qso.TARGETID.to_numpy()))[0]
    if len(rows) == 0:
        f.close()
        return pd.DataFrame(), time.time() - t0
    w, flux, ivar, mask, percam = Q.passB.combine_cameras(f)
    f.close()
    flux, ivar, mask = flux[rows], ivar[rows], mask[rows]
    percam = {c: v[rows] for c, v in percam.items()}
    fm = fm[rows]
    df = pd.DataFrame({c: (np.char.strip(fm[c].astype(str)) if fm[c].dtype.kind in 'SU' else fm[c].astype(fm[c].dtype.newbyteorder('='))) for c in fm.dtype.names})
    df['PROGRAM'] = program
    df = df.merge(qso, on=['TARGETID', 'PROGRAM'], how='left')
    assert len(df) == len(rows)
    df['SNR_B'], df['SNR_R'], df['SNR_Z'] = percam['B'], percam['R'], percam['Z']
    t_read = time.time() - t0
    parts = []
    for i0 in range(0, len(df), Q.passB.CHUNK):
        sl = slice(i0, min(i0 + Q.passB.CHUNK, len(df)))
        parts.append(Q.measure_qso_chunk_desi(df.iloc[sl].copy(), w, flux[sl], ivar[sl], mask[sl]))
        gc.collect()
    out = pd.concat(parts, ignore_index=True)
    out['FILE'] = os.path.basename(path)
    out['T_READ'] = t_read
    out['PROC_SEC'] = time.time() - t0
    return out, t_read


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
    zl = dict(np.load(Q.passB.ZLOOKUP_COMPACT))
    offsets = {int(k): (int(s), int(e)) for k, s, e in zip(zl['OFF_KEY'], zl['OFF_START'], zl['OFF_STOP'])}
    print(f'zlookup loaded: {len(zl["TARGETID"]):,} rows, {time.time()-t0:.0f}s', flush=True)
    STYPE = np.array(['GALAXY', 'QSO', 'STAR', 'OTHER'])
    done = 0; missing = []; nq = 0
    for k, p in enumerate(paths):
        prog = 'dark' if '/dark/' in p else 'bright'
        hpx = os.path.basename(p).split('-')[-1].replace('.fits', '')
        out = os.path.join(a.outdir, prog, f'hpx{hpx}.parquet')
        if os.path.exists(out) and not a.overwrite:
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
            df, t_read = process_file(lp, zsub, prog)
            df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
            done += 1; nq += len(df)
            print(f'{k+1}/{len(paths)} {os.path.basename(lp)} nqso={len(df)} read {t_read:.1f}s total {time.time()-t0:.0f}s', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print('ERROR', p, repr(e), flush=True)
    print(f'done: {done} files, {nq} quasars, missing {len(missing)}, {time.time()-t0:.0f}s')
    if missing:
        with open(os.path.join(a.outdir, f'missing_shard{si}of{sn}.txt'), 'w') as fh:
            fh.write('\n'.join(missing) + '\n')


if __name__ == '__main__':
    main()
