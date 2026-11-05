#!/usr/bin/env python3
"""Stack the non-emission-line galaxies selected by the coefficients of the basis fit (Nao 2026-11-04): a galaxy is passive when
its star-formation coefficient C1/C0 and its line coefficient C2/C0 of the p1 fit are consistent with zero, with ZWARN 0,
S/N >= --min-snr, a good fit (chi2/NPIX < --max-rchi2) and a magnitude bin (MB_5450, or MB_G where the 5450 window is not covered, z > 0.6). The selected galaxies are stacked per
(ZB, MB_5450) cell from the v2 rest-frame stores (flux already in the normalized, chained units: C0 ~ 1), inverse variance
weighted, and the cell stacks are combined into one composite with sqrt N weights (as the passive e0 of 31_passive_empca.py).
Output: <out>/stacks_select_<tag>.npz (WAVE, STACK [ncell, npix], IVAR, N, ZB, MB, E0, E0_COV) and a CSV of the cells.
Usage: python3 34_restack_select.py --fit ~/data/sdssdesi_galfit_p1 --out ~/data/sdssdesi_galfit_p1 --tag p1 [--r1-max 5 --r2-abs 5]
"""
import os, glob, json, time, argparse, importlib
import numpy as np, pandas as pd
from multiprocessing import Pool

STORE = '/Volumes/exdisk3/sdssdesi_restframe_v2'; TABLES = os.path.expanduser('~/data/sdssdesi_passC_v3t/tables')
STATUS = os.path.expanduser('~/data/sdssdesi_status'); ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
G = {}


def init(a):
    G['pb'] = importlib.import_module('02_passB_measure'); G['a'] = a


def one(args):
    prog, fn = args; a = G['a']; pb = G['pb']; base = os.path.basename(fn)
    f = pd.read_parquet(os.path.join(a.fit, prog, base.replace('.npz', '.parquet')))
    t = pd.read_parquet(os.path.join(TABLES, prog, base.replace('.npz', '.parquet')), columns=['TARGETID', 'ZWARN', 'ZB', 'MB_5450', 'MB_G']).drop_duplicates('TARGETID').set_index('TARGETID')
    f = f.join(t, on='TARGETID'); f['MB_5450'] = np.where(f.MB_5450 > -90, f.MB_5450, f.MB_G)   # z > 0.6: the 5450 window is not covered, bin by M_g
    sel = ((f.C1 / f.C0 < a.r1_max) & ((f.C2 / f.C0).abs() < a.r2_abs) & (f.ZWARN == 0) & (f.SNR >= a.min_snr) & (f.CHI2 / np.maximum(f.NPIX, 1) < a.max_rchi2)
           & (f.MB_5450 > -90) & (f.C0 > 0)).to_numpy()
    if a.norm == 'none':   # store units: unit mean in the galaxy's own window (5400-5500, or 3000-3100 at high z where C0 ~ 5 for passive galaxies)
        sel &= ((f.C0 > 0.5) & (f.C0 < 2.0)).to_numpy()
    c0 = f.C0.to_numpy(float)
    idx = np.where(sel)[0]
    if len(idx) == 0:
        return len(f), 0, {}
    zz = np.load(fn); ids, st, off, rid0 = zz['TARGETID'], zz['START'], zz['OFFSET'], int(zz['RID_MIN']); FL, IV = zz['FLUX'], zz['IVAR']
    assert np.array_equal(ids, f.TARGETID.to_numpy())
    out = {}
    for i in idx:
        key = (int(f.ZB.iat[i]), int(f.MB_5450.iat[i]))
        if key not in out:
            out[key] = [np.zeros(pb.NREST), np.zeros(pb.NREST), 0]
        j0 = rid0 + int(st[i]) - pb.RID_MIN; fl = FL[off[i]:off[i + 1]].astype(np.float64); iv = IV[off[i]:off[i + 1]].astype(np.float64)
        ok = np.isfinite(fl) & (iv > 0); iv = np.where(ok, iv, 0.0); fl = np.where(ok, fl, 0.0)
        if a.norm == 'c0':   # normalize by the fitted e0 amplitude instead of the window mean: every galaxy on the scale e0 = 1 at 5450 A
            fl = fl / c0[i]; iv = iv * c0[i]**2
        out[key][0][j0:j0 + len(fl)] += iv * fl; out[key][1][j0:j0 + len(fl)] += iv; out[key][2] += 1
    return len(f), len(idx), out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', required=True); ap.add_argument('--out', required=True); ap.add_argument('--tag', default='p1')
    ap.add_argument('--r1-max', type=float, default=5.0); ap.add_argument('--r2-abs', type=float, default=5.0); ap.add_argument('--min-snr', type=float, default=2.0)
    ap.add_argument('--max-rchi2', type=float, default=1.5); ap.add_argument('--workers', type=int, default=6); ap.add_argument('--min-cell', type=int, default=50)
    ap.add_argument('--norm', default='c0', choices=['c0', 'none'], help='c0: divide each galaxy by its fitted e0 coefficient (model-based normalization, no window bias, works at z > 0.6); none: store units, C0 in 0.5-2 required')
    a = ap.parse_args(); t0 = time.time(); pb = importlib.import_module('02_passB_measure')
    jobs = [(p, f) for p in ['bright', 'dark'] for f in sorted(glob.glob(os.path.join(STORE, p, 'hpx*.npz')))]
    sj = os.path.join(STATUS, f'step_gal_restack_{a.tag}.json')
    def st(done, status, note):
        json.dump(dict(done=done, total=len(jobs), status=status, started=t0, note=note), open(sj, 'w'))
    st(0, 'running', 'selecting and stacking'); acc = {}; ntot = nsel = 0
    with Pool(a.workers, initializer=init, initargs=(a,)) as pool:
        for i, (n, ns, out) in enumerate(pool.imap_unordered(one, jobs, chunksize=2)):
            ntot += n; nsel += ns
            for k, (sx, sw, nn) in out.items():
                if k not in acc:
                    acc[k] = [np.zeros(pb.NREST), np.zeros(pb.NREST), 0]
                acc[k][0] += sx; acc[k][1] += sw; acc[k][2] += nn
            if i % 200 == 0:
                st(i, 'running', f'{nsel:,} of {ntot:,} galaxies selected, {time.time()-t0:.0f}s'); print(i, nsel, ntot, flush=True)
    keys = sorted(k for k in acc if acc[k][2] >= a.min_cell)
    S = np.array([acc[k][0] / np.maximum(acc[k][1], 1e-30) for k in keys]); IV = np.array([acc[k][1] for k in keys]); N = np.array([acc[k][2] for k in keys])
    S[IV <= 0] = 0.0
    w = np.sqrt(N / np.median(N))[:, None] * (IV > 0); cov = (IV > 0).sum(0)
    e0 = np.where(cov >= 5, (w * S).sum(0) / np.maximum(w.sum(0), 1e-30), 0.0)
    os.makedirs(a.out, exist_ok=True); fn = os.path.join(a.out, f'stacks_select_{a.tag}.npz')
    np.savez_compressed(fn, WAVE=pb.RWAVE, STACK=S.astype(np.float32), IVAR=IV.astype(np.float32), N=N, ZB=np.array([k[0] for k in keys]), MB=np.array([k[1] for k in keys]),
                        E0=e0, E0_COV=cov, SEL=np.array([a.r1_max, a.r2_abs, a.min_snr, a.max_rchi2]))
    pd.DataFrame(dict(ZB=[k[0] for k in keys], MB=[k[1] for k in keys], N=N)).to_csv(os.path.join(ROOT, 'csv_files', f'gal_select_cells_{a.tag}.csv'), index=False)
    jd = np.where(cov >= 5)[0]
    print(f'{nsel:,} of {ntot:,} galaxies selected ({nsel/ntot:.1%}); {len(keys)} cells with >= {a.min_cell}; composite defined {pb.RWAVE[jd[0]]:.0f}-{pb.RWAVE[jd[-1]]:.0f} A; {time.time()-t0:.0f}s', flush=True)
    print('wrote', fn, flush=True); st(len(jobs), 'done', f'{nsel:,} of {ntot:,} galaxies, {len(keys)} cells')


if __name__ == '__main__':
    main()
