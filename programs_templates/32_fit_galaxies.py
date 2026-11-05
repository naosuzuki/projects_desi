#!/usr/bin/env python3
"""Fit a galaxy basis to every DESI galaxy in the v2 rest-frame stores (Nao 2026-11-04): the validation table and the
coefficient-based selection of non-emission-line galaxies.
Per galaxy (store /Volumes/exdisk3/sdssdesi_restframe_v2/<prog>/hpxN.npz, rest grid RWAVE, Z from the v3t pass table):
weighted least squares of all vectors of the basis (reference-window form, EIGVEC) on the basis pixels (lines included),
CHI2 and NPIX of that fit, CHI2_E0 of the fit with e0 alone (the passive composite for the p1 basis), the median S/N per pixel,
and the coefficients C0..C10. Written per store file to <out>/<prog>/hpxN.parquet (skipped when present: resumable).
Usage: python3 32_fit_galaxies.py --basis FITS --out DIR [--workers 4] [--progs bright dark]
Status: ~/data/sdssdesi_status/step_gal_fit_<tag>.json
"""
import os, sys, glob, json, time, argparse, importlib
import numpy as np
import pandas as pd
import fitsio
from multiprocessing import Pool

STORE = '/Volumes/exdisk3/sdssdesi_restframe_v2'
TABLES = os.path.expanduser('~/data/sdssdesi_passC_v3t/tables')
STATUS = os.path.expanduser('~/data/sdssdesi_status')
G = {}


def init(basis):
    pb = importlib.import_module('02_passB_measure')
    F = fitsio.FITS(basis); P = F['EIGVEC'].read().astype(np.float64); w = F['WAVE'].read(); ok = F['PIXMASK'].read().astype(bool); F.close()
    assert len(w) == pb.NREST and np.allclose(w, pb.RWAVE, rtol=1e-6)
    G['P'] = P; G['ok'] = ok; G['pb'] = pb


def fit_file(args):
    prog, fn, out = args
    dst = os.path.join(out, prog, os.path.basename(fn).replace('.npz', '.parquet'))
    if os.path.exists(dst):
        return os.path.basename(fn), 0, 0.0
    t0 = time.time(); pb = G['pb']; P = G['P']; ok = G['ok']; nv = len(P)
    zz = np.load(fn); ids, st, off, rid0 = zz['TARGETID'], zz['START'], zz['OFFSET'], int(zz['RID_MIN'])
    n = len(ids); full = np.zeros((n, pb.NREST), np.float32); fiv = np.zeros_like(full)
    FL, IV = zz['FLUX'], zz['IVAR']
    for i in range(n):
        j0 = rid0 + int(st[i]) - pb.RID_MIN; fl = FL[off[i]:off[i + 1]]; full[i, j0:j0 + len(fl)] = fl; fiv[i, j0:j0 + len(fl)] = IV[off[i]:off[i + 1]]
    try:
        tab = pd.read_parquet(os.path.join(TABLES, prog, os.path.basename(fn).replace('.npz', '.parquet')), columns=['TARGETID', 'Z']).set_index('TARGETID')
        z = tab.Z.reindex(ids).to_numpy(float)
    except Exception:
        z = np.full(n, np.nan)
    X = full.astype(np.float64); W = np.where(ok[None, :] & np.isfinite(X) & (fiv > 0), fiv.astype(np.float64), 0.0); X = np.where(W > 0, X, 0.0)
    A = np.einsum('kp,np,lp->nkl', P, W, P); tr = np.einsum('nkk->n', A)
    A += (1e-9 * tr + (tr <= 0))[:, None, None] * np.eye(nv)[None]   # empty galaxies (no usable pixel) get A = I, C = 0
    b = np.einsum('kp,np->nk', P, W * X)
    try:
        C = np.linalg.solve(A, b[..., None])[..., 0]
    except np.linalg.LinAlgError:   # a galaxy with fewer pixels than vectors: least squares row by row
        C = np.stack([np.linalg.lstsq(A[i], b[i], rcond=None)[0] for i in range(n)])
    chi2 = (W * (X - C @ P)**2).sum(1); npix = (W > 0).sum(1)
    c0 = b[:, 0] / np.maximum(A[:, 0, 0], 1e-30); chi2_0 = (W * (X - c0[:, None] * P[0][None, :])**2).sum(1)
    with np.errstate(all='ignore'):
        snr = np.array([np.median((X[i] * np.sqrt(W[i]))[W[i] > 0]) if npix[i] > 0 else np.nan for i in range(n)])
    d = pd.DataFrame({'TARGETID': ids, 'Z': z, 'NPIX': npix.astype(np.int32), 'SNR': snr.astype(np.float32),
                      'CHI2': chi2.astype(np.float32), 'CHI2_E0': chi2_0.astype(np.float32)})
    for k in range(nv):
        d[f'C{k}'] = C[:, k].astype(np.float32)
    os.makedirs(os.path.dirname(dst), exist_ok=True); d.to_parquet(dst + '.tmp'); os.replace(dst + '.tmp', dst)
    return os.path.basename(fn), n, time.time() - t0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--basis', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--tag', default='p1'); ap.add_argument('--workers', type=int, default=4); ap.add_argument('--progs', nargs='+', default=['bright', 'dark'])
    a = ap.parse_args(); t0 = time.time()
    jobs = [(prog, fn, a.out) for prog in a.progs for fn in sorted(glob.glob(os.path.join(STORE, prog, 'hpx*.npz')))]
    sj = os.path.join(STATUS, f'step_gal_fit_{a.tag}.json'); tot = len(jobs)
    done = sum(os.path.exists(os.path.join(a.out, p, os.path.basename(f).replace('.npz', '.parquet'))) for p, f, _ in jobs); ngal = 0
    def st(status, note):
        json.dump(dict(done=done, total=tot, status=status, started=t0, note=note), open(sj, 'w'))
    st('running', f'fit of {os.path.basename(a.basis)}; {done} files present at start')
    print(f'{tot} store files, {done} already done; basis {a.basis}', flush=True)
    with Pool(a.workers, initializer=init, initargs=(a.basis,)) as pool:
        for name, n, dt in pool.imap_unordered(fit_file, jobs, chunksize=1):
            if n:
                done += 1; ngal += n
                print(f'{done}/{tot} {name} {n} galaxies {dt:.1f}s elapsed {time.time()-t0:.0f}s', flush=True)
                if done % 20 == 0:
                    st('running', f'{ngal:,} galaxies fitted this run')
    st('done', f'{ngal:,} galaxies fitted this run, {time.time()-t0:.0f}s')
    print(f'done: {done} files, {ngal} galaxies, {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
