#!/usr/bin/env python3
"""Quasar eigenspectra: weighted EM-PCA (Bailey 2012) of the quasar stacks (05_stack.py --merge -> qso_stacks*.npz).
Same conventions as the galaxy thread (programs_templates/10_galaxy_empca.py), routine programs_stars/step7_empca.empca:
- inputs: the ivar-weighted mean stacks WMEAN (unit mean in the normalization window of their z slice), cells with
  at least --nmin quasars at their best-covered pixel;
- no mean subtraction (coverage depends on z): the mean shape is the first of nvec+1 vectors;
- each stack counts once: its weights 1/WERR^2 are rescaled to a median of 1 over its good pixels;
- a pixel enters the basis where >= --min-cover stacks cover it.
Output: ~/data/sdssdesi_products/quasar/qso_empca_<tag>.fits: EIGVEC [nvec+1, npix], COEFF [nstack, nvec+1], PIXMASK,
VARIANCE (weighted R^2 vs number of vectors), WAVE, CELLS (+ qso_empca_<tag>.png quick look).
Usage: python3 06_empca.py --stacks ~/data/sdssdesi_products/quasar/qso_stacks.npz --nvec 10 --niter 25 --tag v1
       python3 06_empca.py ... --zrange 0.8 2.0 --tag v1_z0820     (one basis per normalization window)
"""
import os, sys, time, json, argparse
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'programs_stars'))
from step7_empca import empca, explained          # weighted EM-PCA with missing data (stellar thread, stable interface)

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
STATUS_V1 = os.path.expanduser('~/data/sdssdesi_status/step_qso_empca.json')
STATUS_V2 = os.path.expanduser('~/data/sdssdesi_status/step_qso_empca_v2.json')   # own gauge row, so v2 never overwrites v1
TAG = ''


def status(done, total, st, note, t0):
    STATUS = STATUS_V2 if TAG.startswith('v2') else STATUS_V1
    json.dump(dict(done=done, total=total, status=st, note=note, started=t0), open(STATUS + '.tmp', 'w')); os.replace(STATUS + '.tmp', STATUS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', default=f'{PROD}/qso_stacks.npz')
    ap.add_argument('--nvec', type=int, default=10)
    ap.add_argument('--niter', type=int, default=25)
    ap.add_argument('--nmin', type=int, default=300, help='minimum quasars per stack (at its best-covered pixel)')
    ap.add_argument('--min-cover', type=int, default=20)
    ap.add_argument('--min-pix', type=int, default=500, help='minimum good pixels per stack')
    ap.add_argument('--tag', default='v1')
    ap.add_argument('--zrange', type=float, nargs=2, default=None, help='only stacks whose z slice centre is in [lo, hi): one basis per normalization window')
    a = ap.parse_args()
    global TAG; TAG = a.tag
    t0 = time.time()
    status(0, 3, 'running', f'{a.tag}: loading stacks', t0)
    z = np.load(a.stacks, allow_pickle=True)
    wave, wm, we, n, cells = z['WAVE'], z['WMEAN'], z['WERR'], z['N'], z['STACK_CELL']
    keep = n.max(axis=1) >= a.nmin
    if a.zrange:
        import importlib; ZE = np.asarray(importlib.import_module('04_cells').ZEDGES, float)
        zb = np.array([int(str(c).split('_')[0][1:]) for c in cells]); zc = 0.5 * (ZE[zb] + ZE[zb + 1])
        keep &= (zc >= a.zrange[0]) & (zc < a.zrange[1])
    X = np.where(np.isfinite(wm) & np.isfinite(we) & (we > 0), wm, 0.0)[keep].astype(np.float64)
    W = np.where(np.isfinite(wm) & np.isfinite(we) & (we > 0), 1.0 / np.where(we > 0, we, 1.0)**2, 0.0)[keep]
    good = (W > 0).sum(axis=1) >= a.min_pix
    X, W, cells_k = X[good], W[good], cells[keep][good]
    for i in range(len(W)):
        W[i] /= np.median(W[i][W[i] > 0])
    print(f'{len(X)} stacks (n >= {a.nmin}, >= {a.min_pix} good pixels) of {len(cells)}; {X.shape[1]} pixels', flush=True)
    status(1, 3, 'running', f'{a.tag}: EM-PCA on {len(X)} stacks', t0)
    nv = a.nvec + 1
    P, C, pixok = empca(X, W, nv, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    R2 = explained(X, W, P)
    print('weighted R^2 vs number of vectors:', ' '.join(f'{r:.4f}' for r in R2), flush=True)
    out = f'{PROD}/qso_empca_{a.tag}.fits'
    with fitsio.FITS(out + '.tmp', 'rw', clobber=True) as f:
        f.write(P.astype(np.float32), extname='EIGVEC')
        f.write(C.astype(np.float32), extname='COEFF')
        f.write(pixok.astype(np.int8), extname='PIXMASK')
        f.write(np.asarray(R2, np.float64), extname='VARIANCE')
        f.write(wave.astype(np.float64), extname='WAVE')
        f.write({'STACK_CELL': np.asarray(cells_k).astype('U40')}, extname='CELLS')
    os.replace(out + '.tmp', out)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(5, 1, figsize=(11, 10), sharex=True)
        for k in range(5):
            ax[k].plot(wave[pixok], P[k][pixok], lw=0.5); ax[k].set_ylabel(f'e{k}')
        ax[-1].set_xscale('log'); ax[-1].set_xlabel('rest wavelength [A]'); ax[0].set_title(f'quasar EM-PCA {a.tag}: {len(X)} stacks')
        fig.tight_layout(); fig.savefig(f'{PROD}/qso_empca_{a.tag}.png', dpi=110)
    except Exception as e:
        print('plot skipped:', e)
    last = os.environ.get('EMPCA_SERIES_LAST', '1') == '1'       # after_walk_v2.sh runs 4 bases; 'done' only after the last
    status(3 if last else 2, 3, 'done' if last else 'running', f'{a.tag}: {len(X)} stacks, R2(5)={R2[min(4, len(R2)-1)]:.4f}', t0)
    print(f'wrote {out}, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
