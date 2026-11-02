#!/usr/bin/env python3
"""Quasar eigenspectra: weighted EM-PCA (Bailey 2012) of the quasar stacks (05_stack.py --merge -> qso_stacks*.npz).
Same conventions as the galaxy thread (programs_templates/10_galaxy_empca.py), routine programs_stars/step7_empca.empca:
- inputs: the ivar-weighted mean stacks WMEAN (unit mean in the normalization window of their z slice), cells with
  at least --nmin quasars at their best-covered pixel;
- no mean subtraction (coverage depends on z): the mean shape is the first of nvec+1 vectors;
- each stack counts once: its weights 1/WERR^2 are rescaled to a median of 1 over its good pixels;
- a pixel enters the basis where >= --min-cover stacks cover it and, over the kept stacks, >= --min-total member spectra
  (default 10 000; Nao 2026-11-01: the sparse blue and red edges are dropped);
- v3 (no z axis, members of all redshifts): a stack pixel is dropped where fewer than max(--min-member,
  --min-member-frac x the stack's best-covered count) members cover it (sparse edges, which also wreck the median weight).
Pivot (Nao 2026-10-31): every stack is scaled to unit mean in the pivot window (default 3000-3050 A, the v3 chaining
window; a stack that does not cover it keeps its chained normalization), e0 is the weighted mean spectrum scaled to unity
there, the EM-PCA runs on the residuals X - e0, and each eigenvector is then shifted, e_k -> e_k - <e_k>_pivot e0, so that
all e_k (k >= 1) are exactly zero at the pivot (the mean over the window) while e0 is exactly one. This is a change of basis
inside the same span (c0 absorbs the shift), so the fits are unchanged; the e_k are no longer orthonormal. --no-pivot gives
the old basis (first vector of the un-centred data). PIVOT_LO/PIVOT_HI in the EIGVEC header.
Output: ~/data/sdssdesi_products/quasar/qso_empca_<tag>.fits: EIGVEC [nvec+1, npix], COEFF [nstack, nvec+1], PIXMASK,
VARIANCE (weighted R^2 vs number of vectors), WAVE, CELLS (+ qso_empca_<tag>.png quick look).
Usage: python3 06_empca.py --stacks ~/data/sdssdesi_products/quasar/qso_stacks.npz --nvec 10 --niter 25 --tag v1
       python3 06_empca.py ... --zrange 0.8 2.0 --tag v1_z0820     (one basis per normalization window)
"""
import os, sys, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qsolib as Q
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'programs_stars'))
from step7_empca import empca, explained          # weighted EM-PCA with missing data (stellar thread, stable interface)

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
STATUS_V1 = os.path.expanduser('~/data/sdssdesi_status/step_qso_empca.json')
STATUS_V2 = os.path.expanduser('~/data/sdssdesi_status/step_qso_empca_v2.json')   # own gauge row, so v2 never overwrites v1
STATUS_V3 = os.path.expanduser('~/data/sdssdesi_status/step_qso_empca_v3.json')
TAG = ''
NO_STATUS = False


def status(done, total, st, note, t0):
    if NO_STATUS: return
    STATUS = STATUS_V2 if TAG.startswith('v2') else (STATUS_V3 if TAG.startswith('v3') else STATUS_V1)
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
    ap.add_argument('--min-member', type=int, default=0, help='drop stack pixels covered by fewer members (v3: 30)')
    ap.add_argument('--min-member-frac', type=float, default=0.0, help="... or by fewer than this fraction of the stack's best-covered count (v3: 0.1)")
    ap.add_argument('--zrange', type=float, nargs=2, default=None, help='only stacks whose z slice centre is in [lo, hi): one basis per normalization window')
    ap.add_argument('--pivot', type=float, nargs=2, default=[3000.0, 3050.0], help='pivot window [A]: e0 = 1 and e_k = 0 there (window mean)')
    ap.add_argument('--no-pivot', action='store_true', help='old basis: nvec+1 vectors of the un-centred stacks, no pivot')
    ap.add_argument('--min-total', type=int, default=10000, help='drop the blue and red edges: pixels covered by fewer member spectra in all (sum over the kept stacks of their good members); 0 = off')
    ap.add_argument('--span', type=float, nargs=2, default=None, help='test: only stacks that cover [lo, hi] A (both sides of the pivot)')
    ap.add_argument('--no-status', action='store_true', help='test runs: do not touch the monitor gauge')
    ap.add_argument('--snmin', type=float, default=5.0, help='per stack, drop the blue end until the running median S/N per pixel (100 px) reaches this; 0 = off (Nao 2026-11-01)')
    ap.add_argument('--drop', type=float, default=0.9, help='... and until the running mean flux (100 px) is >= this x the 200 px redward; 0 = off')
    ap.add_argument('--max-mag', type=float, default=-19.0, help='drop stacks whose luminosity bin lies entirely fainter than this M_AB(3025) (label L; v4, Nao 2026-11-01: the faintest bins are host galaxies, not quasars; -19 since v5, Nao 2026-11-02)')
    a = ap.parse_args()
    global TAG, NO_STATUS; TAG = a.tag; NO_STATUS = a.no_status
    t0 = time.time()
    status(0, 3, 'running', f'{a.tag}: loading stacks', t0)
    z = np.load(a.stacks, allow_pickle=True)
    wave, wm, we, n, cells = z['WAVE'], z['WMEAN'], z['WERR'], z['N'], z['STACK_CELL']
    keep = n.max(axis=1) >= a.nmin
    if a.max_mag is not None and str(cells[0]).startswith('L'):      # v3 labels L_E_S_O_B, L = floor(log F_abs(3025) / 0.4)
        lmin = int(np.floor(-(a.max_mag + 19.81) / 2.5 / 0.4)) + 1       # M_AB(3025) = -2.5 log F - 19.81
        lb = np.array([int(str(c).split('_')[0][1:]) for c in cells]); keep &= lb >= lmin
        print(f'stacks fainter than M_AB(3025) = {a.max_mag:g} (L < {lmin}) dropped: {(lb < lmin).sum()}', flush=True)
    if a.zrange:
        import importlib; ZE = np.asarray(importlib.import_module('04_cells').ZEDGES, float)
        zb = np.array([int(str(c).split('_')[0][1:]) for c in cells]); zc = 0.5 * (ZE[zb] + ZE[zb + 1])
        keep &= (zc >= a.zrange[0]) & (zc < a.zrange[1])
    cov = np.isfinite(wm) & np.isfinite(we) & (we > 0) & (n >= np.maximum(a.min_member, a.min_member_frac * n.max(axis=1))[:, None])
    if a.snmin > 0 or a.drop > 0:                          # Nao 2026-11-01: per-stack blue ends dropped where S/N < snmin or the flux falls toward the blue
        before = cov.sum(axis=1); cov = Q.blue_trim(wave, wm, we, cov, a.snmin, a.drop)
        print(f'blue ends trimmed in {((cov.sum(axis=1) < before) & keep).sum()} of {keep.sum()} stacks', flush=True)
    X = np.where(cov, wm, 0.0)[keep].astype(np.float64)
    W = np.where(cov, 1.0 / np.where(we > 0, we, 1.0)**2, 0.0)[keep]
    good = (W > 0).sum(axis=1) >= a.min_pix
    if a.span:
        sp = (wave >= a.span[0]) & (wave <= a.span[1]); good &= (W[:, sp] > 0).mean(axis=1) >= 0.95
    X, W, cells_k = X[good], W[good], cells[keep][good]
    if a.min_total > 0:                                   # Nao 2026-11-01: drop the sparse blue and red edges
        tot = (n[keep][good] * (W > 0)).sum(axis=0); edge = tot < a.min_total
        X[:, edge] = 0.0; W[:, edge] = 0.0
        print(f'edges with < {a.min_total:,} spectra per pixel dropped: basis spans {wave[~edge].min():.0f}-{wave[~edge].max():.0f} A', flush=True)
    for i in range(len(W)):
        W[i] /= np.median(W[i][W[i] > 0])
    print(f'{len(X)} stacks (n >= {a.nmin}, >= {a.min_pix} good pixels) of {len(cells)}; {X.shape[1]} pixels', flush=True)
    status(1, 3, 'running', f'{a.tag}: EM-PCA on {len(X)} stacks', t0)
    nv = a.nvec + 1; hdr = {}
    if a.no_pivot:
        P, C, pixok = empca(X, W, nv, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    else:
        piv = (wave >= a.pivot[0]) & (wave <= a.pivot[1])
        for i in range(len(X)):                                          # exact unit mean in the pivot window where covered
            ok = piv & (W[i] > 0)
            if ok.sum() >= 0.5 * piv.sum():
                X[i] /= X[i][ok].mean()
        with np.errstate(invalid='ignore'):
            e0 = np.where(W.sum(axis=0) > 0, (W * X).sum(axis=0) / W.sum(axis=0), 0.0)
        e0 /= e0[piv & (W.sum(axis=0) > 0)].mean()                          # e0 = 1 at the pivot
        Pr, Cr, pixok = empca(X - e0, W, a.nvec, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
        e0 = np.where(pixok, e0, 0.0)
        sh = np.array([Pr[k][piv & pixok].mean() for k in range(a.nvec)])  # <e_k>_pivot
        P = np.vstack([e0, Pr - sh[:, None] * e0[None, :]])                # e_k = 0 at the pivot, same span
        C = np.hstack([(1.0 + Cr @ sh)[:, None], Cr])                       # X ~ e0 + sum c_k e_k^old = (1 + c.sh) e0 + sum c_k e_k
        hdr = dict(PIVOT_LO=a.pivot[0], PIVOT_HI=a.pivot[1], MINTOTAL=a.min_total, SNMIN=a.snmin, DROP=a.drop, MAXMAG=a.max_mag)
        print(f'pivot {a.pivot[0]:.0f}-{a.pivot[1]:.0f} A: e0 mean {P[0][piv].mean():.6f}, |e_k| there <= {np.abs(P[1:][:, piv].mean(axis=1)).max():.1e}', flush=True)
    R2 = explained(X, W, P)
    print('weighted R^2 vs number of vectors:', ' '.join(f'{r:.4f}' for r in R2), flush=True)
    out = f'{PROD}/qso_empca_{a.tag}.fits'
    with fitsio.FITS(out + '.tmp', 'rw', clobber=True) as f:
        f.write(P.astype(np.float32), extname='EIGVEC', header=hdr)
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
            if not a.no_pivot:
                ax[k].axvspan(*a.pivot, color='C0', alpha=0.25); ax[k].axhline(1.0 if k == 0 else 0.0, color='k', lw=0.3)
        ax[-1].set_xscale('log'); ax[-1].set_xlabel('rest wavelength [A]'); ax[0].set_title(f'quasar EM-PCA {a.tag}: {len(X)} stacks')
        fig.tight_layout(); fig.savefig(f'{PROD}/qso_empca_{a.tag}.png', dpi=110)
    except Exception as e:
        print('plot skipped:', e)
    last = os.environ.get('EMPCA_SERIES_LAST', '1') == '1'       # after_walk_v2.sh runs 4 bases; 'done' only after the last
    status(3 if last else 2, 3, 'done' if last else 'running', f'{a.tag}: {len(X)} stacks, R2(5)={R2[min(4, len(R2)-1)]:.4f}', t0)
    print(f'wrote {out}, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
