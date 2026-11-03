#!/usr/bin/env python3
"""Galaxy eigenspectra: weighted EM-PCA (Bailey 2012) of the pass-C cell stacks (ms.tex Sec. 3.6 items 3-6).

Uses the EM-PCA routine of the stellar thread (programs_stars/step7_empca.empca) with galaxy settings:
- the stacks cover different rest-frame windows (1400-9900 A over 0<z<1.6), so a pixel enters the basis when at
  least --min-cover stacks cover it (absolute count), not a fraction of all stacks;
- the data are NOT mean-subtracted: with redshift-dependent coverage a per-pixel mean is biased by which cells
  cover which pixel, so the mean shape is fitted jointly as the first component (nvec+1 vectors), as in Redrock;
- inputs are the normalized stacks (unit mean F_lambda in the line-free window of their z slice, normwin.py);
  each stack counts once: its ivar is rescaled to a median of 1 over its good pixels (as for the stars), so the
  richest low-z cells do not dominate the basis.
- reference window (--ref-win, default 5400-5500 A): the basis is changed within its span so that e0 has mean 1 and
  e1..eN mean 0 there; COEFF[:, 0] is then each stack's model flux in the window (its normalization).
Output: <out>/galaxy_empca_<tag>.fits  HDUs EIGVEC [nvec+1, npix] (reference-window basis), COEFF [nstack, nvec+1],
        EIGVEC_ORTHO / COEFF_ORTHO (the orthonormal EM-PCA basis), PIXMASK, VARIANCE (weighted R^2 vs number of
        vectors), WAVE, CELLS (EMPCA_SCALE, TO_STACK); and a quick-look PNG.
Usage: python3 10_galaxy_empca.py --stacks ~/data/sdssdesi_passC/stacks --nvec 10 --niter 25 --tag v1
"""
import os, sys, glob, argparse, time
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'programs_stars'))
from step7_empca import empca, explained          # weighted EM-PCA with missing data (stellar thread)


# telluric bands (observed A): the rest pixels a stack's z slice maps onto them are dropped from that stack (--mask-telluric).
# Narrow sky lines (OI 5577, 6300, Na D) are not masked: within a dz = 0.05 slice each galaxy puts them on a different rest pixel,
# so a stack carries ~2 % of them at most (v3 check: no observed-frame excess, 25_obs_residuals.py).
TELLURIC = ((6860.0, 6960.0), (7590.0, 7700.0), (9300.0, 9650.0))
DZ = 0.05


def load_stacks(sdirs, min_snr, min_member_frac=0.0, mask_telluric=False):
    X, W, cells = [], [], []
    for sdir in sdirs:
      src = 'sdss' if 'passG_sdss' in sdir else 'desi'
      for fn in sorted(glob.glob(os.path.join(sdir, 'stacks_z*.fits'))):
          f = fitsio.FITS(fn)
          fl, iv, c, wave = f['FLUX'].read(), f['IVAR'].read(), f['CELLS'].read(), f['WAVE'].read()
          ns = f['NSPEC'].read() if min_member_frac > 0 else None
          for i in range(len(c)):
              good = np.isfinite(fl[i]) & (iv[i] > 0)
              if ns is not None:                      # drop the stack's edge pixels covered by few of its members
                  good &= ns[i] >= min_member_frac * max(int(c['N_STACKED'][i]), 1)
              if mask_telluric:                       # rest range of each band over the stack's z slice
                  zlo, zhi = int(c['ZB'][i]) * DZ, (int(c['ZB'][i]) + 1) * DZ
                  for lo, hi in TELLURIC:
                      good &= ~((wave >= lo / (1 + zhi)) & (wave <= hi / (1 + zlo)))
              if good.sum() < 200:
                  continue
              snr = np.median(fl[i][good] * np.sqrt(iv[i][good]))
              if not np.isfinite(snr) or snr < min_snr:
                  continue
              w = np.where(good, iv[i], 0.0); w = w / np.median(w[good])
              X.append(np.where(good, fl[i], 0.0)); W.append(w)
              cells.append({k: c[k][i] for k in c.dtype.names} | {'STACK_SNR': snr, 'NPIX': int(good.sum()), 'SOURCE': src})
    return np.array(X), np.array(W), pd.DataFrame(cells), wave


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', nargs='+', default=[os.path.expanduser('~/data/sdssdesi_passC/stacks')],
                    help='one or more stack directories (09_merge_stacks.py output); each stack enters once')
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC/empca'))
    ap.add_argument('--tag', default='v1'); ap.add_argument('--nvec', type=int, default=10)
    ap.add_argument('--niter', type=int, default=25); ap.add_argument('--min-cover', type=int, default=20)
    ap.add_argument('--min-snr', type=float, default=5.0)
    ap.add_argument('--mask-telluric', action='store_true', help='drop each stack\'s rest pixels on the A, B and 9300-9650 A water bands')
    ap.add_argument('--min-cover-frac', type=float, default=0.0, help='also require this fraction of the best-covered pixel\'s stack count (drops thin edges)')
    ap.add_argument('--min-member-frac', type=float, default=0.0, help='per stack, keep pixels with NSPEC >= this x N_STACKED')
    ap.add_argument('--renorm', type=int, default=1, help='passes that rescale each stack to the mean shape e0 (0 = off)')
    ap.add_argument('--ref-win', type=float, nargs=2, default=[5400.0, 5500.0],
                    help='reference window: e0 has mean 1 there and e1..eN mean 0 (same span, so the same fits)')
    ap.add_argument('--stack-weight', default='equal', choices=['equal', 'sqrtn', 'n'],
                    help='relative weight of the stacks: equal (each stack at median weight 1; v1-v3t), x sqrt(N_STACKED / median) '
                         '(the Quasar thread\'s recipe, Nao 2026-11-03), or x N_STACKED / median')
    ap.add_argument('--center', action='store_true',
                    help='classic PCA (Nao 2026-11-03): e0 = the weighted mean of the stacks (coefficient 1), e1..eN = EM-PCA of the '
                         'residuals about it; VARIANCE = R^2 of the variance about the mean, VARIANCE_TOTAL = R^2 with c0 free')
    a = ap.parse_args()
    t0 = time.time()
    X, W, cells, wave = load_stacks(a.stacks, a.min_snr, a.min_member_frac, a.mask_telluric)
    if a.stack_weight != 'equal':                      # Nao 2026-11-03: as in programs_quasars/06_empca.py --stack-weight sqrtn
        nst = cells.N_STACKED.to_numpy(float); f = nst / np.median(nst)
        if a.stack_weight == 'sqrtn':
            f = np.sqrt(f)
        W = W * f[:, None]
        print(f'stack weight {a.stack_weight}: factor 5/50/95 % {np.percentile(f, [5, 50, 95]).round(3)} (N_STACKED median {np.median(nst):.0f})', flush=True)
    ncov = (W > 0).sum(0)
    if a.min_cover_frac > 0:
        a.min_cover = max(a.min_cover, int(np.ceil(a.min_cover_frac * ncov.max())))
        okc = np.where(ncov >= a.min_cover)[0]
        print(f'min cover {a.min_cover} stacks ({a.min_cover_frac:.0%} of {ncov.max()}): {wave[okc[0]]:.0f}-{wave[okc[-1]]:.0f} A', flush=True)
    print(f'stacks used: {len(X)} (S/N >= {a.min_snr}), pixels {X.shape[1]}, load {time.time()-t0:.0f}s', flush=True)
    nv = a.nvec + 1                                   # component 0 = mean shape (no centering)
    P, C, pixok = empca(X, W, nv, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    # common amplitude: the z-slice windows (normwin.py) give each slice its own flux scale, which leaks into the
    # coefficients (c0 jumps x3 at z=0.60 in desi_v1). Rescale each stack so that its weighted projection on the mean
    # shape e0 is 1 over the pixels it covers, then redo the decomposition.
    scale = np.ones(len(X))
    for r in range(a.renorm):
        e0 = P[0]; s = (W * X * e0).sum(1) / np.maximum((W * e0**2).sum(1), 1e-30)
        s = s * np.sign(np.median(s)); s = np.where(np.isfinite(s) & (s > 0), s, np.nan)
        keep = np.isfinite(s); X, W, cells, s = X[keep], W[keep], cells[keep].reset_index(drop=True), s[keep]
        X = X / s[:, None]; scale = scale[keep] * s          # W unchanged: each stack keeps median weight 1 (counts once)
        print(f'renorm pass {r+1}: stacks {len(X)}, scale 5/50/95% {np.percentile(s, [5, 50, 95]).round(3)}', flush=True)
        P, C, pixok = empca(X, W, nv, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    # every factor is kept: EMPCA_SCALE here; RENORM and CHAIN from the chaining step (19/23_chain) are left untouched.
    # model of the stack as stored by 09_merge_stacks = (C @ P) * TO_STACK, TO_STACK = EMPCA_SCALE / (RENORM * CHAIN)
    cells['EMPCA_SCALE'] = scale
    cells['TO_STACK'] = scale / (cells.RENORM * cells.CHAIN) if {'RENORM', 'CHAIN'} <= set(cells) else scale
    R2c = None
    if a.center:                                      # Nao 2026-11-03: classic PCA about the weighted mean of the (rescaled) stacks
        m = (W * X).sum(0) / np.maximum(W.sum(0), 1e-30); m[~pixok] = 0.0
        Xc = np.where(W > 0, X - m[None, :], 0.0)
        Pr, Cr, pixok_r = empca(Xc, W, a.nvec, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
        pixok &= pixok_r
        R2c = explained(Xc[:, pixok], W[:, pixok], Pr[:, pixok])
        P = np.vstack([m[None, :], Pr]); C = np.hstack([np.ones((len(X), 1)), Cr])
        print(f'centered: R^2 of the variance about the mean vs k: ' + ' '.join(f'{k+1}:{r:.4f}' for k, r in enumerate(R2c)), flush=True)
    R2 = explained(X[:, pixok], W[:, pixok], P[:, pixok])        # over the basis pixels (dropped edges are not 'unexplained')
    print(f'EM-PCA done {time.time()-t0:.0f}s; pixels in basis {pixok.sum()} ({wave[pixok][0]:.0f}-{wave[pixok][-1]:.0f} A); '
          f'R^2 vs k: ' + ' '.join(f'{k+1}:{r:.4f}' for k, r in enumerate(R2)), flush=True)
    # reference point (Nao 2026-10-31): change basis inside the same span so that e0 = 1 and e1..eN = 0 (window means) in
    # --ref-win. Then COEFF[:, 0] is the model flux in the window (the normalization) and e1..eN change the shape only.
    # P_k = e'_k + v_k e'_0 with v = window means of P, so C'_0 = C @ v and C'_k = C_k (k >= 1). The new vectors are not
    # orthogonal; EIGVEC_ORTHO / COEFF_ORTHO keep the EM-PCA ones.
    ref = (wave >= a.ref_win[0]) & (wave <= a.ref_win[1]) & pixok
    if ref.sum() < 5:
        raise SystemExit(f'reference window {a.ref_win} has {ref.sum()} basis pixels')
    v = P[:, ref].mean(1)
    Pn = P.copy(); Pn[0] = P[0] / v[0]; Pn[1:] = P[1:] - v[1:, None] * Pn[0]
    Cn = C.copy(); Cn[:, 0] = C @ v
    print(f'reference window {a.ref_win[0]:.0f}-{a.ref_win[1]:.0f} A ({ref.sum()} pixels): e0 mean {Pn[0, ref].mean():.4f}, '
          f'max |ek mean| {np.abs(Pn[1:, ref].mean(1)).max():.1e}; model change {np.abs(Cn @ Pn - C @ P).max():.1e}', flush=True)
    os.makedirs(a.out, exist_ok=True)
    fn = os.path.join(a.out, f'galaxy_empca_{a.tag}.fits')
    hdr = {'REF_LO': a.ref_win[0], 'REF_HI': a.ref_win[1], 'STACKWT': a.stack_weight, 'CENTER': int(a.center)}
    with fitsio.FITS(fn + '.tmp', 'rw', clobber=True) as fo:
        fo.write(Pn.astype(np.float32), extname='EIGVEC', header=hdr); fo.write(Cn.astype(np.float32), extname='COEFF', header=hdr)
        fo.write(P.astype(np.float32), extname='EIGVEC_ORTHO'); fo.write(C.astype(np.float32), extname='COEFF_ORTHO')
        fo.write(pixok.astype(np.int16), extname='PIXMASK')
        fo.write(R2c if R2c is not None else R2, extname='VARIANCE'); fo.write(R2, extname='VARIANCE_TOTAL')
        fo.write(wave, extname='WAVE'); fo.write(cells.to_records(index=False), extname='CELLS')
    os.replace(fn + '.tmp', fn)
    try:                                              # quick look with the normalization window marked (24_plot_empca.py)
        import importlib; importlib.import_module('24_plot_empca').plot_basis(fn)
    except Exception as e:
        print('plot skipped:', e)
    print('wrote', fn)


if __name__ == '__main__':
    main()
