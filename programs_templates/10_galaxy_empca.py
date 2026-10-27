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
Output: <out>/galaxy_empca_<tag>.fits  HDUs EIGVEC [nvec+1, npix], COEFF [nstack, nvec+1], PIXMASK, VARIANCE
        (weighted R^2 vs number of vectors), WAVE, CELLS; and a quick-look PNG.
Usage: python3 10_galaxy_empca.py --stacks ~/data/sdssdesi_passC/stacks --nvec 10 --niter 25 --tag v1
"""
import os, sys, glob, argparse, time
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'programs_stars'))
from step7_empca import empca, explained          # weighted EM-PCA with missing data (stellar thread)


def load_stacks(sdir, min_snr):
    X, W, cells = [], [], []
    for fn in sorted(glob.glob(os.path.join(sdir, 'stacks_z*.fits'))):
        f = fitsio.FITS(fn)
        fl, iv, c, wave = f['FLUX'].read(), f['IVAR'].read(), f['CELLS'].read(), f['WAVE'].read()
        for i in range(len(c)):
            good = np.isfinite(fl[i]) & (iv[i] > 0)
            if good.sum() < 200:
                continue
            snr = np.median(fl[i][good] * np.sqrt(iv[i][good]))
            if not np.isfinite(snr) or snr < min_snr:
                continue
            w = np.where(good, iv[i], 0.0); w = w / np.median(w[good])
            X.append(np.where(good, fl[i], 0.0)); W.append(w)
            cells.append({k: c[k][i] for k in c.dtype.names} | {'STACK_SNR': snr, 'NPIX': int(good.sum())})
    return np.array(X), np.array(W), pd.DataFrame(cells), wave


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', default=os.path.expanduser('~/data/sdssdesi_passC/stacks'))
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC/empca'))
    ap.add_argument('--tag', default='v1'); ap.add_argument('--nvec', type=int, default=10)
    ap.add_argument('--niter', type=int, default=25); ap.add_argument('--min-cover', type=int, default=20)
    ap.add_argument('--min-snr', type=float, default=5.0)
    a = ap.parse_args()
    t0 = time.time()
    X, W, cells, wave = load_stacks(a.stacks, a.min_snr)
    print(f'stacks used: {len(X)} (S/N >= {a.min_snr}), pixels {X.shape[1]}, load {time.time()-t0:.0f}s', flush=True)
    nv = a.nvec + 1                                   # component 0 = mean shape (no centering)
    P, C, pixok = empca(X, W, nv, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    R2 = explained(X, W, P)
    print(f'EM-PCA done {time.time()-t0:.0f}s; pixels in basis {pixok.sum()} ({wave[pixok][0]:.0f}-{wave[pixok][-1]:.0f} A); '
          f'R^2 vs k: ' + ' '.join(f'{k+1}:{r:.4f}' for k, r in enumerate(R2)), flush=True)
    os.makedirs(a.out, exist_ok=True)
    fn = os.path.join(a.out, f'galaxy_empca_{a.tag}.fits')
    with fitsio.FITS(fn + '.tmp', 'rw', clobber=True) as fo:
        fo.write(P.astype(np.float32), extname='EIGVEC'); fo.write(C.astype(np.float32), extname='COEFF')
        fo.write(pixok.astype(np.int16), extname='PIXMASK'); fo.write(R2, extname='VARIANCE')
        fo.write(wave, extname='WAVE'); fo.write(cells.to_records(index=False), extname='CELLS')
    os.replace(fn + '.tmp', fn)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(min(nv, 6), 1, figsize=(10, 1.6 * min(nv, 6)), sharex=True)
        for k in range(min(nv, 6)):
            ax[k].plot(wave, np.where(pixok, P[k], np.nan), lw=0.5); ax[k].set_ylabel(f'e{k}')
        ax[-1].set_xlabel('rest wavelength [A, vacuum]'); ax[-1].set_xscale('log')
        fig.tight_layout(); fig.savefig(fn.replace('.fits', '.png'), dpi=110)
    except Exception as e:
        print('plot skipped:', e)
    print('wrote', fn)


if __name__ == '__main__':
    main()
