#!/usr/bin/env python3
"""Rotate e0 of a galaxy EM-PCA basis to the mean galaxy (Nao 2026-11-03).
Input: galaxy_empca_<tag>.fits from 10_galaxy_empca.py in the reference-window basis (EIGVEC: e0 = 1 and e1..eN = 0 in
REF_LO-REF_HI; COEFF). The weighted mean of the stack models, m = sum_i s_i c_i P / sum_i s_i with the same per-stack weights as
the decomposition (s_i = sqrt(N_i / median N) for STACKWT = sqrtn, N_i / median N for 'n', 1 for 'equal'), lies in the span of the
basis. Since e1..eN have zero window mean, m has window mean a_0 (the weighted mean of c_0), and the change of basis is
    e0' = m / a_0 = e0 + sum_k (a_k / a_0) e_k,   e_k' = e_k (k >= 1),
    c0' = c0,                                       c_k' = c_k - c0 a_k / a_0,
so every model c' P' = c P is unchanged, e0' is the mean galaxy with mean 1 in the window, and the weighted mean of each c_k' is 0.
EIGVEC / COEFF are rewritten, EIGVEC_ORTHO / COEFF_ORTHO (the EM-PCA basis) are kept, header E0ROT = 'mean'.
Usage: python3 29_rotate_e0.py IN.fits OUT.fits
"""
import sys
import numpy as np
import fitsio


def main(fin, fout):
    F = fitsio.FITS(fin)
    P = F['EIGVEC'].read().astype(np.float64); C = F['COEFF'].read().astype(np.float64); hdr = F['EIGVEC'].read_header()
    cells = F['CELLS'].read(); wave = F['WAVE'].read(); pixok = F['PIXMASK'].read().astype(bool)
    wt = str(hdr.get('STACKWT', 'equal')).strip()
    n = cells['N_STACKED'].astype(float); f = n / np.median(n)
    s = np.sqrt(f) if wt == 'sqrtn' else f if wt == 'n' else np.ones(len(n))
    a = (s[:, None] * C).sum(0) / s.sum()
    ref = (wave >= hdr['REF_LO']) & (wave <= hdr['REF_HI']) & pixok
    m = a @ P
    print(f'{fin}: {len(C)} stacks, weights {wt}; mean-galaxy window mean {m[ref].mean():.4f} (a0 {a[0]:.4f}); '
          f'a_k / a_0 for k >= 1: ' + ' '.join(f'{x:+.4f}' for x in a[1:] / a[0]))
    Pn = P.copy(); Pn[0] = m / a[0]
    Cn = C.copy(); Cn[:, 1:] = C[:, 1:] - C[:, [0]] * (a[1:] / a[0])[None, :]
    dm = np.abs(Cn @ Pn - C @ P)[:, pixok].max()
    print(f'e0 window mean {Pn[0, ref].mean():.6f}; max |ek mean| {np.abs(Pn[1:, ref].mean(1)).max():.1e}; model change {dm:.1e}; '
          f'weighted mean of c_k\' (k >= 1) {np.abs((s[:, None] * Cn[:, 1:]).sum(0) / s.sum()).max():.1e}; '
          f'max |e0\' - e0| {np.abs(Pn[0] - P[0])[pixok].max():.4f}')
    hdr['E0ROT'] = 'mean'
    with fitsio.FITS(fout + '.tmp', 'rw', clobber=True) as fo:
        fo.write(Pn.astype(np.float32), extname='EIGVEC', header=hdr); fo.write(Cn.astype(np.float32), extname='COEFF', header=hdr)
        for ext in ['EIGVEC_ORTHO', 'COEFF_ORTHO', 'PIXMASK', 'VARIANCE', 'WAVE', 'CELLS']:
            fo.write(F[ext].read(), extname=ext)
    import os; os.replace(fout + '.tmp', fout); print('wrote', fout)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
