#!/usr/bin/env python3
"""Rotate a galaxy basis so that the emission lines sit in a few vectors and the rest are line-free continuum shapes
(Nao 2026-11-03). A change of basis inside the span of the 11 vectors: every stack model is unchanged.
Input: a basis in the reference-window form (EIGVEC e0 = 1 and e1..eN = 0 in REF_LO-REF_HI, COEFF), e.g. galaxy_eigenspectra_v2.fits.
1. Line pixels: the 02_passB_measure.LINES table within +-LINE_KMS.
2. SVD of the line-pixel block of e1..eN, L = E[1:, lines] = U S V^T; the orthogonal rotation U^T sorts the line power
   (S^2) into the leading rotated vectors. n_line = the number of vectors holding LINE_FRAC of the line power (at most MAX_LINE).
3. e0 (the mean galaxy) loses its lines: its line pixels are fitted with the n_line line vectors and the fit is subtracted.
4. Inside the line block and inside the continuum block the vectors are re-ordered by the stack-weighted variance of their
   coefficients (SVD of the coefficients), so e1.. are the line vectors by variance, then the continuum vectors by variance.
Coefficients follow (C' = C T^-1 for E' = T E); the window property (e0 = 1, ek = 0) is kept. Header LINEROT = n_line.
Usage: python3 30_rotate_lines.py IN.fits OUT.fits [--kms 600] [--frac 0.95] [--max-line 4]
"""
import sys, os, argparse, importlib
import numpy as np
import fitsio

CKMS = 299792.458


def line_mask(wave, kms):
    pb = importlib.import_module('02_passB_measure')
    lam = np.array([x[1] for g in pb.LINES.values() for x in g])
    lw = np.log10(wave)[:, None]
    return (np.abs(lw - np.log10(lam)[None, :]) < np.log10(1 + kms / CKMS)).any(1)


def content(E, mask, pixok):
    """rms of each vector over the line pixels / rms over all basis pixels."""
    return np.sqrt((E[:, mask & pixok]**2).mean(1) / np.maximum((E[:, pixok]**2).mean(1), 1e-30))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fin'); ap.add_argument('fout')
    ap.add_argument('--kms', type=float, default=600.0); ap.add_argument('--frac', type=float, default=0.95)
    ap.add_argument('--max-line', type=int, default=4)
    ap.add_argument('--by', default='variance', choices=['variance', 'pixels'],
                    help='variance (default): rank by the line variance of the stacks (coefficients whitened); pixels: SVD of the line pixels')
    a = ap.parse_args()
    F = fitsio.FITS(a.fin)
    E = F['EIGVEC'].read().astype(np.float64); C = F['COEFF'].read().astype(np.float64); hdr = F['EIGVEC'].read_header()
    wave = F['WAVE'].read(); pixok = F['PIXMASK'].read().astype(bool); cells = F['CELLS'].read()
    nv = len(E); mask = line_mask(wave, a.kms) & pixok
    ref = (wave >= hdr['REF_LO']) & (wave <= hdr['REF_HI']) & pixok
    wt = str(hdr.get('STACKWT', 'equal')).strip(); n = cells['N_STACKED'].astype(float); f = n / np.median(n)
    s = np.sqrt(f) if wt == 'sqrtn' else f if wt == 'n' else np.ones(len(n))
    print(f'{a.fin}: {nv} vectors, {mask.sum()} line pixels (+-{a.kms:.0f} km/s), line content before: ' +
          ' '.join(f'{x:.2f}' for x in content(E, mask, pixok)))
    # 2. rotate e1..eN so the line power is concentrated
    if a.by == 'pixels':                              # orthogonal rotation: SVD of the line-pixel block of the vectors
        L = E[1:, mask]; U, S, Vt = np.linalg.svd(L, full_matrices=False); pw = S**2 / (S**2).sum(); T1 = U.T
    else:                                             # by the line VARIANCE of the stacks: whiten the coefficients (Sigma_C = Q Q^T,
        C1 = C[:, 1:]; mu = (s[:, None] * C1).sum(0) / s.sum()      # C'' = C1 Q^-T has unit covariance), then diagonalize
        SigC = ((s[:, None] * (C1 - mu)).T @ (C1 - mu)) / s.sum()   # the line Gram in the whitened coordinates
        lam_c, Vc = np.linalg.eigh(SigC); Q = Vc * np.sqrt(np.maximum(lam_c, 1e-30))[None, :]
        G = E[1:, mask] @ E[1:, mask].T; Gw = Q.T @ G @ Q
        lam, U = np.linalg.eigh(Gw); o = np.argsort(lam)[::-1]; lam, U = lam[o], U[:, o]
        pw = lam / lam.sum(); T1 = U.T @ Q.T          # E' = T1 E1: coefficients C1 T1^-1 have unit covariance, line variance lam_j
    cum = np.cumsum(pw); n_line = min(int(np.searchsorted(cum, a.frac)) + 1, a.max_line)
    print(f'line {"power" if a.by == "pixels" else "variance"} by rotated vector: ' + ' '.join(f'{x:.3f}' for x in pw) +
          f' -> {n_line} line vectors ({cum[n_line-1]:.3f} of the total)')
    T = np.eye(nv); T[1:, 1:] = T1
    # 3. e0 without lines: fit its line pixels with the line vectors
    E1 = T @ E
    b, *_ = np.linalg.lstsq(E1[1:1 + n_line, mask].T, E[0, mask], rcond=None)
    T0 = np.eye(nv); T0[0, 1:1 + n_line] = -b
    T = T0 @ T
    # 4. order each block by the stack-weighted variance of its coefficients
    Ct = C @ np.linalg.inv(T)
    for lo, hi in [(1, 1 + n_line), (1 + n_line, nv)]:
        if hi - lo < 2:
            continue
        if a.by == 'pixels':                          # orthonormal vectors: order by the coefficient variance (SVD of the coefficients)
            B = np.sqrt(s)[:, None] * Ct[:, lo:hi]; _, _, V = np.linalg.svd(B - B.mean(0), full_matrices=False)
            R = np.eye(nv); R[lo:hi, lo:hi] = V        # E'' = R E' (rows lo..hi mixed), C'' = C' V^T
        else:                                         # unit-variance coefficients: PCA of the block's vectors (SVD of E'[lo:hi])
            Eb = (T @ E)[lo:hi][:, pixok]; _, _, V = np.linalg.svd(Eb, full_matrices=False)
            Ub, Sb, _ = np.linalg.svd(Eb, full_matrices=False); R = np.eye(nv); R[lo:hi, lo:hi] = Ub.T
        T = R @ T
    En = T @ E; Cn = C @ np.linalg.inv(T)
    sgn = np.sign((np.sqrt(s)[:, None] * Cn).sum(0)); sgn[sgn == 0] = 1; sgn[0] = 1   # weighted-mean coefficient positive
    En = En * sgn[:, None]; Cn = Cn * sgn[None, :]
    dm = np.abs(Cn @ En - C @ E)[:, pixok].max()
    print(f'line content after: ' + ' '.join(f'{x:.2f}' for x in content(En, mask, pixok)) +
          f'; e0 window mean {En[0, ref].mean():.6f}, max |ek mean| {np.abs(En[1:, ref].mean(1)).max():.1e}; model change {dm:.1e}')
    var = (s[:, None] * (Cn - (s[:, None] * Cn).sum(0) / s.sum())**2).sum(0) / s.sum()
    print('weighted coefficient variance by vector: ' + ' '.join(f'{x:.2e}' for x in var[1:]))
    hdr['LINEROT'] = n_line; hdr['LINE_KMS'] = a.kms
    with fitsio.FITS(a.fout + '.tmp', 'rw', clobber=True) as fo:
        fo.write(En.astype(np.float32), extname='EIGVEC', header=hdr); fo.write(Cn.astype(np.float32), extname='COEFF', header=hdr)
        for ext in [h.get_extname() for h in F][2:]:
            fo.write(F[ext].read(), extname=ext)
    os.replace(a.fout + '.tmp', a.fout); print('wrote', a.fout)


if __name__ == '__main__':
    main()
