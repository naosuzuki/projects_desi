#!/usr/bin/env python3
"""Galaxy basis with a PASSIVE composite as e0 (Nao 2026-11-03): the reference is the sqrt(N)-weighted mean of the passive stacks
(DB_EFF = 4 and EB <= 2: Dn4000 ~ 2.0, line EW ~ 0), every stack gets one amplitude on it (weighted fit on the line-free pixels),
and e1..eN are the EM-PCA of the residuals. The residual of a star-forming stack is "what star formation adds" (young continuum
plus the lines together), so the vectors come out as nature correlates them instead of a rotation forcing the lines apart.
Same stacks, weights and member cut as 10_galaxy_empca.py; the coverage cut is lower (--min-cover-frac 0.1) to widen the range.
Output in the v2 format: EIGVEC (e0 = passive composite with mean 1 in the reference window, e1..eN with window mean 0),
COEFF (c0 = the passive amplitude), EIGVEC_ORTHO / COEFF_ORTHO (e0 raw + the orthonormal residual vectors), PIXMASK,
VARIANCE (R^2 of the residuals about the passive reference), VARIANCE_TOTAL (R^2 with c0 free), WAVE, CELLS.
Header: STACKWT, PASSIVE = number of passive stacks, E0 = 'passive', REF_LO / REF_HI.
Usage: python3 31_passive_empca.py --stacks DIR DIR --tag v3t_5450_sqrtn_p --out DIR [--min-cover-frac 0.1] [--max-step 0.1 --wave-min 2000]
p2 (Nao 2026-11-04): --max-step drops the stacks with a step on the 3000-3100 A normalization window (faint galaxies normalized per
spectrum in that window at low S/N; the box showed up in e6-e8 of p1), --wave-min 2000 drops the far UV.
"""
import os, sys, time, argparse, importlib
import numpy as np
import fitsio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'programs_stars'))
from step7_empca import empca, explained

G = importlib.import_module('10_galaxy_empca'); RL = importlib.import_module('30_rotate_lines')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', nargs='+', required=True); ap.add_argument('--out', required=True); ap.add_argument('--tag', default='v3t_5450_sqrtn_p')
    ap.add_argument('--nvec', type=int, default=10); ap.add_argument('--niter', type=int, default=25)
    ap.add_argument('--min-cover', type=int, default=150); ap.add_argument('--min-cover-frac', type=float, default=0.1)
    ap.add_argument('--min-member-frac', type=float, default=0.5); ap.add_argument('--min-snr', type=float, default=5.0)
    ap.add_argument('--stack-weight', default='sqrtn', choices=['equal', 'sqrtn', 'n'])
    ap.add_argument('--passive-cover', type=int, default=20, help='e0 is defined where at least this many passive stacks cover the pixel')
    ap.add_argument('--line-kms', type=float, default=600.0); ap.add_argument('--ref-win', type=float, nargs=2, default=[5400.0, 5500.0])
    ap.add_argument('--wave-min', type=float, default=0.0, help='drop pixels bluer than this (Nao 2026-11-04: 2000 A, the far UV is noise)')
    ap.add_argument('--max-step', type=float, default=0.0,
                    help='drop stacks whose weighted mean over the 3000-3100 A normalization window differs from the mean over the 100 A '
                         'on each side by more than this (the window step of a noisy per-spectrum normalization; 0 = keep all)')
    ap.add_argument('--e0', default='', help='npz of 34_restack_select.py: use its E0 (composite of the coefficient-selected passive '
                    'galaxies, on the same rest grid) as the reference instead of the composite of the passive stacks')
    ap.add_argument('--amp', default='fit', choices=['fit', 'unit'],
                    help='fit: one amplitude per stack on the line-free pixels; unit: amplitude 1 for every stack (the stacks are '
                         'normalized to unit window flux and chained), i.e. the residuals are X - e0 as about a mean')
    a = ap.parse_args(); t0 = time.time()
    X, W, cells, wave = G.load_stacks(a.stacks, a.min_snr, a.min_member_frac)
    nst = cells.N_STACKED.to_numpy(float); f = nst / np.median(nst)
    if a.stack_weight == 'sqrtn':
        f = np.sqrt(f)
    elif a.stack_weight == 'equal':
        f = np.ones_like(f)
    W = W * f[:, None]
    if a.max_step > 0:   # window step: stacks of faint galaxies normalized per spectrum in 3000-3100 A show a box on the window
        win = (wave >= 3000) & (wave <= 3100); side = ((wave >= 2900) & (wave < 3000)) | ((wave > 3100) & (wave <= 3200))
        def wm(m):
            return (W[:, m] * X[:, m]).sum(1) / np.maximum(W[:, m].sum(1), 1e-30), (W[:, m] > 0).sum(1)
        mw, nw = wm(win); ms, ns = wm(side); step = np.where((nw > 50) & (ns > 100), mw - ms, 0.0)
        bad = np.abs(step) > a.max_step; cells['WIN_STEP'] = step
        print(f'window step |F(3000-3100) - F(sides)| > {a.max_step}: dropping {bad.sum()} of {len(X)} stacks '
              f'({nst[bad].sum():,.0f} galaxies, {f[bad].sum()/f.sum():.1%} of the weight); median N {np.median(nst[bad]):.0f}, '
              f'z 5/50/95 % {np.percentile(cells.Z_MED[bad], [5, 50, 95]).round(2)}, M_B max {cells.MB[bad].max()}', flush=True)
        X, W, f, nst, cells = X[~bad], W[~bad], f[~bad], nst[~bad], cells[~bad].reset_index(drop=True)
    if a.wave_min > 0:
        W[:, wave < a.wave_min] = 0.0
    ncov = (W > 0).sum(0); a.min_cover = max(a.min_cover, int(np.ceil(a.min_cover_frac * ncov.max())))
    okc = np.where(ncov >= a.min_cover)[0]
    print(f'stacks {len(X)}, min cover {a.min_cover} ({a.min_cover_frac:.0%} of {ncov.max()}): {wave[okc[0]]:.0f}-{wave[okc[-1]]:.0f} A; load {time.time()-t0:.0f}s', flush=True)
    # passive composite
    p = ((cells.DB_EFF == 4) & (cells.EB <= 2)).to_numpy()
    num = (W[p] * X[p]).sum(0); den = W[p].sum(0); pcov = (W[p] > 0).sum(0)
    e0 = np.where(pcov >= a.passive_cover, num / np.maximum(den, 1e-30), 0.0); e0def = pcov >= a.passive_cover
    if a.e0:
        zz = np.load(a.e0); assert len(zz['WAVE']) == len(wave) and np.allclose(zz['WAVE'], wave, rtol=1e-6)
        e0 = zz['E0'].astype(np.float64); e0def = (zz['E0_COV'] >= 5) & (e0 != 0); e0[~e0def] = 0.0
        print(f'external e0 from {a.e0}: {len(zz["N"])} cells, {zz["N"].sum():,} galaxies, selection {zz["SEL"]}', flush=True)
    ref = (wave >= a.ref_win[0]) & (wave <= a.ref_win[1]) & e0def
    e0 = e0 / e0[ref].mean()
    jd = np.where(e0def)[0]
    print(f'passive composite: {p.sum()} stacks, {nst[p].sum():,.0f} galaxies, z 5/50/95 % {np.percentile(cells.Z_MED[p], [5, 50, 95]).round(2)}, '
          f'defined {wave[jd[0]]:.0f}-{wave[jd[-1]]:.0f} A ({e0def.sum()} pixels)', flush=True)
    # one amplitude per stack on the line-free pixels
    lines = RL.line_mask(wave, a.line_kms)
    Wa = np.where((~lines & e0def)[None, :], W, 0.0)
    amp = (Wa * X * e0).sum(1) / np.maximum((Wa * e0**2).sum(1), 1e-30)
    if a.amp == 'unit':
        print(f'fitted passive amplitude (not used) 5/50/95 % {np.percentile(amp, [5, 50, 95]).round(3)}', flush=True); amp = np.ones(len(X))
    print(f'passive amplitude 5/50/95 % {np.percentile(amp, [5, 50, 95]).round(3)}; by DB_EFF: ' +
          ' '.join(f'{d}:{np.median(amp[cells.DB_EFF == d]):.3f}' for d in sorted(cells.DB_EFF.unique())), flush=True)
    R = np.where(W > 0, X - amp[:, None] * e0[None, :], 0.0)
    Pr, Cr, pixok = empca(R, W, a.nvec, a.niter, min_cover=a.min_cover, min_frac=0.0, min_wfrac=0.0)
    R2p = explained(R[:, pixok], W[:, pixok], Pr[:, pixok])
    P = np.vstack([e0[None, :], Pr]); C = np.hstack([amp[:, None], Cr])
    R2 = explained(X[:, pixok], W[:, pixok], P[:, pixok])
    print(f'EM-PCA done {time.time()-t0:.0f}s; pixels in basis {pixok.sum()} ({wave[pixok][0]:.0f}-{wave[pixok][-1]:.0f} A)', flush=True)
    print('R^2 of the residuals about the passive reference vs k: ' + ' '.join(f'{k+1}:{r:.4f}' for k, r in enumerate(R2p)), flush=True)
    print('R^2 with c0 free vs k: ' + ' '.join(f'{k+1}:{r:.4f}' for k, r in enumerate(R2)), flush=True)
    # reference-window basis (as 10_galaxy_empca): e0 = 1, ek = 0 in the window
    refw = (wave >= a.ref_win[0]) & (wave <= a.ref_win[1]) & pixok
    v = P[:, refw].mean(1); Pn = P.copy(); Pn[0] = P[0] / v[0]; Pn[1:] = P[1:] - v[1:, None] * Pn[0]; Cn = C.copy(); Cn[:, 0] = C @ v
    print(f'reference window: e0 mean {Pn[0, refw].mean():.4f}, max |ek mean| {np.abs(Pn[1:, refw].mean(1)).max():.1e}; model change {np.abs(Cn @ Pn - C @ P).max():.1e}', flush=True)
    cells['EMPCA_SCALE'] = 1.0
    cells['TO_STACK'] = 1.0 / (cells.RENORM * cells.CHAIN) if {'RENORM', 'CHAIN'} <= set(cells) else 1.0
    os.makedirs(a.out, exist_ok=True); fn = os.path.join(a.out, f'galaxy_empca_{a.tag}.fits')
    hdr = {'REF_LO': a.ref_win[0], 'REF_HI': a.ref_win[1], 'STACKWT': a.stack_weight, 'CENTER': 0, 'PASSIVE': int(p.sum()), 'E0': 'select' if a.e0 else 'passive',
           'LINE_KMS': a.line_kms, 'MINCOVF': a.min_cover_frac, 'AMP': a.amp, 'WAVEMIN': a.wave_min, 'MAXSTEP': a.max_step, 'E0FILE': os.path.basename(a.e0) if a.e0 else ''}
    with fitsio.FITS(fn + '.tmp', 'rw', clobber=True) as fo:
        fo.write(Pn.astype(np.float32), extname='EIGVEC', header=hdr); fo.write(Cn.astype(np.float32), extname='COEFF', header=hdr)
        fo.write(P.astype(np.float32), extname='EIGVEC_ORTHO'); fo.write(C.astype(np.float32), extname='COEFF_ORTHO')
        fo.write(pixok.astype(np.int16), extname='PIXMASK'); fo.write(R2p, extname='VARIANCE'); fo.write(R2, extname='VARIANCE_TOTAL')
        fo.write(wave, extname='WAVE'); fo.write(cells.to_records(index=False), extname='CELLS')
    os.replace(fn + '.tmp', fn)
    try:
        importlib.import_module('24_plot_empca').plot_basis(fn)
    except Exception as e:
        print('plot skipped:', e)
    print('wrote', fn, flush=True)


if __name__ == '__main__':
    main()
