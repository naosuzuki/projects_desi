#!/usr/bin/env python3
"""Plot a galaxy EM-PCA basis (10_galaxy_empca.py output), one panel per vector, with the normalization point marked:
a red bar over the reference window (header REF_LO/REF_HI of EIGVEC; e0 = 1 and e1..eN = 0 there), the target level as a
dashed line (1 for e0, 0 for the others) and the window mean as a red dot. Bases without a reference window are drawn as is.
--two-point: the basis is changed within its
span so that both windows are pinned and both are drawn (blue bar and dot: 3000-3100 A) (e0 = 1 / 0, e1 = 0 / 1, e2..eN = 0 / 0 in 5400-5500 / 3000-3100), written to *_2pt.png.
y range -0.1..0.1 for e1..eN (e0, and e1 of the two-point basis: -0.1..1.5).
Usage: python3 24_plot_empca.py BASIS.fits [--nshow 11] [--out PNG]
"""
import argparse
import numpy as np
import fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

YLIM = (-0.1, 0.1)                                         # e1..eN (Nao 2026-11-01)
E0LIM = (-0.1, 1.5)                                        # e0 is 1 in the window, so it keeps its own range


WIN2 = (3000.0, 3100.0)                                    # second normalization window (z >= 0.30 galaxies and stacks)


def two_point(P, w, ok, ref, win2=WIN2):
    """Change of basis within the same span: e0 = 1 / 0, e1 = 0 / 1, e2..eN = 0 / 0 (means in ref / win2).
    P must already be the reference-window basis (e0 = 1, e1..eN = 0 in ref). Returns the new basis."""
    s2 = (w >= win2[0]) & (w <= win2[1]) & ok
    a = P[:, s2].mean(1)
    j = 1 + int(np.argmax(np.abs(a[1:])))                 # vector with the largest 3050 value carries the 3050 point
    e1 = P[j] / a[j]
    Q = P - a[:, None] * e1                                # every vector minus its 3050 value times e1 (ref value unchanged: e1 = 0 there)
    Q[j] = e1
    order = [0, j] + [k for k in range(1, len(P)) if k != j]
    return Q[order]


def plot_basis(fn, out=None, nshow=11, both=False):
    with fitsio.FITS(fn) as f:
        P = f['EIGVEC'].read(); h = f['EIGVEC'].read_header(); w = f['WAVE'].read(); ok = f['PIXMASK'].read().astype(bool)
    ref = (h['REF_LO'], h['REF_HI']) if 'REF_LO' in h else None
    if both and ref:
        P = two_point(P, w, ok, ref)
    n = min(nshow, len(P))
    fig, ax = plt.subplots(n, 1, figsize=(11, 1.5 * n + 0.6), sharex=True)
    for k in range(n):
        ax[k].plot(w, np.where(ok, P[k], np.nan), lw=0.5, color='k'); ax[k].set_ylabel(f'e{k}')
        ax[k].set_ylim(*(E0LIM if k == 0 or (both and k == 1) else YLIM))
        if ref:
            sel = (w >= ref[0]) & (w <= ref[1]) & ok
            ax[k].axvspan(*ref, color='red', alpha=0.35, lw=0)
            ax[k].axhline(1.0 if k == 0 else 0.0, color='red', ls='--', lw=0.7)
            ax[k].plot(np.mean(ref), P[k, sel].mean(), 'o', color='red', ms=4)
            if both:                                       # the second window only where it is pinned (two-point basis)
                s2 = (w >= WIN2[0]) & (w <= WIN2[1]) & ok
                ax[k].axvspan(*WIN2, color='blue', alpha=0.3, lw=0)
                ax[k].plot(np.mean(WIN2), P[k, s2].mean(), 'o', color='blue', ms=4)
                ax[k].axhline(1.0 if k == 1 else 0.0, color='blue', ls=':', lw=0.7)
    ax[-1].set_xscale('log'); ax[-1].set_xlabel('rest wavelength [A, vacuum]')
    if not ref:
        t = '   (no reference window)'
    elif both:
        t = (f'   two-point basis: red {ref[0]:.0f}-{ref[1]:.0f} A (e0 = 1, others 0), blue {WIN2[0]:.0f}-{WIN2[1]:.0f} A (e1 = 1, others 0)')
    else:
        t = (f'   red: normalization window {ref[0]:.0f}-{ref[1]:.0f} A (e0 = 1, e1..eN = 0)')
    fig.suptitle(fn.split('/')[-1] + t, fontsize=10)
    fig.tight_layout(); out = out or fn.replace('.fits', '_2pt.png' if both else '.png'); fig.savefig(out, dpi=100); plt.close(fig)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('basis'); ap.add_argument('--nshow', type=int, default=11); ap.add_argument('--out')
    ap.add_argument('--two-point', action='store_true', help='also pin 3000-3100 A: e0 = 1/0, e1 = 0/1, others 0/0 (5450/3050)')
    a = ap.parse_args(); print(plot_basis(a.basis, a.out, a.nshow, a.two_point))
