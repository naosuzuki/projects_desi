#!/usr/bin/env python3
"""Eigenspectra figure of a quasar basis, e0 to e<nvec-1> (Nao 2026-11-02: 'eigenspectra plots from e0 to e9').
Usage: python3 15_plot_eigvec.py [--basis ../data_eigenspectra/quasar_eigenspectra_v2.fits] [--nvec 10] [--out PNG]
One panel per vector on a log wavelength axis, the pivot window (header PIVOT_LO/HI) shaded, the variance share of each
vector in its label (VARIANCE = cumulative weighted R^2 with 1..k vectors; the label shows the gain of e_k); e0 is the mean shape (unity in the pivot window). Default output: next to the basis, <basis>_e0-e<n-1>.png."""
import os, sys, argparse, numpy as np, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LINES = {'Lya': 1215.67, 'CIV': 1549.06, 'CIII]': 1908.73, 'MgII': 2798.75, '[OII]': 3728.48, 'Hb': 4862.68, '[OIII]': 5008.24, 'Ha': 6564.61}

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--basis', default=os.path.join(HERE, '..', 'data_eigenspectra', 'quasar_eigenspectra_v3.fits'))
    ap.add_argument('--nvec', type=int, default=10); ap.add_argument('--out', default=None); ap.add_argument('--title', default=None)
    a = ap.parse_args()
    E = fitsio.read(a.basis, ext='EIGVEC'); h = fitsio.read_header(a.basis, ext='EIGVEC'); wave = fitsio.read(a.basis, ext='WAVE')
    ok = fitsio.read(a.basis, ext='PIXMASK') > 0; var = fitsio.read(a.basis, ext='VARIANCE'); nst = len(fitsio.read(a.basis, ext='COEFF'))
    n = min(a.nvec, len(E)); piv = (h.get('PIVOT_LO'), h.get('PIVOT_HI')); gain = np.diff(np.r_[0.0, var])      # VARIANCE is the cumulative R^2 with 1..k+1 vectors
    fig, ax = plt.subplots(n, 1, figsize=(11, 1.9 * n), sharex=True)
    for k in range(n):
        ax[k].plot(wave[ok], E[k][ok], lw=0.5, color='C0')
        ax[k].axhline(1.0 if k == 0 else 0.0, color='k', lw=0.3)
        if piv[0] is not None: ax[k].axvspan(piv[0], piv[1], color='C0', alpha=0.25)
        for name, l in LINES.items():
            if wave[ok].min() < l < wave[ok].max():
                ax[k].axvline(l, color='0.7', lw=0.3, zorder=0)
                if k == 0: ax[k].text(l, 0.97, name, transform=ax[k].get_xaxis_transform(), ha='center', va='top', fontsize=7, color='0.4')
        ax[k].set_ylabel(f'e{k}\nR$^2$ {100 * var[k]:.1f} %' + (f'\n+{100 * gain[k]:.2f}' if k > 0 else ''), fontsize=8)
        ax[k].margins(x=0.01)
    ax[-1].set_xscale('log'); ax[-1].set_xlabel('rest wavelength [A]')
    ticks = [800, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000]; ticks = [t for t in ticks if wave[ok].min() <= t <= wave[ok].max()]
    ax[-1].set_xticks(ticks); ax[-1].set_xticklabels([str(t) for t in ticks]); ax[-1].minorticks_off()
    ax[0].set_title(a.title or f'{os.path.basename(a.basis)}: {nst} stacks, {wave[ok].min():.0f}-{wave[ok].max():.0f} A; e0 = mean shape, R$^2$ = variance explained by e0..e$_k$')
    out = a.out or os.path.splitext(a.basis)[0] + f'_e0-e{n - 1}.png'
    fig.tight_layout(); fig.savefig(out, dpi=100); print(out)
