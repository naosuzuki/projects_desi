#!/usr/bin/env python3
"""Figures of the released stellar bases, telluric-corrected with sqrt(N) stack weights (Nao 2026-11-03), sized for ms.tex.
Luminosity bins: AB magnitude per class, DESI top, SDSS below (plot_lum.binned of the corrected stacks). Eigenspectra: e0 (mean, G flux = 1) to e9, one column per class, telluric bands shaded; the label is the gain in R^2 of that vector (EM-PCA does not sort by it).
Output: plots_luminosity/plots_star_lum_tell.png, plots_star_empca_v4_<survey>.png; the latter copied to data_eigenspectra/star_eigenspectra_v4_<survey>_e0-e9.png by hand.
"""
import os, sys
import numpy as np, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path[:0] = [HERE, os.path.join(ROOT, 'programs_telluric')]
import plot_lum, tellapply
OUT = os.path.join(ROOT, 'plots_luminosity'); REL = os.path.join(ROOT, 'data_eigenspectra'); WAVE = 10**(3.5502 + 1e-4 * np.arange(4447)); CLASSES = {'desi': ['WD', 'HOT', 'AF', 'GK', 'M'], 'sdss': ['WD', 'AF', 'GK', 'M']}
BLUE, GRAY, BAND, INK, MUTED = '#2a78d6', '#a8a79f', '#ecebe6', '#0b0b0b', '#52514e'


def style(q, size=6):
    q.tick_params(labelsize=size, colors=MUTED, length=2, pad=1.5, width=0.5); [q.spines[s].set_visible(False) for s in ('top', 'right')]; [q.spines[s].set_color(GRAY) for s in ('left', 'bottom')]; [q.spines[s].set_linewidth(0.5) for s in ('left', 'bottom')]


def lum():
    all_ = ['WD', 'HOT', 'AF', 'GK', 'M']; cm = plt.get_cmap('rainbow_r'); norm = plt.Normalize(*plot_lum.MRANGE)
    fig, ax = plt.subplots(2, 5, figsize=(7.3, 4.6), sharex=True)
    for r, s in enumerate(('desi', 'sdss')):
        w, edges, nb, acc = plot_lum.binned(f'{s}_0.5x0.1_gspphot', s, 1.0, True)
        for j, g in enumerate(all_):
            q = ax[r, j]; style(q); ks = sorted(k for k in acc if k[0] == g); lo, hi = np.inf, -np.inf
            if not ks: q.axis('off'); continue
            q.set_title(f'{s.upper()} {g}', fontsize=7, loc='left', color=INK, pad=2)
            for k in ks:
                sf, sn = acc[k][:2]; y = plot_lum.smooth(np.where(sn > 0, sf / np.where(sn > 0, sn, 1), np.nan))
                with np.errstate(all='ignore'): m = -2.5 * np.log10(y * w**2 / plot_lum.C_AA) - 48.60
                q.plot(w, m, color=cm(norm(0.5 * (edges[k[1]] + edges[k[1] + 1]))), lw=0.4); f = m[np.isfinite(m)]; lo, hi = min(lo, np.percentile(f, 0.2)), max(hi, np.percentile(f, 99.9))
            q.set_ylim(hi + 0.3, lo - 0.3); q.set_xlim(3550, 9900); q.locator_params(axis='x', nbins=3)
        ax[r, 0].set_ylabel('M_AB(λ)', fontsize=7, color=MUTED)
    for j in range(5): ax[1, j].set_xlabel('rest wavelength [Å]', fontsize=6.5, color=MUTED); ax[1, j].tick_params(labelbottom=True)
    fig.tight_layout(pad=0.4, w_pad=0.3, h_pad=0.8, rect=(0, 0, 0.93, 1)); cax = fig.add_axes([0.945, 0.12, 0.012, 0.76]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=norm), cax=cax); cb.set_label('M_G0 of the bin', fontsize=7, color=MUTED); cb.ax.tick_params(labelsize=6, colors=MUTED, length=2); cb.outline.set_visible(False)
    f = os.path.join(OUT, 'plots_star_lum_tell.png'); fig.savefig(f, dpi=220); plt.close(fig); print(f)


def eigen(s):
    cs = CLASSES[s]; fig, ax = plt.subplots(10, len(cs), figsize=(7.3, 8.6), sharex=True); band = [(n[0], n[-1]) for b in (['desi'] if s == 'desi' else ['dr8', 'dr17']) for n in tellapply.basis(b).bands]
    for j, c in enumerate(cs):
        p = os.path.join(REL, f'star_eigenspectra_v4_{s}_{c}.fits'); M, E, V, P = (fitsio.read(p, ext=e).astype(float) for e in ('MEAN', 'EIGVEC', 'VARIANCE', 'PIXMASK')); A = np.where(P > 0, np.vstack([M, E[:9]]), np.nan); r2 = np.diff(np.r_[0, V]); ok = np.isfinite(A[0]); ax[0, j].set_title(f'{s.upper()} {c}', fontsize=7.5, loc='left', color=INK, pad=3)
        for k in range(10):
            q = ax[k, j]; style(q, 5)
            for lo, hi in band: q.axvspan(lo, hi, color=BAND, lw=0, zorder=0)
            if k: q.axhline(0, color=GRAY, lw=0.3, zorder=1)
            q.plot(WAVE, A[k], color=BLUE if k else INK, lw=0.25, zorder=2); lo, hi = np.nanpercentile(A[k], [0.2, 99.8]); q.set_ylim(lo - 0.12 * (hi - lo), hi + 0.25 * (hi - lo)); q.set_xlim(WAVE[ok][0], WAVE[ok][-1]); q.locator_params(axis='y', nbins=3); q.locator_params(axis='x', nbins=3)
            q.text(0.985, 0.96, 'e0' if k == 0 else f'e{k}  +{r2[k - 1]:.4f}', transform=q.transAxes, fontsize=5, color=INK, ha='right', va='top')
        ax[-1, j].set_xlabel('rest wavelength [Å]', fontsize=6, color=MUTED)
    fig.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.15); f = os.path.join(OUT, f'plots_star_empca_v4_{s}.png'); fig.savefig(f, dpi=220); plt.close(fig); print(f)


if __name__ == '__main__':
    lum()
    for s in CLASSES: eigen(s)
