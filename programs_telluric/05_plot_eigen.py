#!/usr/bin/env python3
"""Stellar eigenspectra e0 to e9 after the telluric correction, one figure per survey and class, with v2 underneath (Nao 2026-11-02).
e0 is the weighted mean (step7_empca.py keeps it as component 0, G flux = 1), e1..e9 are the first nine EIGVEC; v2 vectors are flipped to the sign of the corrected ones. Telluric bands shaded.
--sqrtn (Nao 2026-11-03): the sqrt(N) stack weights of step7_empca.py over the equal weights, both telluric corrected.
Output: plots_luminosity/plots_eigen_v2tell_<survey>_<class>.png, plots_eigen_sqrtn_<survey>_<class>.png
"""
import os, sys
import numpy as np, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import tellapply
TELL = '/Volumes/exdisk3/sdssdesi_stars/v2_tell'; V2 = os.path.expanduser('~/data/sdssdesi_stars/v2'); OUT = os.path.join(ROOT, 'plots_luminosity')
WAVE = 10**(3.5502 + 1e-4 * np.arange(4447)); CLASSES = {'sdss': ['WD', 'AF', 'GK', 'M'], 'desi': ['WD', 'HOT', 'AF', 'GK', 'M']}; NSHOW = 10
BLUE, GRAY, BAND, INK, MUTED = '#2a78d6', '#a8a79f', '#ecebe6', '#0b0b0b', '#52514e'


def read(d, survey, c, sfx=''):
    p = f'{d}/empca_{survey}_0.5x0.1_gspphot_{c}{sfx}.fits'
    return {e: fitsio.read(p, ext=e).astype(float) for e in ('MEAN', 'EIGVEC', 'VARIANCE', 'PIXMASK')}


def vectors(b):
    """e0..e9 with the pixels outside the fit blanked."""
    return np.where(b['PIXMASK'] > 0, np.vstack([b['MEAN'], b['EIGVEC'][:NSHOW - 1]]), np.nan)


def bands_of(survey):
    return [(n[0], n[-1]) for s in (['desi'] if survey == 'desi' else ['dr8', 'dr17']) for n in tellapply.basis(s).bands]


def plot(survey, c, sqrtn=False):
    a, b = (read(TELL, survey, c, '_sqrtn'), read(TELL, survey, c)) if sqrtn else (read(TELL, survey, c), read(V2, survey, c)); labs = ('equal weights', 'sqrt(N) weights') if sqrtn else ('v2', 'telluric corrected'); A, B = vectors(a), vectors(b); r2 = np.diff(np.r_[0, a['VARIANCE']]); w = np.isfinite(A[0])
    fig, ax = plt.subplots(NSHOW, 1, figsize=(12, 15), sharex=True)
    for k in range(NSHOW):
        o = np.isfinite(A[k]) & np.isfinite(B[k]); y = B[k] * (np.sign(np.dot(A[k][o], B[k][o])) or 1); q = ax[k]
        for lo, hi in bands_of(survey): q.axvspan(lo, hi, color=BAND, lw=0, zorder=0)
        q.axhline(0, color=GRAY, lw=0.5, zorder=1); q.plot(WAVE, y, color=GRAY, lw=1.1, zorder=2, label=labs[0]); q.plot(WAVE, A[k], color=BLUE, lw=0.6, zorder=3, label=labs[1])
        lo, hi = np.nanpercentile(A[k], [0.2, 99.8]); pad = 0.12 * (hi - lo); q.set_ylim(lo - pad, hi + pad)     # the y range follows the vector, not the strongest line or edge pixel
        tag = 'e0  mean' if k == 0 else f'e{k}  R² +{r2[k - 1]:.4f}'
        q.set_title(tag, fontsize=8.5, loc='left', color=INK, pad=2)
        q.tick_params(labelsize=7, colors=MUTED, length=2); [q.spines[s].set_visible(False) for s in ('top', 'right')]; [q.spines[s].set_color(GRAY) for s in ('left', 'bottom')]
    ax[0].set_xlim(WAVE[w][0], WAVE[w][-1]); ax[0].legend(loc='lower right', bbox_to_anchor=(1, 1.0), fontsize=8, frameon=False, ncol=2, borderaxespad=0); ax[-1].set_xlabel('wavelength [Å]; shaded: telluric bands', fontsize=9, color=MUTED)
    fig.suptitle(f'{survey.upper()} {c}: eigenspectra e0 to e9, {labs[1]} over {labs[0]}{", telluric corrected" if sqrtn else ""}', fontsize=11, x=0.045, ha='left', color=INK)
    fig.tight_layout(h_pad=0.3); f = os.path.join(OUT, f'plots_eigen_{"sqrtn" if sqrtn else "v2tell"}_{survey}_{c}.png'); fig.savefig(f, dpi=120); plt.close(fig); print(f, flush=True)


if __name__ == '__main__':
    for s, cs in CLASSES.items():
        for c in cs: plot(s, c, '--sqrtn' in sys.argv)
