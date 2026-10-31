#!/usr/bin/env python3
"""Luminosity-binned stellar spectra from the step-6 HR-cell stacks (Nao 2026-10-31; same layout as the galaxy thread's
programs_templates/21_plot_lum.py and the quasar thread's programs_quasars/12_lum_stacks.py).
One column per temperature class (step7_empca.star_class of the stack's median M_G0, (BP-RP)_0: WD, HOT, AF, GK, M). Per 1-mag
bin of the stack's median dereddened absolute Gaia G (MG0_MED, Vega) the star-weighted mean over all stacks of that class and
luminosity (all colours and [M/H] splits together). Per pixel a stack enters only where NSPEC >= 0.5 N (as step7
--min-member-frac 0.5), weighted by NSPEC; 9-pixel boxcar. SDSS is cut to 3810-9200 A, the range of the SDSS v2 bases (the
stacks jump at the edges of the SDSS-I / BOSS member subsets, HANDOVER_stars).
No anchoring is needed: the stacks are already dereddened absolute F_lambda at 10 pc, so
M_AB(lambda) = -2.5 log10(F_lambda lambda^2 / c) - 48.60 directly. Bottom row: stars contributing per pixel.
--flam (Nao 2026-10-31; as the quasar thread's 12_lum_stacks.py --simple): y = the F_lambda magnitude (ST system,
M_ST = -2.5 log10 F_lambda,abs - 21.10) instead of the F_nu (AB) magnitude, so hot stars brighten to the blue.
Usage: python3 plot_lum.py --survey desi|sdss [--flam] [--dm 1.0] [--out PNG]
Outputs: PNG (default plots_stars_by_lum_<survey>[_flam].png in the repo root) and .csv (one row per class x magnitude bin).
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import fitsio
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from step7_empca import star_class, CLASSES
import starlib

V2 = os.path.expanduser('~/data/sdssdesi_stars/v2')
SMOOTH = 9; MINFRAC = 0.5; WRANGE = {'desi': (0, 1e9), 'sdss': (3810, 9200)}; C_AA = 2.99792458e18; MRANGE = (2.0, 15.0)
CLASS_LAB = {'WD': 'WD  M_G0 > 9, (BP-RP)_0 < 0.8', 'HOT': 'HOT  (BP-RP)_0 < 0.3', 'AF': 'AF  (BP-RP)_0 0.3-0.7',
             'GK': 'GK  (BP-RP)_0 0.7-2.0', 'M': 'M  (BP-RP)_0 > 2.0'}
LINES = {'Ca K': 3934.8, 'Ca H': 3969.6, 'Hd': 4102.9, 'G': 4305.0, 'Hg': 4341.7, 'Hb': 4862.7, 'Mg b': 5175.0, 'Na D': 5895.0,
         'Ha': 6564.6, 'TiO': 7054.0, 'Ca T': 8542.1}


def smooth(y):
    ok = np.isfinite(y); num = uniform_filter1d(np.where(ok, y, 0.0), SMOOTH); den = uniform_filter1d(ok.astype(float), SMOOTH)
    return np.where(ok & (den > 0.5), num / np.where(den > 0, den, 1), np.nan)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', default='desi', choices=['desi', 'sdss'])
    ap.add_argument('--tag', default=None); ap.add_argument('--dm', type=float, default=1.0); ap.add_argument('--out', default=None)
    ap.add_argument('--flam', action='store_true', help='y = F_lambda (ST) magnitude instead of F_nu (AB)')
    a = ap.parse_args(); tag = a.tag or f'{a.survey}_0.5x0.1_gspphot'
    out = a.out or os.path.join(os.path.dirname(HERE), f'plots_stars_by_lum_{a.survey}{"_flam" if a.flam else ""}.png')
    F = fitsio.FITS(os.path.join(V2, f'stacks_{tag}.fits'))
    fl = F['FLUX'].read().astype(float); ns = F['NSPEC'].read().astype(float); arr = F['CELLS'].read(); F.close()
    c = pd.DataFrame({k: arr[k].astype(arr[k].dtype.newbyteorder('=')) for k in ('N', 'MG0_MED', 'COL0_MED')})
    c['CLS'] = [star_class(m, k) for m, k in zip(c.MG0_MED, c.COL0_MED)]
    w = starlib.SWAVE
    edges = np.arange(MRANGE[0], MRANGE[1] + 1e-9, a.dm); nb = len(edges) - 1
    acc = {}                                               # (class, mbin) -> [sum n f, sum n, sum N M, sum N, nstacks]
    for i, r in c.iterrows():
        b = int(np.digitize(r.MG0_MED, edges) - 1)
        if not (0 <= b < nb):
            continue
        good = (ns[i] >= MINFRAC * r.N) & np.isfinite(fl[i]) & (fl[i] > 0) & (w >= WRANGE[a.survey][0]) & (w <= WRANGE[a.survey][1])
        if not good.any():
            continue
        k = (r.CLS, b)
        if k not in acc:
            acc[k] = [np.zeros(len(w)), np.zeros(len(w)), 0.0, 0.0, 0]
        acc[k][0] += np.where(good, ns[i] * fl[i], 0.0); acc[k][1] += np.where(good, ns[i], 0.0)
        acc[k][2] += r.N * r.MG0_MED; acc[k][3] += r.N; acc[k][4] += 1
    cols = [g for g in CLASSES if any(k[0] == g for k in acc)]
    fig, ax = plt.subplots(2, len(cols), figsize=(5.4 * len(cols), 10.5), sharex=True, squeeze=False,
                           gridspec_kw=dict(height_ratios=[3, 1]))
    cm = plt.get_cmap('rainbow_r'); norm = plt.Normalize(*MRANGE)          # bright red -> faint purple
    rows = []; yr = {g: [np.inf, -np.inf] for g in cols}
    for (g, b), (sf, sn, sm_, sN, nst) in sorted(acc.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        col = cols.index(g); mg = sm_ / sN
        y = smooth(np.where(sn > 0, sf / np.where(sn > 0, sn, 1), np.nan))
        with np.errstate(all='ignore'):
            m = -2.5 * np.log10(y) - 21.10 if a.flam else -2.5 * np.log10(y * w**2 / C_AA) - 48.60
        ax[0, col].plot(w, m, color=cm(norm(0.5 * (edges[b] + edges[b + 1]))), lw=0.8,
                        label=f'M_G {edges[b]:.0f} to {edges[b+1]:.0f}  (<M_G> {mg:.2f}, {nst} stacks, max {int(sn.max()):,})')
        ax[1, col].plot(w, np.where(sn > 0, sn, np.nan), color=cm(norm(0.5 * (edges[b] + edges[b + 1]))), lw=0.8)
        fin = m[np.isfinite(m)]
        if len(fin):
            yr[g] = [min(yr[g][0], np.percentile(fin, 0.2)), max(yr[g][1], np.percentile(fin, 99.9))]
        rows.append(dict(CLASS=g, M_lo=edges[b], M_hi=edges[b + 1], MG0_mean=round(mg, 3), n_stacks=nst, n_stars=int(sN),
                         max_stars_per_pixel=int(sn.max()), wave_lo=float(w[sn > 0].min()), wave_hi=float(w[sn > 0].max())))
    for col, g in enumerate(cols):
        a0 = ax[0, col]; lo, hi = yr[g]; a0.set_ylim(hi + 0.5, lo - 1.2); a0.set_xlim(3550, 9900); a0.grid(alpha=0.3)
        for nm, lam in LINES.items():
            a0.text(lam, lo - 1.1, nm, fontsize=7, rotation=90, color='0.4', ha='center', va='top')
        a0.set_title(CLASS_LAB[g], fontsize=11); a0.legend(fontsize=7, loc='best' if a.flam else 'lower right')
        ax[1, col].set_yscale('log'); ax[1, col].set_xlabel('rest wavelength [A]')
    ax[0, 0].set_ylabel('absolute F_lambda magnitude  M_ST = -2.5 log F_lambda,abs - 21.1  (10 pc, dereddened)' if a.flam
                        else 'absolute AB magnitude M_AB(lambda)  (F_nu at 10 pc, dereddened)'); ax[1, 0].set_ylabel('stars per pixel')
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm); fig.colorbar(sm, ax=ax, label='M_G0 (stack median, Vega)', fraction=0.015)
    fig.suptitle(f'Stellar spectra by luminosity ({a.dm:.0f}-mag bins of M_G0) and temperature class, {a.survey.upper()} v2 HR-cell stacks'
                 '\nstar-weighted mean over all colour / [M/H] stacks of the bin (pixels with >= half the members); '
                 + ('y = absolute F_lambda (ST) magnitude' if a.flam else 'y = absolute AB magnitude') + ' straight from F_lambda at 10 pc (no anchoring)', fontsize=12)
    fig.savefig(out, dpi=75); plt.close(fig)
    pd.DataFrame(rows).to_csv(os.path.splitext(out)[0] + '.csv', index=False); print(out)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
