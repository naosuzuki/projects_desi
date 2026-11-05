#!/usr/bin/env python3
"""Luminosity-binned quasar spectra from the v3 stacks (Nao 2026-10-31; same layout as the galaxy thread's
programs_templates/21_plot_lum.py). v3 has no redshift cells: every member is normalized at 3000-3050 A (measured where the window is
covered, own window x the DESI-fitted ratio for z<0.19 and z>2.25; 04_cells.py --chain, 05_stack.py --v3), so a stack joins all redshifts.
One column per class (--by ew: primary-line EW class from EBIN, pairs merged: weak < 1/2, typical 1/2-2, strong > 2 x the
line's median; Halpha / Hbeta / Mg II / C IV by redshift. --by slope: continuum slope alpha_lambda from SBIN).
Per 1-mag bin of the stack's median M_AB(3025) (LOGFABS_3025 of its members) the member-weighted mean over all stacks of that
class and luminosity. Per pixel a stack enters only where >= max(NMIN, FRAC x its best-covered count) members cover it, with
its weight (members per pixel) tapered linearly to 0 over the last TAPER pixels at each end of that range; 9-pixel boxcar.
Shown in F_lambda (Nao 2026-10-31: a quasar brightens to the blue in F_lambda, while in F_nu / AB magnitude the same
spectrum fades to the blue): --y logflam (default) is log10 F_lambda,abs [erg/s/cm2/A at 10 pc], anchored at 3000-3050 A to
the bin's member-weighted mean log F_lambda,abs(3025) (from M_AB(3025)); --y mag is the same in the ST magnitude
M_ST = -2.5 log10 F_lambda,abs - 21.10; --y ab is the absolute AB magnitude of F_nu = F_lambda lambda^2.
Pixels with fewer than --min-total (10 000) member spectra over all stacks are dropped (the sparse blue and red edges).
Bottom row: quasars contributing per pixel. --simple: one panel, all classes together, one curve per 1-mag bin.
Stacks are the ivar-weighted means WMEAN (what the EM-PCA fits); the unweighted MEAN carries single-member spikes.
Usage: python3 12_lum_stacks.py [--stacks quasar/v3/qso_stacks.npz] [--by ew|slope] [--simple] [--y logflam|mag|ab] [--dm 1.0] [--out PNG]
Outputs: PNG (--out; the launcher puts it in plots_luminosity/) and CSV (--csv; csv_files/), one row per class x magnitude bin.
"""
import os, argparse
import numpy as np
import qsolib as Q
import pandas as pd
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{P}/checks'
NMIN = 30; FRAC = 0.1; TAPER = 200; MIN_TOTAL = 10000; SMOOTH = 9; C_AA = 2.99792458e18; REF = (3000.0, 3050.0); MRANGE = (-30.0, -15.0)
GROUPS = {'ew': ('STACK_EBIN', {0: 'weak lines  EW < 1/2 median', 1: 'typical  1/2-2 x median', 2: 'strong lines  EW > 2 x median',
                                -1: 'line EW not measured'},
                 lambda v: -1 if v < 0 else min(2, max(0, (int(v) - 1) // 2))),   # EBIN 0,1 / 2,3 / 4,5
          'slope': ('STACK_SBIN', {0: 'alpha_lam < -1.75', 1: '-1.75 .. -1.25', 2: '> -1.25', -1: 'slope unknown'},
                    lambda v: -1 if v < 0 else min(2, int(v) // 2)),             # SBIN 0,1 / 2,3 / 4,5
          'bal': ('STACK_BBIN', {0: 'non-BAL', 1: 'BAL', -1: 'BAL class collapsed'},
                  lambda v: -1 if v < 0 else min(1, int(v)))}                   # BBIN 0 / 1-3 (two-class tables: 0 / 1; Nao 2026-11-04)
LINES = {'Lya': 1215.67, 'CIV': 1549.06, 'CIII]': 1908.73, 'MgII': 2798.75, '[OII]': 3728.4, 'Hb': 4862.7, '[OIII]': 5008.2, 'Ha': 6564.6}


def mab(logf, lam=3025.0):
    return -2.5 * (logf + 2 * np.log10(lam) - np.log10(C_AA)) - 48.60


def smooth(y):
    ok = np.isfinite(y); num = uniform_filter1d(np.where(ok, y, 0.0), SMOOTH); den = uniform_filter1d(ok.astype(float), SMOOTH)
    return np.where(ok & (den > 0.5), num / np.where(den > 0, den, 1), np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', default=f'{P}/v3/qso_stacks.npz'); ap.add_argument('--by', default='ew', choices=list(GROUPS))
    ap.add_argument('--dm', type=float, default=1.0); ap.add_argument('--out', default=None)
    ap.add_argument('--csv', default=None, help='table path (default: next to the PNG; Nao 2026-11-02: PNGs in plots_luminosity/, CSVs in csv_files/)')
    ap.add_argument('--simple', action='store_true', help='one panel, all classes together')
    ap.add_argument('--min-total', type=int, default=MIN_TOTAL, help='drop the blue and red edges: pixels with fewer member spectra in all (over the stacks with >= NMIN members)')
    ap.add_argument('--cells', default=f'{P}/qso_cells_v3_all.parquet')
    ap.add_argument('--snmin', type=float, default=5.0, help='per stack, drop the blue end until the running median S/N per pixel (100 px) reaches this; 0 = off')
    ap.add_argument('--drop', type=float, default=0.9, help='... and until the running mean flux (100 px) is >= this x the 200 px redward')
    ap.add_argument('--max-mag', type=float, default=-19.0, help='drop luminosity bins fainter than this M_AB(3025) (Nao 2026-11-01: the faintest bins are host galaxies, not quasars; -19 since v5, Nao 2026-11-02: the -19..-18 bin still dips at 3000 A)')
    ap.add_argument('--y', default='logflam', choices=['logflam', 'mag', 'ab'], help='y axis: log F_lambda,abs, ST magnitude of F_lambda, or AB magnitude of F_nu')
    a = ap.parse_args(); out = a.out or f'{OUT}/qso_lum_stacks_v3_{"simple" if a.simple else a.by}.png'
    gcol, glab, gfun = GROUPS[a.by]; cols = [g for g in [0, 1, 2, -1] if g in glab]
    if a.simple:
        glab = {0: 'all quasars'}; gfun = lambda v: 0; cols = [0]
    s = np.load(a.stacks, allow_pickle=True); w = s['WAVE']; F = s['WMEAN'].astype(float); N = s['N'].astype(float); cells = s['STACK_CELL']
    E = s['WERR'].astype(float)
    # WMEAN (ivar-weighted), as the EM-PCA uses: the unweighted MEAN carries single-member spikes (e.g. 7173 A rest in L1_E2_S5_O2_B-1,
    # one DESI member with a ~1e7 flux pixel near 9650 A observed and ~zero ivar)
    # per stack: median M_AB(3025) and class of its members
    c = pd.read_parquet(a.cells, columns=['TRAIN', 'STACK_CELL', 'LOGFABS_3025', 'STACK_EBIN', 'STACK_SBIN', 'STACK_BBIN'])
    c = c[c.TRAIN & (c.STACK_CELL != '')]
    info = c.groupby('STACK_CELL').agg(M=('LOGFABS_3025', lambda x: mab(np.median(x))), G=(gcol, 'first'))
    edges = np.arange(MRANGE[0], a.max_mag + 1e-9, a.dm); nb = len(edges) - 1
    ref = (w >= REF[0]) & (w <= REF[1]); j = np.arange(len(w)); acc = {}
    big = N.max(axis=1) >= NMIN; edge = (N[big] * (N[big] >= np.maximum(NMIN, FRAC * N[big].max(axis=1))[:, None])).sum(axis=0) < a.min_total
    F[:, edge] = np.nan; N[:, edge] = 0                   # Nao 2026-11-01: the sparse blue and red edges are dropped
    print(f'edges with < {a.min_total:,} spectra per pixel dropped: {w[~edge].min():.0f}-{w[~edge].max():.0f} A')
    goodall = (N >= np.maximum(NMIN, FRAC * N.max(axis=1))[:, None]) & np.isfinite(F)
    if a.snmin > 0:
        trimmed = Q.blue_trim(w, F, E, goodall, a.snmin, a.drop)
        cut = [w[goodall[i]].min() for i in range(len(F)) if goodall[i].any() and trimmed[i].sum() < goodall[i].sum()]
        print(f'blue ends trimmed (S/N < {a.snmin:g} or flux falling below {a.drop:g} x the 200 px redward) in {len(cut)} of {int((N.max(axis=1) >= NMIN).sum())} stacks', flush=True)
        goodall = trimmed
    for i, cell in enumerate(cells):
        if cell not in info.index or N[i].max() < NMIN:
            continue
        Mi, g = info.at[cell, 'M'], gfun(info.at[cell, 'G']); b = int(np.digitize(Mi, edges) - 1)
        if not (0 <= b < nb):
            continue
        good = goodall[i]
        if not good.any():
            continue
        g0, g1 = j[good].min(), j[good].max()
        wt = N[i] * np.clip(np.minimum(j - g0, g1 - j) / TAPER, 0, 1) * good
        k = (g, b)
        if k not in acc:
            acc[k] = [np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0.0, 0.0]
        acc[k][0] += np.where(good, wt * F[i], 0.0); acc[k][1] += wt; acc[k][2] += np.where(good, N[i], 0.0)
        nref = N[i][ref].mean(); acc[k][3] += nref * Mi; acc[k][4] += nref
    fig, ax = plt.subplots(2, len(cols), figsize=(13 if a.simple else 5.4 * len(cols), 10.5), sharex=True, squeeze=False,
                           gridspec_kw=dict(height_ratios=[3, 1]))
    cm = plt.get_cmap('rainbow_r'); norm = plt.Normalize(-29, -17); rows = []
    for (g, b), (sf, swt, sn, sm_, sN) in sorted(acc.items()):
        col = cols.index(g); mc = 0.5 * (edges[b] + edges[b + 1])
        y = smooth(np.where(swt > 0, sf / np.where(swt > 0, swt, 1), np.nan)); fnu = y * w**2
        r = np.nanmean(fnu[ref])
        if sN <= 0:
            continue
        mg = sm_ / sN
        logf3025 = -(mg + 48.60) / 2.5 - 2 * np.log10(3025.0) + np.log10(C_AA)   # mean M_AB(3025) -> log F_lambda,abs(3025)
        with np.errstate(all='ignore'):
            lflam = logf3025 + np.log10(y)      # every member is normalized to unit mean at 3000-3050 A (its template's), so the stack is in units of F(3025)
            m = {'logflam': lflam, 'mag': -2.5 * lflam - 21.10, 'ab': mg - 2.5 * np.log10(fnu / (r if np.isfinite(r) and r > 0 else np.nanmean(fnu)))}[a.y]
        ax[0, col].plot(w, m, color=cm(norm(mc)), lw=0.8, label=f'M {edges[b]:.0f} to {edges[b+1]:.0f}  (<M> {mg:.2f}, max {int(sn.max()):,})')
        ax[1, col].plot(w, np.where(sn > 0, sn, np.nan), color=cm(norm(mc)), lw=0.8)
        rows.append(dict(GROUP=g, LABEL=glab[g], M_lo=edges[b], M_hi=edges[b + 1], M3025_mean=round(mg, 3), max_qso_per_pixel=int(sn.max()),
                         wave_lo=float(w[sn > 0].min()), wave_hi=float(w[sn > 0].max())))
    for col, g in enumerate(cols):
        a0 = ax[0, col]; a0.grid(alpha=0.3)
        a0.set_xlim(700, 9900) if a.simple else (a0.set_xlim(600, 10000), a0.set_xscale('log'))
        a0.set_ylim(*{'logflam': (-2.6, 4.6), 'mag': (-15, -32), 'ab': (-14, -31)}[a.y])
        a0.axvspan(*REF, color='C0', alpha=0.25); a0.axvspan(5100, 5130, color='C3', alpha=0.2); a0.axvspan(1275, 1290, color='C2', alpha=0.2)
        for nm, lam in LINES.items():
            a0.text(lam, a0.get_ylim()[1], nm, fontsize=7, rotation=90, color='0.4', ha='center', va='top')
        a0.set_title(glab[g], fontsize=11); a0.legend(fontsize=7, loc='lower right')
        ax[1, col].set_yscale('log'); ax[1, col].set_xlabel('rest wavelength [A]')
        if not a.simple:
            ax[1, col].set_xticks([700, 1000, 1500, 2000, 3000, 5000, 7000, 10000]); ax[1, col].set_xticklabels(['700', '1000', '1500', '2000', '3000', '5000', '7000', '10000'])
    ax[0, 0].set_ylabel({'logflam': 'log F_lambda,abs  [erg/s/cm2/A at 10 pc]', 'mag': 'absolute F_lambda magnitude  M_ST = -2.5 log F_lambda,abs - 21.1',
                         'ab': 'absolute AB magnitude M_AB(lambda)  (F_nu)'}[a.y] + '  (anchored at 3000-3050 A)')
    ax[1, 0].set_ylabel('quasars per pixel')
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm); fig.colorbar(sm, ax=ax, label='M_AB(3025) (stack median)', fraction=0.015)
    fig.suptitle(f'Quasar spectra by luminosity ({a.dm:.0f}-mag bins of M_AB(3025))' + ('' if a.simple else ' and ' +
                 ('primary-line strength (Halpha / Hbeta / Mg II / C IV by z)' if a.by == 'ew' else 'continuum slope')) +
                 '\nv3 stacks: no redshift cells, every member normalized by the amplitude of its stack template (unity at 3000-3050 A, blue band); blue ends trimmed where S/N < ' + f'{a.snmin:g} or the flux falls toward the blue' +
                 '; member-weighted, edges tapered; red / green = the other normalization windows', fontsize=11)
    fig.savefig(out, dpi=75); plt.close(fig)
    pd.DataFrame(rows).to_csv(a.csv or os.path.splitext(out)[0] + '.csv', index=False); print(out, a.csv or '')
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
