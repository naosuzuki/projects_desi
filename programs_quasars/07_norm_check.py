#!/usr/bin/env python3
"""Normalization check (Nao 2026-10-30): v1 quasar stacks by luminosity bin, per redshift slice.
Each stack is unit mean in its window (5115 A z<0.8, 3025 A 0.8-2.0, 1282 A >=2.0); LBIN = floor(log10 F_lam,abs/0.4)
measured AT that window, so it is a different wavelength on each side of z=0.8 and z=2.0.
Fig 1 (qso_norm_by_L_shape.png): per ZBIN, the n-weighted stack of each LBIN as stacked (unit window mean).
Fig 2 (qso_norm_by_L_abs.png): the same scaled to absolute level, log F_lam,abs (10 pc) vs rest wavelength, one panel per
LBIN, coloured by z; across a window boundary the curves of one luminosity population should join.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
P = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{P}/checks'
LOG4PI = 40.078
d = np.load(f'{P}/qso_stacks.npz', allow_pickle=True); W = d['WAVE']; M = d['MEAN']; NP = d['N']
t = pd.read_csv(f'{P}/qso_stacks_cells.csv'); t = t[t.ok]
ZE = np.r_[np.arange(0, 1.01, 0.1), np.arange(1.2, 3.01, 0.2), [3.5, 4.0, 4.5, 5.0]]
WIN = lambda z: 5115.0 if z < 0.8 else (3025.0 if z < 2.0 else 1282.5)
rows = []
for (zb, lb), g in t.groupby(['ZBIN', 'LBIN']):
    i = g.CELL_I.to_numpy(); w = np.where(np.isfinite(M[i]), NP[i], 0) * 1.0
    m = (np.nan_to_num(M[i]) * w).sum(0) / np.maximum(w.sum(0), 1); cov = w.sum(0) > 30
    m = np.where(cov, uniform_filter1d(m, 5), np.nan)
    logfa = np.average(g.logL_med - LOG4PI, weights=g.n)
    rows.append(dict(zb=zb, lb=lb, n=g.n.sum(), ncell=len(g), z=np.average(g.z_med, weights=g.n), logfa=logfa, m=m))
R = pd.DataFrame(rows)
print(f'{len(t)} stacks (n>=300), {R.zb.nunique()} z slices, {R.lb.nunique()} luminosity bins, {len(R)} (z, L) combinations')
print(R.pivot(index='lb', columns='zb', values='n').fillna(0).astype(int).to_string())
cm = plt.cm.viridis; lbs = sorted(R.lb.unique()); col = {l: cm(k / max(1, len(lbs) - 1)) for k, l in enumerate(lbs)}
fig, ax = plt.subplots(6, 4, figsize=(20, 22)); ax = ax.ravel()
for zb in range(24):
    a = ax[zb]
    if not (R.zb == zb).any():
        a.set_visible(False); continue
    for r in R[R.zb == zb].itertuples():
        a.plot(W, np.where(r.m > 0, r.m, np.nan), color=col[r.lb], lw=0.7, label=f'L{r.lb} n={r.n}')
    a.axvline(WIN(0.5 * (ZE[zb] + ZE[zb + 1])), color='r', ls=':', lw=0.8)
    a.set_xscale('log'); a.set_xlim(600, 10000); a.set_ylim(0, 6); a.set_title(f'z {ZE[zb]:.1f}-{ZE[zb+1]:.1f}', fontsize=9)
    a.legend(fontsize=6, ncol=2); a.tick_params(labelsize=7)
fig.suptitle('v1 quasar stacks by luminosity bin (as stacked: unit mean in the red-dotted window)'); fig.tight_layout()
fig.savefig(f'{OUT}/qso_norm_by_L_shape.png', dpi=80); plt.close(fig)
fig, ax = plt.subplots(4, 6, figsize=(24, 15)); ax = ax.ravel()
zc = plt.cm.turbo
for k, lb in enumerate([l for l in lbs if l >= -4][:24]):
    a = ax[k]
    for r in R[R.lb == lb].itertuples():
        a.plot(W, np.log10(np.where(r.m > 0, r.m, np.nan)) + r.logfa, color=zc(r.zb / 23), lw=0.8)
    a.set_xscale('log'); a.set_xlim(600, 10000); a.set_title(f'LBIN {lb}: log F_abs(window) {0.4*lb:.1f}-{0.4*lb+0.4:.1f}', fontsize=9); a.tick_params(labelsize=7)
sm = plt.cm.ScalarMappable(cmap=zc, norm=plt.Normalize(0, 5)); fig.colorbar(sm, ax=ax.tolist(), label='z')
fig.suptitle('v1 quasar stacks at absolute level, log10 F_lam,abs (erg/s/cm2/A at 10 pc), per luminosity bin, coloured by z')
fig.savefig(f'{OUT}/qso_norm_by_L_abs.png', dpi=80); plt.close(fig)
R.drop(columns='m').to_csv(f'{OUT}/qso_norm_by_L.csv', index=False)
