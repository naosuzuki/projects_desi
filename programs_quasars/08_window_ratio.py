#!/usr/bin/env python3
"""Window flux ratios for a continuous normalization (Nao 2026-10-31): R1 = F(3000-3050)/F(5100-5130) and
R2 = F(3000-3050)/F(1275-1290), measured on the stacks that cover both windows (>= COVER of the stack's members at
each window; every member is 1 in its own window, so the stack mean at the other window is the ratio).
Tests the dependence on luminosity (log F_lam,abs at 3025 A, 10 pc), continuum slope bin and z; fits
log R = a + b (logF3025 - X0) per ratio, n-weighted. Output: checks/qso_window_ratio.csv (one row per stack used),
checks/qso_window_ratio_fit.json, checks/qso_window_ratio.png.
"""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{P}/checks'; COVER = 0.8; LOG4PI = 40.078; X0 = 2.0
d = np.load(f'{P}/qso_stacks.npz', allow_pickle=True); W = d['WAVE']; M = d['MEAN']; N = d['N']
t = pd.read_csv(f'{P}/qso_stacks_cells.csv'); t = t[t.ok].copy()
for nm, (lo, hi) in {'5115': (5100, 5130), '3025': (3000, 3050), '1282': (1275, 1290)}.items():
    s = (W >= lo) & (W <= hi); n = N[:, s].astype(float); m = np.nan_to_num(M[:, s])
    t['M' + nm] = ((m * n).sum(1) / np.maximum(n.sum(1), 1))[t.CELL_I]; t['F' + nm] = (n.mean(1))[t.CELL_I] / t.N_MAX
z = t.z_med.to_numpy(); lf = t.logL_med - LOG4PI
rows = []
for r in t.itertuples():
    if r.z_med < 0.8 and r.F3025 >= COVER and r.M3025 > 0:           # unit at 5115
        rows.append(dict(ratio='R1', src='z<0.8 (1 at 5115)', logR=np.log10(r.M3025), logF3025=r.logL_med - LOG4PI + np.log10(r.M3025), cell=r))
    if 0.8 <= r.z_med < 2.0 and r.F5115 >= COVER and r.M5115 > 0:     # unit at 3025
        rows.append(dict(ratio='R1', src='0.8-2.0 (1 at 3025)', logR=-np.log10(r.M5115), logF3025=r.logL_med - LOG4PI, cell=r))
    if 0.8 <= r.z_med < 2.0 and r.F1282 >= COVER and r.M1282 > 0:
        rows.append(dict(ratio='R2', src='0.8-2.0 (1 at 3025)', logR=-np.log10(r.M1282), logF3025=r.logL_med - LOG4PI, cell=r))
    if r.z_med >= 2.0 and r.F3025 >= COVER and r.M3025 > 0:           # unit at 1282
        rows.append(dict(ratio='R2', src='z>=2.0 (1 at 1282)', logR=np.log10(r.M3025), logF3025=r.logL_med - LOG4PI + np.log10(r.M3025), cell=r))
R = pd.DataFrame([{**{k: v for k, v in x.items() if k != 'cell'}, 'STACK_CELL': x['cell'].STACK_CELL, 'n': x['cell'].n,
                   'z': x['cell'].z_med, 'ZBIN': x['cell'].ZBIN, 'LBIN': x['cell'].LBIN, 'SBIN': x['cell'].SBIN,
                   'EBIN': x['cell'].EBIN, 'alpha_med': x['cell'].alpha_med} for x in rows])
R.to_csv(f'{OUT}/qso_window_ratio.csv', index=False)
fit = {}
for k in ('R1', 'R2'):
    s = R[R.ratio == k]; w = s.n.to_numpy(float); x = s.logF3025.to_numpy() - X0; y = s.logR.to_numpy()
    mean = np.average(y, weights=w); sd = np.sqrt(np.average((y - mean)**2, weights=w))
    A = np.c_[np.ones_like(x), x]; sw = np.sqrt(w)
    (a, b), *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    res = y - (a + b * x); sdr = np.sqrt(np.average(res**2, weights=w))
    # with slope (alpha) and z as well
    A3 = np.c_[np.ones_like(x), x, s.alpha_med.to_numpy() + 1.5, s.z.to_numpy() - np.average(s.z, weights=w)]
    c3, *_ = np.linalg.lstsq(A3 * sw[:, None], y * sw, rcond=None); sd3 = np.sqrt(np.average((y - A3 @ c3)**2, weights=w))
    byL = s.groupby('LBIN').apply(lambda g: pd.Series(dict(n=int(g.n.sum()), logR=np.average(g.logR, weights=g.n)))).round(4)
    byS = s.groupby('SBIN').apply(lambda g: pd.Series(dict(n=int(g.n.sum()), logR=np.average(g.logR, weights=g.n)))).round(4)
    bysrc = s.groupby('src').apply(lambda g: pd.Series(dict(nstack=len(g), n=int(g.n.sum()), logR=np.average(g.logR, weights=g.n)))).round(4)
    fit[k] = dict(nstack=len(s), nqso=int(w.sum()), logR_mean=mean, R_mean=10**mean, scatter_dex=sd,
                  fit_L=dict(a=a, b_per_dex=b, X0=X0, scatter_dex=sdr),
                  fit_L_alpha_z=dict(a=c3[0], b_logF=c3[1], c_alpha_plus1p5=c3[2], d_z=c3[3], scatter_dex=sd3),
                  by_LBIN=byL.reset_index().to_dict('records'), by_SBIN=byS.reset_index().to_dict('records'),
                  by_source=bysrc.reset_index().to_dict('records'))
    print(f'{k}: {len(s)} stacks, {int(w.sum()):,} quasars; log R = {mean:+.4f} (R = {10**mean:.4f}), scatter {sd:.4f} dex')
    print(f'    vs luminosity: log R = {a:+.4f} {b:+.4f} (logF3025 - {X0}), scatter {sdr:.4f}')
    print(f'    + slope, z:    a {c3[0]:+.4f}  b(logF) {c3[1]:+.4f}  c(alpha+1.5) {c3[2]:+.4f}  d(z) {c3[3]:+.4f}  scatter {sd3:.4f}')
    print('    by source:\n' + bysrc.to_string()); print('    by LBIN:\n' + byL.T.to_string()); print('    by SBIN:\n' + byS.T.to_string())
json.dump(fit, open(f'{OUT}/qso_window_ratio_fit.json', 'w'), indent=1, default=float)
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
for i, k in enumerate(('R1', 'R2')):
    s = R[R.ratio == k]
    for j, (col, lab) in enumerate([('logF3025', 'log F_lam,abs(3025) at 10 pc'), ('alpha_med', 'median alpha_lambda of the stack'), ('z', 'z')]):
        a_ = ax[i, j]; sc = a_.scatter(s[col], s.logR, s=np.sqrt(s.n), c=s.SBIN, cmap='coolwarm', alpha=0.7)
        if col == 'logF3025':
            xx = np.linspace(s[col].min(), s[col].max(), 10); f = fit[k]['fit_L']; a_.plot(xx, f['a'] + f['b_per_dex'] * (xx - X0), 'k-')
        a_.axhline(fit[k]['logR_mean'], color='gray', ls=':'); a_.set_xlabel(lab); a_.set_ylabel(f'log {k}')
    ax[i, 0].set_title({'R1': 'R1 = F(3000-3050)/F(5100-5130)', 'R2': 'R2 = F(3000-3050)/F(1275-1290)'}[k])
fig.colorbar(sc, ax=ax.ravel().tolist(), label='SBIN (-1 = slope collapsed)'); fig.savefig(f'{OUT}/qso_window_ratio.png', dpi=80)
