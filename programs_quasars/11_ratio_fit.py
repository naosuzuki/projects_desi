#!/usr/bin/env python3
"""Smooth fit of log window-flux ratio vs M_AB(3025), DESI only, all overlap redshifts pooled (Nao 2026-10-31: the redshift
evolution is small). Same per-quasar ratios as 10_ratio_vs_mag.py (COV >= 0.8 in both windows, NORM_SRC 0).
Quasars whose ID appears more than once in the norm table are dropped (float-TARGETID collisions, see HANDOVER_quasars).
Fit to the 0.25-mag binned medians, weighted by their standard error (1.2533 sigma/sqrt n):
  R1 = F(3025)/F(5115): softplus, two straight lines joined smoothly at M0 over width w:
       logR = a + s_bright (M - M0) + (s_faint - s_bright) w ln(1 + exp((M - M0)/w))
       (nearly flat when bright; host-galaxy light at 5115 A makes it fall steeply when faint)
  R2 = F(3025)/F(1282): straight line  logR = a + b (M + 25)
Outputs: checks/qso_ratio_fit.png, checks/qso_ratio_fit.json (parameters, covariance, binned medians and residuals).
"""
import os, json, numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

P = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{P}/checks'
COVER = 0.8; DM = 0.25; NMIN = 100; C_AA = 2.99792458e18
t = pd.read_parquet(f'{P}/v2/qso_norm_table.parquet')
t = t[(t.STORE_SURVEY == 'desi') & (t.NORM_SRC == 0) & t.Z.notna()]
t = t[~t.ID.duplicated(keep=False)]


def softplus(M, a, s_faint, s_bright, M0, w):
    return a + s_bright * (M - M0) + (s_faint - s_bright) * w * np.logaddexp(0, (M - M0) / w)


def line(M, a, b):
    return a + b * (M + 25)


FITS = {'R1': ('3025', '5115', softplus, [0.4, -0.35, -0.02, -21, 1]), 'R2': ('3025', '1282', line, [-0.54, 0.02])}
res = {}
fig, ax = plt.subplots(2, 2, figsize=(14, 9), gridspec_kw=dict(height_ratios=[3, 1]), sharex='col')
for i, (k, (wn, wd, fn, p0)) in enumerate(FITS.items()):
    ok = (t['COV' + wn] >= COVER) & (t['COV' + wd] >= COVER) & (t['W' + wn] > 0) & (t['W' + wd] > 0)
    s = t[ok].copy(); s['logR'] = np.log10(s['W' + wn] / s['W' + wd])
    lf = np.where((s.W3025 - 1).abs() < 1e-3, s.LOGFABS_NORM, s.LOGFABS_NORM + s.logR)
    s['M'] = -2.5 * (lf + 2 * np.log10(3025.0) - np.log10(C_AA)) - 48.60
    e = np.arange(np.floor(s.M.quantile(0.002) / DM) * DM, s.M.quantile(0.998) + DM, DM)
    b = np.digitize(s.M, e) - 1
    rows = []
    for bi in range(len(e) - 1):
        y = s.logR[b == bi].to_numpy()
        if len(y) >= NMIN:
            q16, q50, q84 = np.percentile(y, [16, 50, 84])
            rows.append((0.5 * (e[bi] + e[bi + 1]), q50, q16, q84, len(y), 1.2533 * (q84 - q16) / 2 / np.sqrt(len(y))))
    B = pd.DataFrame(rows, columns=['M', 'med', 'p16', 'p84', 'n', 'err'])
    p, C = curve_fit(fn, B.M, B.med, p0=p0, sigma=B.err, absolute_sigma=True, maxfev=20000)
    B['fit'] = fn(B.M, *p); B['resid'] = B.med - B.fit
    chi2 = float(((B.resid / B.err)**2).sum()); dof = len(B) - len(p)
    names = ['a', 's_faint', 's_bright', 'M0', 'w'] if k == 'R1' else ['a', 'b']
    res[k] = dict(model='softplus: a + s_bright(M-M0) + (s_faint-s_bright) w ln(1+exp((M-M0)/w))' if k == 'R1' else 'line: a + b*(M+25)', n_quasars=int(len(s)),
                  z_range=[round(float(s.Z.min()), 3), round(float(s.Z.max()), 3)], params=dict(zip(names, map(float, p))),
                  errors=dict(zip(names, map(float, np.sqrt(np.diag(C))))), cov=C.tolist(), chi2=chi2, dof=dof,
                  rms_resid_dex=float(np.sqrt(np.average(B.resid**2, weights=B.n))), max_abs_resid_dex=float(B.resid.abs().max()),
                  M_range=[float(B.M.min()), float(B.M.max())], bins=B.round(5).to_dict(orient='list'))
    a0 = ax[0, i]
    h = a0.hexbin(s.M, s.logR, gridsize=(90, 60), extent=(e[0], e[-1], -1.6, 1.0) if k == 'R1' else (e[0], e[-1], -1.1, 0.0),
                  bins='log', cmap='Greys', mincnt=1)
    a0.fill_between(B.M, B.p16, B.p84, color='C0', alpha=0.15, label='16-84 %')
    a0.errorbar(B.M, B.med, B.err, fmt='o', ms=3, color='C0', label=f'DESI median per {DM} mag ({len(s):,} quasars)')
    mm = np.linspace(B.M.min(), B.M.max(), 300); a0.plot(mm, fn(mm, *p), 'C3-', lw=2, label='fit')
    ptxt = ', '.join(f'{n}={v:.3f}' for n, v in zip(names, p))
    a0.set_title(('R1 = F(3000-3050)/F(5100-5130)' if k == 'R1' else 'R2 = F(3000-3050)/F(1275-1290)') +
                 f'   z {res[k]["z_range"][0]}-{res[k]["z_range"][1]}\n{res[k]["model"]}\n{ptxt}', fontsize=9)
    a0.set_ylabel(f'log {k}'); a0.legend(fontsize=8, loc='lower left' if k == 'R1' else 'upper right'); a0.invert_xaxis() if i == 0 else None
    a1 = ax[1, i]
    a1.axhline(0, color='k', lw=0.5); a1.errorbar(B.M, B.resid, B.err, fmt='o', ms=3, color='C0')
    a1.set_ylim(-0.12, 0.12); a1.set_xlabel('M_AB at 3025 A (absolute)'); a1.set_ylabel('median - fit (dex)')
    a1.set_title(f'rms {res[k]["rms_resid_dex"]:.4f} dex, max {res[k]["max_abs_resid_dex"]:.4f} dex', fontsize=9)
ax[0, 1].invert_xaxis()
fig.tight_layout(); fig.savefig(f'{OUT}/qso_ratio_fit.png', dpi=90)
json.dump(res, open(f'{OUT}/qso_ratio_fit.json', 'w'), indent=1)
for k, r in res.items():
    print(k, r['model'], r['n_quasars'], {n: f'{v:.4f}+-{r["errors"][n]:.4f}' for n, v in r['params'].items()},
          f'chi2/dof {r["chi2"]:.1f}/{r["dof"]}', f'rms {r["rms_resid_dex"]:.4f} max {r["max_abs_resid_dex"]:.4f}', r['M_range'])
