#!/usr/bin/env python3
"""Window flux ratios vs absolute magnitude, per quasar (Nao 2026-10-31), from quasar/v2/qso_norm_table.parquet (09_norm_table.py).
R1 = F_lam(3000-3050)/F_lam(5100-5130), R2 = F_lam(3000-3050)/F_lam(1275-1290), measured on each quasar's own normalized
spectrum when it covers both windows (good-pixel fraction >= COVER in each, window-mean normalization NORM_SRC 0).
Magnitude axis: M_AB at 3025 A = -2.5 log10(F_nu,abs) - 48.60, F_nu = F_lam * lam^2 / c, F_lam,abs at 10 pc = 10**LOGFABS_NORM
moved to 3025 A with the measured ratio. Bins of 0.5 mag where every bin of the populated range holds >= NMIN quasars, else 1.0.
Outputs: checks/qso_ratio_vs_mag.png, checks/qso_ratio_vs_mag.csv (median, 16/84 %, n per survey x z range x mag bin).
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

P = os.path.expanduser('~/data/sdssdesi_products/quasar'); OUT = f'{P}/checks'
COVER = 0.8; NMIN = 100; C_AA = 2.99792458e18
t = pd.read_parquet(f'{P}/v2/qso_norm_table.parquet')
t = t[(t.NORM_SRC == 0) & t.Z.notna()]
t['SV'] = np.where(t.STORE_SURVEY == 'desi', 'DESI', np.where(t.STORE_SURVEY == 'dr17', 'SDSS DR17', 'SDSS DR8'))


def mab(logf, lam):
    return -2.5 * (logf + 2 * np.log10(lam) - np.log10(C_AA)) - 48.60


rows = []
for k, (wnum, wden) in {'R1': ('3025', '5115'), 'R2': ('3025', '1282')}.items():
    ok = (t['COV' + wnum] >= COVER) & (t['COV' + wden] >= COVER) & (t['W' + wnum] > 0) & (t['W' + wden] > 0)
    s = t[ok].copy()
    s['logR'] = np.log10(s['W' + wnum] / s['W' + wden])
    # the window the quasar was normalized in: LOGFABS_NORM is there; move to 3025 with the measured ratio
    norm3025 = (s['W3025'] - 1).abs() < 1e-3
    s['logF3025'] = np.where(norm3025, s.LOGFABS_NORM, s.LOGFABS_NORM + s.logR)
    s['M3025'] = mab(s.logF3025, 3025.0)
    s['ratio'] = k
    rows.append(s[['ratio', 'SV', 'Z', 'logR', 'M3025', 'ALPHA_LAM', 'SBIN', 'LBIN']])
S = pd.concat(rows, ignore_index=True)


def choose_dm(x):
    lo, hi = np.percentile(x, [1, 99])
    for dm in (0.5, 1.0):
        e = np.arange(np.floor(lo / dm) * dm, hi + dm, dm); n, _ = np.histogram(x, e)
        if (n >= NMIN).all():
            return dm, e
    return 1.0, np.arange(np.floor(lo), hi + 1.0, 1.0)


stats = []
fig, ax = plt.subplots(2, 2, figsize=(15, 11))
zsplit = {'R1': [(0.2, 0.5), (0.5, 0.8), (0.8, 1.05)], 'R2': [(1.8, 2.0), (2.0, 2.25), (2.25, 2.45)]}
cols = {'DESI': 'C0', 'SDSS DR17': 'C3', 'SDSS DR8': 'C2'}; ls = ['-', '--', ':']
for i, k in enumerate(('R1', 'R2')):
    s = S[S.ratio == k]; dm, e = choose_dm(s.M3025)
    print(f'{k}: {len(s):,} quasars, M_AB(3025) bins of {dm} mag ({e[0]:.1f} .. {e[-1]:.1f})')
    for sv in cols:
        for j, (z0, z1) in enumerate(zsplit[k]):
            g = s[(s.SV == sv) & (s.Z >= z0) & (s.Z < z1)]
            if len(g) < NMIN:
                continue
            b = np.digitize(g.M3025, e) - 1
            for bi in range(len(e) - 1):
                y = g.logR[b == bi]
                if len(y) >= 30:
                    q16, q50, q84 = np.percentile(y, [16, 50, 84])
                    stats.append(dict(ratio=k, survey=sv, z_lo=z0, z_hi=z1, M_lo=e[bi], M_hi=e[bi + 1], n=len(y),
                                      logR_med=q50, logR_16=q16, logR_84=q84, dmag_med=-2.5 * q50))
            st = pd.DataFrame([x for x in stats if x['ratio'] == k and x['survey'] == sv and x['z_lo'] == z0])
            if len(st):
                xm = 0.5 * (st.M_lo + st.M_hi)
                ax[0, i].plot(xm, st.logR_med, ls[j], color=cols[sv], marker='o', ms=3, label=f'{sv} z {z0}-{z1}')
                if j == 1 or (k == 'R1' and j == 0):
                    ax[0, i].fill_between(xm, st.logR_16, st.logR_84, color=cols[sv], alpha=0.08)
    ax[0, i].set_xlabel('M_AB at 3025 A (absolute)'); ax[0, i].set_ylabel(f'log {k}  (median, 16-84 % band)')
    ax[0, i].invert_xaxis(); ax[0, i].legend(fontsize=8)
    sec = ax[0, i].secondary_yaxis('right', functions=(lambda y: -2.5 * y, lambda m: -m / 2.5)); sec.set_ylabel('m(3025) - m(other), f_lam mag')
    ax[0, i].set_title({'R1': 'R1 = F(3000-3050) / F(5100-5130)', 'R2': 'R2 = F(3000-3050) / F(1275-1290)'}[k] + f'   ({dm} mag bins)')
    # number per bin
    for sv in cols:
        g = s[s.SV == sv]
        if len(g):
            ax[1, i].hist(g.M3025, e, histtype='step', color=cols[sv], label=f'{sv} {len(g):,}')
    ax[1, i].set_yscale('log'); ax[1, i].invert_xaxis(); ax[1, i].set_xlabel('M_AB at 3025 A'); ax[1, i].set_ylabel('quasars per bin'); ax[1, i].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f'{OUT}/qso_ratio_vs_mag.png', dpi=90)
pd.DataFrame(stats).to_csv(f'{OUT}/qso_ratio_vs_mag.csv', index=False)
print(pd.DataFrame(stats).groupby(['ratio', 'survey', 'z_lo']).apply(lambda g: pd.Series(dict(n=g.n.sum(), logR=np.average(g.logR_med, weights=g.n)))).round(3).to_string())

# pooled lookup (all surveys, all overlap z): the ratio to chain each stack with, by M_AB(3025) in the chosen bins
look = []
for k in ('R1', 'R2'):
    s = S[S.ratio == k]; dm, e = choose_dm(s.M3025); b = np.digitize(s.M3025, e) - 1
    for bi in range(len(e) - 1):
        y = s.logR[b == bi]
        if len(y) >= NMIN:
            q16, q50, q84 = np.percentile(y, [16, 50, 84])
            look.append(dict(ratio=k, M_lo=e[bi], M_hi=e[bi + 1], n=len(y), logR_med=round(q50, 4), logR_16=round(q16, 4), logR_84=round(q84, 4),
                             R_med=round(10**q50, 4), sem=round(1.2533 * y.std() / np.sqrt(len(y)), 4)))
L = pd.DataFrame(look); L.to_csv(f'{OUT}/qso_ratio_lookup.csv', index=False)
print(L.to_string(index=False))
