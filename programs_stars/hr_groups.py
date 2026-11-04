#!/usr/bin/env python3
"""Stellar groups chosen in 2D on the Gaia DR3 HR diagram, (BP-RP)_0 and M_G0 of gaia_desi_stars_hr.parquet (Nao 2026-11-03: a 1D colour cut on the cell median let contaminants in; the DESI fibermap BP, RP are not used).
Ridges: the WD sequence, M_G0 against colour, and the main sequence (mode per colour bin); a star's offset from the ridge is DM = M_G0 - ridge.
WD: M_G0 > 7 + 3.5 colour, M_G0 > 8, colour < 1.3, |DM| < 1. WD_HOT colour < 0.1, one branch. WD_A 0.1-0.75 and DM < 0.1, the bright branch, mostly DA. WD_B the same colours and DM >= 0.1, the faint branch, about half H-weak (DB, DC). WD_COOL >= 0.75.
Main sequence, colour >= 0.55: MS |DM| band -0.45 to 0.6 by colour, MS_F 0.55-0.75, MS_G 0.75-1.05, MS_K 1.05-1.9, MS_M0 1.9-2.7, MS_M5 >= 2.7. BINARY DM -1.1 to -0.45 at colour >= 0.85. SD DM > 0.6. SGB above those with M_G0 >= 3. GIANT M_G0 < 3. BLUE colour < 0.55, M_G0 < 8. OTHER the rest (WD outside |DM| < 1).
Usage: from hr_groups import assign, GROUPS. __main__ draws the groups: plots_luminosity/plots_star_hr_groups.png, csv_files/plots_star_hr_groups.csv
"""
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); CAT = os.path.expanduser('~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet')
WD_POLY = np.array([1.79, -4.558, 3.304, -1.279, 3.53, 11.806]); MS_COL = np.arange(0.5, 4.01, 0.1)
GROUPS = ['WD_HOT', 'WD_A', 'WD_B', 'WD_COOL', 'MS_F', 'MS_G', 'MS_K', 'MS_M0', 'MS_M5', 'BINARY', 'SD', 'SGB', 'GIANT', 'BLUE', 'OTHER']


def load():
    return pd.read_parquet(CAT, columns=['targetid', 'MG0', 'BPRP0', 'parallax_over_error']).drop_duplicates('targetid').rename(columns=dict(targetid='TARGETID')).dropna(subset=['MG0', 'BPRP0'])


def ms_ridge(h):
    """mode of M_G0 per colour bin, stars with parallax S/N > 20, on the lower sequence (M_G0 > 3, not the WD locus)."""
    g = h[(h.parallax_over_error > 20) & (h.MG0 > 3) & (h.MG0 <= 7 + 3.5 * h.BPRP0) ]; y = []
    for c in MS_COL:
        x = g.MG0[(g.BPRP0 >= c - 0.05) & (g.BPRP0 < c + 0.05)]
        if len(x) < 30: y.append(np.nan); continue
        n, e = np.histogram(x, bins=np.arange(3, 18, 0.1)); s = np.convolve(n, np.ones(5) / 5, 'same'); y.append(e[np.argmax(s)] + 0.05)
    y = np.array(y); k = np.isfinite(y); return np.polyfit(MS_COL[k], y[k], 7)


def assign(h, pms):
    c, m = h.BPRP0.to_numpy(float), h.MG0.to_numpy(float); out = np.full(len(h), 'OTHER', dtype='U8'); dw = m - np.polyval(WD_POLY, c); dm = m - np.polyval(pms, np.clip(c, 0.5, 4.0))
    wd = (m > 7 + 3.5 * c) & (m > 8) & (c < 1.3) & (np.abs(dw) < 1); notwd = ~((m > 7 + 3.5 * c) & (m > 8) & (c < 1.3))
    out[wd & (c < 0.1)] = 'WD_HOT'; out[wd & (c >= 0.1) & (c < 0.75) & (dw < 0.1)] = 'WD_A'; out[wd & (c >= 0.1) & (c < 0.75) & (dw >= 0.1)] = 'WD_B'; out[wd & (c >= 0.75)] = 'WD_COOL'
    u = notwd & (c >= 0.55); ms = u & (dm > -0.45) & (dm <= 0.6)
    for n, a, b in (('MS_F', 0.55, 0.75), ('MS_G', 0.75, 1.05), ('MS_K', 1.05, 1.9), ('MS_M0', 1.9, 2.7), ('MS_M5', 2.7, 9)): out[ms & (c >= a) & (c < b)] = n
    out[u & (dm > 0.6)] = 'SD'; up = u & (dm <= -0.45); out[up & (dm > -1.1) & (c >= 0.85) & (m >= 3)] = 'BINARY'; out[up & (out == 'OTHER') & (m >= 3)] = 'SGB'; out[u & (m < 3)] = 'GIANT'; out[notwd & (c < 0.55) & (m < 8)] = 'BLUE'
    return out, dw, dm


if __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    h = load(); p = ms_ridge(h); g, dw, dm = assign(h, p); h['GROUP'] = g; print('MS ridge', np.round(p, 5).tolist()); z = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_cache/zall_stars_main.parquet'), columns=['TARGETID', 'PROGRAM', 'SUBTYPE', 'ZCAT_PRIMARY', 'ZWARN']); z = z[z.ZCAT_PRIMARY & (z.ZWARN == 0)].drop_duplicates('TARGETID'); t = h.merge(z, on='TARGETID', how='left')
    pd.set_option('display.width', 180); pd.set_option('display.max_columns', 20); tab = pd.crosstab(t.GROUP, t.PROGRAM.fillna('none')).reindex(GROUPS).fillna(0).astype(int); print(tab.to_string()); s = pd.crosstab(t.GROUP, t.SUBTYPE.astype(str).str.strip().replace('nan', 'none')).reindex(GROUPS).fillna(0).astype(int); print((100 * s.div(s.sum(1), axis=0)).round(1).to_string()); tab.to_csv(os.path.join(ROOT, 'csv_files', 'plots_star_hr_groups.csv'))
    fig, ax = plt.subplots(figsize=(8.5, 10)); cols = dict(WD_HOT='#2a78d6', WD_A='#1f4fa8', WD_B='#d6402a', WD_COOL='#8a6bbf', MS_F='#e0a100', MS_G='#b8860b', MS_K='#2e8b57', MS_M0='#7a3b1d', MS_M5='#c2185b', BINARY='#555555', SD='#00a3a3', SGB='#9e9ac8', GIANT='#e6550d', BLUE='#3182bd', OTHER='#c8c8c8')
    for n in GROUPS[::-1]:
        x = h[h.GROUP == n]; x = x.sample(min(len(x), 40000), random_state=1); ax.scatter(x.BPRP0, x.MG0, s=.4, color=cols[n], rasterized=True, label=f'{n} {int((h.GROUP == n).sum()):,}', linewidths=0)
    ax.set_xlim(-1, 5); ax.set_ylim(17, -5); ax.set_xlabel('(BP-RP)0'); ax.set_ylabel('M_G0'); ax.legend(markerscale=12, fontsize=8, frameon=False, loc='upper right'); fig.tight_layout(); f = os.path.join(ROOT, 'plots_luminosity', 'plots_star_hr_groups.png'); fig.savefig(f, dpi=110); print(f)
