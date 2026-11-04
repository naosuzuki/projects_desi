#!/usr/bin/env python3
"""Numbers quoted in ms.tex for the stellar bases v6 (hr2d) and v7 (hr2d_single): v7 against v6 per group, Redrock subtypes per 2D group and per v5 class, and the line strengths of the WD stars.
Overlap k: the smallest cosine between the spans of the first k vectors of the two bases (1 = one contains the other). Noise: 1.4826 x MAD of the first differences of a vector, median over e1-e10. Shape: rms of the ratio of a stack's G-normalised flux in v7 and v6, for stacks with the same name.
WD lines (S/N > 10, wd_line_ew.parquet): He-I-like is EW(5876) > 1.5 A at > 4 sigma, EW(4471) > 1 A and |EW(Hb)| < 3 A; DA-like is EW(Hb) > 5 A at > 4 sigma; the rest weak.
Output: csv_files/plots_star_v6v7_compare.csv, plots_star_hr_redrock.csv, plots_star_wd_lines.csv
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import hr_groups as H
D = '/Volumes/exdisk3/sdssdesi_stars'; REL = os.path.join(ROOT, 'data_eigenspectra'); CSV = os.path.join(ROOT, 'csv_files'); ZALL = os.path.expanduser('~/data/sdssdesi_cache/zall_stars_main.parquet'); G = ['WDHOT', 'WDA', 'WDB', 'WDCOOL', 'MSF', 'MSG', 'MSK', 'MSM0', 'MSM5', 'BINARY', 'SD', 'SGB']


def basis(v, g):
    p = os.path.join(REL, f'star_eigenspectra_{v}_desi_{g}.fits'); N = fitsio.read(p, ext='NORM'); k = np.char.strip(N['KIND'].astype(str)) == 'stack'
    return dict(M=fitsio.read(p, ext='MEAN').astype(float), E=fitsio.read(p, ext='EIGVEC').astype(float), V=fitsio.read(p, ext='VARIANCE').astype(float), P=fitsio.read(p, ext='PIXMASK') > 0, stars=int(N['N'][k].sum()), stacks=int(k.sum()))


def noise(E):
    d = np.diff(E, axis=1); return np.median(1.4826 * np.median(np.abs(d - np.median(d, axis=1, keepdims=True)), axis=1))


def stacks(b):
    F = fitsio.FITS(f'{D}/{b}/stacks_desi_hr2d.fits'); fl = F['FLUX'].read().astype(float); n = np.char.strip(F['CELLS'].read()['STACK'].astype(str)); F.close(); return dict(zip(n, fl))


def compare():
    a, b = stacks('hr2d'), stacks('hr2d_single'); rows = []
    for g in G:
        x, y = basis('v6', g), basis('v7', g); o = x['P'] & y['P']; r = dict(group=g, stars_v6=x['stars'], stars_v7=y['stars'], stacks_v6=x['stacks'], stacks_v7=y['stacks'], r2_10_v6=x['V'][-1], r2_10_v7=y['V'][-1]); ee = x['E'][:, o], y['E'][:, o]
        r['cos_e1'] = abs(ee[0][0] @ ee[1][0]) / np.linalg.norm(ee[0][0]) / np.linalg.norm(ee[1][0])
        for k in (3, 10):
            qa = np.linalg.qr(ee[0][:k].T)[0]; qb = np.linalg.qr(ee[1][:k].T)[0]; r[f'overlap_{k}'] = np.linalg.svd(qa.T @ qb, compute_uv=False).min()
        r['noise_ratio'] = noise(ee[1]) / noise(ee[0]); r['mean_pct_rms'] = 100 * np.sqrt(np.mean((y['M'][o] / x['M'][o] - 1)**2)); s = []
        for n in a:
            if n.endswith('_' + g) and n in b:
                u, v = a[n], b[n]; w = np.isfinite(u) & np.isfinite(v) & (u > 0) & (v > 0)
                if w.sum() > 3000: q = v[w] / u[w]; s.append(100 * np.std(q / np.median(q)))
        r['stack_shape_pct_rms'] = np.median(s) if s else np.nan; rows.append(r)
    t = pd.DataFrame(rows); t.to_csv(os.path.join(CSV, 'plots_star_v6v7_compare.csv'), index=False); print(t.round(3).to_string(index=False))


def redrock():
    h = H.load(); h['GROUP'] = H.assign(h, H.ms_ridge(h))[0]; z = pd.read_parquet(ZALL, columns=['TARGETID', 'SUBTYPE', 'ZCAT_PRIMARY', 'ZWARN']); z = z[z.ZCAT_PRIMARY & (z.ZWARN == 0)].drop_duplicates('TARGETID'); t = h.merge(z, on='TARGETID'); s = t.SUBTYPE.astype(str).str.strip().str[0]; s = s.where(s.isin(list('WABFGKM')), 'other').replace({'W': 'WD'}); x = pd.crosstab(t.GROUP, s).reindex(H.GROUPS); p = (100 * x.div(x.sum(1), axis=0)).round(1); p['N'] = x.sum(1); p.to_csv(os.path.join(CSV, 'plots_star_hr_redrock.csv')); print(p.to_string())


def wdlines():
    t = pd.read_parquet(f'{D}/hr2d/wd_line_ew.parquet'); s = t[t.SNR > 10]; he = (s.He5876 > 1.5) & (s.He5876 / s.He5876_s > 4) & (s.He4471 > 1) & (s.Hb.abs() < 3); da = (s.Hb > 5) & (s.Hb / s.Hb_s > 4) & ~he; s = s.assign(TYPE=np.where(he, 'He-I-like', np.where(da, 'DA-like', 'weak'))); x = pd.crosstab(s.GROUP, s.TYPE); x.to_csv(os.path.join(CSV, 'plots_star_wd_lines.csv')); print(x.to_string()); y = s[s.TYPE == 'He-I-like']; print('He-I-like colour 5, 50, 95 %', y.BPRP0.quantile([.05, .5, .95]).round(2).tolist(), 'min max', round(y.BPRP0.min(), 2), round(y.BPRP0.max(), 2)); hot = s[s.BPRP0 < 0.1]; print('colour < 0.1: median DW', hot.groupby('TYPE').DW.median().round(3).to_dict(), '16-84', hot.groupby('TYPE').DW.quantile([.16, .84]).round(2).to_dict())


if __name__ == '__main__':
    pd.set_option('display.width', 200); pd.set_option('display.max_columns', 30); compare(); redrock(); wdlines()
