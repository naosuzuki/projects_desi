#!/usr/bin/env python3
"""DESI stars whose fibermap has no Gaia BP, RP (stored as 0; Nao 2026-11-03). step5 and step6 take (BP-RP)_0 = 0 - E(BP-RP), so they land in the blue cells with their true M_G whatever their colour.
Per stack of the blue classes (WD, HOT): members with BP = RP = 0 in the pass tables, the Gaia catalog colour of those that have one, and the spectral colour 2.5 log10(F8100 / F4600), positive = red.
Reads tables only, no spectra. Output: csv_files/plots_star_missing_bprp.csv
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from step7_empca import star_class
TELL = '/Volumes/exdisk3/sdssdesi_stars/v2_tell'; TAG = 'desi_0.5x0.1_gspphot'; PASS = os.path.expanduser('~/data/sdssdesi_passS'); CACHE = os.path.expanduser('~/data/sdssdesi_cache'); WAVE = 10**(3.5502 + 1e-4 * np.arange(4447))


def main():
    F = fitsio.FITS(f'{TELL}/stacks_{TAG}.fits'); c = F['CELLS'].read(); fl = F['FLUX'].read().astype(float); F.close(); name = np.char.strip(c['STACK'].astype(str)); cls = np.array([star_class(m, k) for m, k in zip(c['MG0_MED'], c['COL0_MED'])])
    w = lambda a, b: np.nanmedian(np.where(fl > 0, fl, np.nan)[:, (WAVE > a) & (WAVE < b)], axis=1); sc = 2.5 * np.log10(w(8000, 8200) / w(4500, 4700)); blue = name[np.isin(cls, ['WD', 'HOT'])]
    m = pd.read_csv(f'{TELL}/members_{TAG}.csv', usecols=['SID', 'STACK']); m = m[m.STACK.isin(blue)].merge(pd.read_parquet(f'{CACHE}/stars/step6_selection_{TAG}.parquet', columns=['SID', 'TARGETID', 'FILE']), on='SID'); h = pd.read_parquet(f'{CACHE}/gaia_desi_stars_hr.parquet', columns=['targetid', 'BPRP0']).drop_duplicates('targetid').rename(columns=dict(targetid='TARGETID')); m = m.merge(h, on='TARGETID', how='left'); print(len(m), 'members of', len(blue), 'stacks;', m.FILE.nunique(), 'files', flush=True)
    rows = []
    def path(f):
        s = os.path.basename(f).replace('.fits', '').split('-'); return os.path.join(PASS, s[-2], f'hpx{s[-1]}.parquet') if len(s) >= 3 else None
    for f, g in m.groupby('FILE'):
        p = path(f)
        if p is None or not os.path.exists(p): continue
        x = pd.read_parquet(p, columns=['FILE', 'TARGETID', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG']); rows.append(x[x.TARGETID.isin(g.TARGETID)])
    x = pd.concat(rows).drop_duplicates('TARGETID'); m = m.merge(x[['TARGETID', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG']], on='TARGETID', how='left'); m['MISS'] = (m.GAIA_PHOT_BP_MEAN_MAG == 0) | (m.GAIA_PHOT_RP_MEAN_MAG == 0); print('pass rows found', m.GAIA_PHOT_BP_MEAN_MAG.notna().sum(), 'of', len(m), flush=True)
    k = dict(zip(name, zip(cls, c['MG0_MED'], c['COL0_MED'], sc))); t = m.groupby('STACK').agg(n=('SID', 'size'), miss=('MISS', 'sum'), cat=('BPRP0', 'count')).reset_index(); t['miss_frac'] = t.miss / t.n; v = m[m.MISS].groupby('STACK').BPRP0; t['miss_cat'] = t.STACK.map(v.count()).fillna(0).astype(int); t['miss_bprp_med'] = t.STACK.map(v.median())
    t['cls'] = t.STACK.map(lambda s: k[s][0]); t['mg'] = t.STACK.map(lambda s: k[s][1]); t['col'] = t.STACK.map(lambda s: k[s][2]); t['spec_col'] = t.STACK.map(lambda s: k[s][3]); t.sort_values(['cls', 'mg', 'STACK']).to_csv(os.path.join(ROOT, 'csv_files', 'plots_star_missing_bprp.csv'), index=False)
    pd.set_option('display.width', 180); pd.set_option('display.max_rows', 200); print(t[t.miss > 0].sort_values(['cls', 'mg'])[['STACK', 'cls', 'n', 'miss', 'miss_frac', 'miss_cat', 'miss_bprp_med', 'mg', 'col', 'spec_col']].round(3).to_string(index=False))
    print(t.groupby('cls')[['n', 'miss']].sum(), 'stacks with any', t[t.miss > 0].groupby('cls').size().to_dict(), 'over half', t[t.miss_frac > 0.5].groupby('cls').size().to_dict()); print('catalog colours of the missing', m[m.MISS].BPRP0.quantile([.05, .25, .5, .75, .95]).round(2).to_dict(), 'in catalog', m[m.MISS].BPRP0.notna().sum(), 'of', int(m.MISS.sum()), 'share bluer than 0.8 by catalog', round((m[m.MISS].BPRP0 < 0.8).sum() / max(m[m.MISS].BPRP0.notna().sum(), 1), 3)); print('with BP, RP:', m[~m.MISS].BPRP0.quantile([.05, .5, .95]).round(2).to_dict())


if __name__ == '__main__':
    main()
