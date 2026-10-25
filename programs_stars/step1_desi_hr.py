#!/usr/bin/env python3
"""Step 0b of the stellar thread: Gaia DR3 x DESI DR2 crossmatch (projects_gaia, 5 CSV chunks, 12.6 M rows) ->
quality sample (parallax S/N > 10, RUWE < 1.4, not non-single, not VARIABLE; ms.tex Sec. 'Stars'), HR coordinates
M_G0 = G + 5 log10(plx/100 mas) - A_G(GSP-Phot), (BP-RP)_0 = BP-RP - E(BP-RP)(GSP-Phot), and the number of HR cells
holding >= floor stars for several grids -> ~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet (3.92 M rows, 30 s).
Parallax zero point: global -0.017 mas (the full Lindegren+21 recipe needs nu_eff / ecl_lat; refine in pass S).
Usage: python3 step1_desi_hr.py
"""
import os, time
import numpy as np, pandas as pd

CSV = os.path.expanduser('~/github/projects_gaia/csvfiles/gaiadr3_desidr2pix_star{i}.csv')
OUT = os.path.expanduser('~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet')
COLS = ['targetid', 'parallax', 'parallax_over_error', 'ruwe', 'non_single_star', 'phot_variable_flag', 'phot_g_mean_mag', 'bp_rp',
        'has_xp_continuous', 'teff_gspphot', 'logg_gspphot', 'mh_gspphot', 'ag_gspphot', 'ebpminrp_gspphot']
FEH_BINS = [-9, -1.0, -0.5, 0.0, 9]       # 4 [M/H] bins; stars without GSP-Phot [M/H] form a 5th ("unknown") bin


def cells(df, dmg, dcol, feh_bins, floor):
    key = pd.DataFrame({'a': np.floor(df.MG0 / dmg), 'b': np.floor(df.BPRP0 / dcol)})
    if feh_bins is not None:
        key['f'] = pd.cut(df.mh_gspphot, feh_bins, labels=False).fillna(-1)
    n = key.value_counts()
    return int((n >= floor).sum()), int(n[n >= floor].sum()), len(n)


def main():
    t0 = time.time()
    g = pd.concat([pd.read_csv(CSV.format(i=i), usecols=COLS) for i in range(1, 6)], ignore_index=True)
    print('rows', f'{len(g):,}', 'unique targetid', f'{g.targetid.nunique():,}', f'{time.time() - t0:.0f}s')
    plx = g.parallax + 0.017
    q = (g.parallax_over_error > 10) & (g.ruwe < 1.4) & (g.non_single_star == 0) & (g.phot_variable_flag != 'VARIABLE') & (plx > 0)
    print('quality:', f'{q.sum():,}', ' plx S/N>50:', f'{(q & (g.parallax_over_error > 50)).sum():,}')
    print('XP (G<17.65):', f'{(g.has_xp_continuous == True).sum():,}', ' and quality:', f'{(q & (g.has_xp_continuous == True)).sum():,}')
    print('GSP-Phot [M/H]:', f'{g.mh_gspphot.notna().sum():,}', ' and quality:', f'{(q & g.mh_gspphot.notna()).sum():,}')
    s = g[q].copy()
    s['MG'] = s.phot_g_mean_mag + 5 * np.log10(plx[q] / 100.0)
    s['MG0'] = s.MG - s.ag_gspphot.fillna(0); s['BPRP0'] = s.bp_rp - s.ebpminrp_gspphot.fillna(0)
    print('\nHR cells with >= floor stars (grid M_G0 x (BP-RP)_0 [x [M/H]]):')
    for dmg, dcol in ((1.0, 0.2), (0.5, 0.1), (0.25, 0.05)):
        for feh in (None, FEH_BINS):
            for floor in (100, 300, 1000):
                nc, ns, nn = cells(s, dmg, dcol, feh, floor)
                print(f'  dMG={dmg} dcol={dcol} [M/H]={"no" if feh is None else "4 bins"} floor={floor:4d}: cells {nc:5d}  stars inside {ns:9,}  (populated cells {nn})')
    for lab, m in (('plx S/N>50', s.parallax_over_error > 50), ('XP available', s.has_xp_continuous == True)):
        nc, ns, _ = cells(s[m], 0.5, 0.1, None, 300)
        print(f'  {lab}: dMG=0.5 dcol=0.1 floor=300: cells {nc}, stars {ns:,}')
    print('M dwarfs (M_G0>10, BP-RP>2):', f'{((s.MG0 > 10) & (s.BPRP0 > 2)).sum():,}', ' WD region (M_G0>9, BP-RP<0.8):',
          f'{((s.MG0 > 9) & (s.BPRP0 < 0.8)).sum():,}', ' giants (M_G0<2, BP-RP>0.9):', f'{((s.MG0 < 2) & (s.BPRP0 > 0.9)).sum():,}')
    s[['targetid', 'MG', 'MG0', 'BPRP0', 'parallax', 'parallax_over_error', 'has_xp_continuous', 'mh_gspphot', 'teff_gspphot',
       'logg_gspphot', 'ag_gspphot']].to_parquet(OUT)
    print('wrote', OUT, f'{len(s):,}', f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
