#!/usr/bin/env python3
"""Galaxy release v3: the two cell tables (version G: MABS_G axis; version 5450: M_AB(5450) axis; v3lib.py), from the per-galaxy
records of the v2 stacks (DESI <passC_v2_resp>/normtab, SDSS <passG_sdss_v2_2win>/normtab), before the re-read of the spectra.
The re-read (08_passC_stack.py --v3, 11_passG_sdss.py --v3) accumulates every cell with >= ACC_MIN_N members here; galaxies of
other cells are measured and tabled but not stacked.
Outputs: ~/data/sdssdesi_products/galaxy_cell_table_v3{G,5450}.csv (ZB, MB, CB, EB, DB, N, N_DESI, N_SDSS, medians),
         ~/data/sdssdesi_products/v3_galaxy_mags.png (distributions of both magnitudes), printed summary.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')
import v3lib, gridcells as gc

P = os.path.expanduser('~/data/sdssdesi_products')
COLS = ['Z', 'ZB', 'CB', 'EB', 'DB', 'EW', 'COLOR', 'DN4000', 'MABS_G', 'MABS_U', 'MABS_NUV', 'NORM_SCALE', 'RESCALE',
        'W5450', 'COV5450', 'W3050', 'COV3050']


def main():
    parts = [pd.read_parquet(f, columns=COLS).assign(SRC='DESI') for f in sorted(glob.glob(os.path.expanduser('~/data/sdssdesi_passC_v2_resp/normtab/norm_shard*.parquet')))]
    parts += [pd.read_parquet(f, columns=COLS).assign(SRC='SDSS') for f in sorted(glob.glob(os.path.expanduser('~/data/sdssdesi_passG_sdss_v2_2win/normtab/norm_dr*.parquet')))]
    t = pd.concat(parts, ignore_index=True); del parts
    print(f'{len(t):,} galaxies ({t.SRC.value_counts().to_dict()})', flush=True)
    dm = pb.COSMO.distmod(np.clip(t.Z.to_numpy(float), 1e-4, None)).value
    sc = (t.NORM_SCALE * t.RESCALE).to_numpy(float)
    f5 = np.where((t.COV5450 >= v3lib.MIN_GOOD) & (t.W5450 > 0), t.W5450 / sc, np.nan)
    f3 = np.where((t.COV3050 >= v3lib.MIN_GOOD) & (t.W3050 > 0), t.W3050 / sc, np.nan)
    # effective type needs a first MABS_G (measured or via the nominal type), then MABS_G is redone with the effective type
    mg0, _ = v3lib.mabs_g_full(t.MABS_G, t.MABS_U, t.MABS_NUV, t.DB)
    t['DB_EFF'] = v3lib.db_eff(t.Z, t.DB, t.EW, mg0)
    t['MABS_G_FULL'], t['MG_SRC'] = v3lib.mabs_g_full(t.MABS_G, t.MABS_U, t.MABS_NUV, t.DB_EFF)
    look = v3lib.load_lookup_5450()
    t['M_AB5450'], t['M5_SRC'] = v3lib.m5450(f5, f3, t.DB_EFF, dm, look)
    print('MABS_G source:', t.MG_SRC.value_counts().sort_index().to_dict(), '  M_AB(5450) source:', t.M5_SRC.value_counts().sort_index().to_dict())
    print('re-typed galaxies:', int((t.DB_EFF != t.DB).sum()))
    for c in ('MABS_G_FULL', 'M_AB5450'):
        q = np.nanpercentile(t[c], [0.5, 1, 50, 99, 99.5])
        print(f'{c}: 0.5/1/50/99/99.5 % ' + ' '.join(f'{v:.2f}' for v in q) + f'; by survey median ' +
              ', '.join(f'{s} {g[c].median():.2f}' for s, g in t.groupby('SRC')))
    d = (t.M_AB5450 - t.MABS_G_FULL)[t.M5_SRC == 0]
    print('M_AB(5450) - MABS_G (fiber light): ' + ', '.join(f'{s} {g.median():+.2f}' for s, g in d.groupby(t.SRC[t.M5_SRC == 0])))

    fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
    bins = np.arange(-26, -10, 0.25)
    for s, g in t.groupby('SRC'):
        ax[0].hist(g.MABS_G_FULL, bins, histtype='step', label=f'{s} MABS_G'); ax[0].hist(g.M_AB5450, bins, histtype='step', ls='--', label=f'{s} M_AB(5450)')
    ax[0].set_yscale('log'); ax[0].legend(fontsize=8); ax[0].set_xlabel('absolute magnitude'); ax[0].invert_xaxis()
    h = ax[1].hexbin(t.Z, t.MABS_G_FULL, gridsize=80, bins='log', cmap='viridis', mincnt=1, extent=(0, 1.6, -26, -12)); ax[1].set_ylabel('MABS_G (full)')
    ax[2].hexbin(t.Z, t.M_AB5450, gridsize=80, bins='log', cmap='viridis', mincnt=1, extent=(0, 1.6, -26, -12)); ax[2].set_ylabel('M_AB(5450)')
    for a_ in ax[1:]:
        a_.invert_yaxis(); a_.set_xlabel('z'); a_.axvline(0.8, color='w', lw=0.8)
    fig.suptitle('v3 luminosity axes: MABS_G (estimated from U / NUV above z ~ 0.75) and M_AB(5450) (from 3000-3100 A above z ~ 0.78)')
    fig.tight_layout(); fig.savefig(os.path.join(P, 'v3_galaxy_mags.png'), dpi=80); plt.close(fig)

    for ver, mcol, (lo, hi) in (('G', 'MABS_G_FULL', (gc.M_MIN, gc.M_MAX)), ('5450', 'M_AB5450', (v3lib.M5450_MIN, v3lib.M5450_MAX))):
        x = t.assign(MB=np.floor(t[mcol]).astype('Int64'))
        ok = x[mcol].notna() & x.COLOR.notna() & (x.MB >= lo) & (x.MB <= hi) & (x.CB >= gc.C_MIN) & (x.CB <= gc.C_MAX)
        x = x[ok]
        cells = x.groupby(gc.KEY).agg(N=('Z', 'size'), N_DESI=('SRC', lambda s_: int((s_ == 'DESI').sum())),
                                      Z_MED=('Z', 'median'), MABS_MED=(mcol, 'median'), COLOR_MED=('COLOR', 'median'),
                                      EW_MED=('EW', 'median'), DN4000_MED=('DN4000', 'median')).reset_index()
        cells['N_SDSS'] = cells.N - cells.N_DESI; cells['MAXIS'] = mcol
        cells = cells.sort_values(gc.KEY).reset_index(drop=True)
        cells.to_csv(os.path.join(P, f'galaxy_cell_table_v3{ver}.csv'), index=False)
        print(f'version {ver} ({mcol}, MB {lo}..{hi}): {ok.sum():,} galaxies in the grid, cells >= 100: {(cells.N >= 100).sum():,}, '
              f'>= 200: {(cells.N >= 200).sum():,} ({100 * cells[cells.N >= 200].N.sum() / ok.sum():.1f}% of galaxies)')


if __name__ == '__main__':
    main()
