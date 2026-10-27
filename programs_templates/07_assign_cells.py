#!/usr/bin/env python3
"""Step 5b: assign every galaxy to its cell of the 5-axis stacking grid (ms.tex Sec. 3.6, Table galgrid).

Axes: z (0.05, 0-1.6) x M (1 mag; band-switched absolute magnitude from the whole-galaxy flux: M_R z<0.35,
M_G 0.35-0.72, M_U 0.72-1.45, M_NUV 1.45-1.6) x color (0.25 mag; g-r for z<0.5, r-z above; dereddened LS
total fluxes) x line EW (log2 bins; Halpha for z<0.46, [OII] above) x Dn4000 (0.25 bins).
Input : ~/data/sdssdesi_products/passB_summary_columns.parquet (from 06_passB_summary.py)
Output: ~/data/sdssdesi_products/galaxy_cells.parquet   (one row per galaxy: TARGETID, PROGRAM, HEALPIX, cell indices, CELL_ID)
        ~/data/sdssdesi_products/galaxy_cell_table.csv  (one row per populated cell: indices, N, medians of the axes)
"""
import os, time
import numpy as np
import pandas as pd

P = os.path.expanduser('~/data/sdssdesi_products')
FLOOR = 300
LS_EXT = {'G': 3.214, 'R': 2.165, 'Z': 1.211}
EW_EDGES = np.array([1, 2, 4, 8, 16, 32, 64])                      # 8 bins: <1 ... >64 A
DN_EDGES = np.array([1.15, 1.40, 1.65, 1.90])                     # 5 bins
t0 = time.time()
d = pd.read_parquet(f'{P}/passB_summary_columns.parquet',
                    columns=['TARGETID', 'PROGRAM', 'HEALPIX', 'SPECTYPE', 'Z', 'ZWARN', 'ZCAT_PRIMARY', 'COADD_FIBERSTATUS', 'DELTACHI2',
                             'SNR_SDSSPIX', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'EBV', 'MABS_R', 'MABS_G', 'MABS_U', 'MABS_NUV',
                             'OII_SUM_EW', 'HALPHA_EW', 'DN4000', 'CAL_OK'])
g = d[(d.SPECTYPE == 'GALAXY') & (d.ZWARN == 0) & (d.COADD_FIBERSTATUS == 0) & d.ZCAT_PRIMARY.fillna(False).astype(bool)
      & (d.DELTACHI2 > 25) & (d.Z > 0.001) & (d.Z < 1.6)].copy()
del d
print(f'galaxies: {len(g):,}  ({time.time()-t0:.0f}s)', flush=True)

z = g.Z.to_numpy()
# axis 2: band-switched absolute magnitude (whole-galaxy flux, K from the spectrum)
M = np.where(z < 0.35, g.MABS_R, np.where(z < 0.72, g.MABS_G, np.where(z < 1.45, g.MABS_U, g.MABS_NUV)))
g['MBAND'] = np.where(z < 0.35, 'R', np.where(z < 0.72, 'G', np.where(z < 1.45, 'U', 'NUV')))
g['MABS'] = M
# axis 3: observed color from dereddened total fluxes
ebv = g.EBV.fillna(0).to_numpy()
with np.errstate(all='ignore'):
    mg = 22.5 - 2.5 * np.log10(g.FLUX_G) - LS_EXT['G'] * ebv
    mr = 22.5 - 2.5 * np.log10(g.FLUX_R) - LS_EXT['R'] * ebv
    mz = 22.5 - 2.5 * np.log10(g.FLUX_Z) - LS_EXT['Z'] * ebv
g['COLOR'] = np.where(z < 0.5, mg - mr, mr - mz)
g['COLOR'] = g.COLOR.replace([np.inf, -np.inf], np.nan)          # zero fluxes -> no color
g['MABS'] = g.MABS.replace([np.inf, -np.inf], np.nan)
# axis 4: line EW (Halpha below z=0.46, [OII] above; one scale, ratio ~0.95 in the overlap)
g['EW'] = np.where(z < 0.46, g.HALPHA_EW, g.OII_SUM_EW)
# bin indices
g['ZB'] = np.floor(z / 0.05).astype(int)
g['MB'] = np.floor(g.MABS).astype('Int64')
g['CB'] = np.floor(g.COLOR / 0.25).astype('Int64')
g['EB'] = pd.Series(np.searchsorted(EW_EDGES, np.nan_to_num(g.EW.to_numpy(), nan=-99)), index=g.index).astype('Int64')  # 0 = <1 A (incl. none/absorption)
g['DB'] = pd.Series(np.searchsorted(DN_EDGES, np.nan_to_num(g.DN4000.to_numpy(), nan=-9)), index=g.index).astype('Int64')
ok = g.MABS.notna() & g.COLOR.notna() & g.CAL_OK.fillna(False).astype(bool) & (g.MB >= -24) & (g.MB <= -16) & (g.CB >= 0) & (g.CB <= 7)
g['ASSIGNED'] = ok
print(f'assignable (finite M, color, calibrated, within the grid): {ok.sum():,} = {100*ok.mean():.1f}%', flush=True)
key = ['ZB', 'MB', 'CB', 'EB', 'DB']
a = g[ok]
cells = a.groupby(key, observed=True).agg(N=('TARGETID', 'size'), Z_MED=('Z', 'median'), MABS_MED=('MABS', 'median'), COLOR_MED=('COLOR', 'median'),
                                          EW_MED=('EW', 'median'), DN4000_MED=('DN4000', 'median'), SNR_MED=('SNR_SDSSPIX', 'median'),
                                          MBAND=('MBAND', 'first')).reset_index()
cells['STACK_SNR_EST'] = cells.SNR_MED * np.sqrt(cells.N)
cells = cells.sort_values(key).reset_index(drop=True)
cells['CELL_ID'] = np.arange(len(cells))
cells['POPULATED'] = cells.N >= FLOOR
g = g.merge(cells[key + ['CELL_ID', 'POPULATED']], on=key, how='left')
g[['TARGETID', 'PROGRAM', 'HEALPIX', 'Z', 'MABS', 'MBAND', 'COLOR', 'EW', 'DN4000', 'SNR_SDSSPIX'] + key + ['ASSIGNED', 'CELL_ID', 'POPULATED']].to_parquet(f'{P}/galaxy_cells.parquet')
cells.to_csv(f'{P}/galaxy_cell_table.csv', index=False)
print(f'cells with >=1: {len(cells):,}; >= {FLOOR}: {cells.POPULATED.sum():,} containing {100*cells[cells.POPULATED].N.sum()/len(a):.1f}% of assigned galaxies')
for fl in (100, 200, 300, 500, 1000):
    print(f'  floor {fl:>4}: {(cells.N >= fl).sum():>6,} cells, {100*cells[cells.N >= fl].N.sum()/len(a):5.1f}% of galaxies')
per = cells[cells.POPULATED].groupby('ZB').size()
print('populated cells per z slice:', ' '.join(f'{k*0.05:.2f}:{v}' for k, v in per.items()))
print(f'done {time.time()-t0:.0f}s')
