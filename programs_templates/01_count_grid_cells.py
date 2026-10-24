"""Count populated cells of the galaxy stacking grid (Sec. 3.6 of the SDSS-DESI paper).

Axes measurable from the catalog today: redshift, total luminosity M_r, observed color.
The line-strength (Halpha / [OII] EW) and Dn4000 axes come from the coadd pass and are
not counted here.

Grid:  z in 0.05 slices (0-1.6) x M_r in 1-mag bins (-24..-16) x color in 0.25-mag bins (0..2)
       color = g-r for z<0.5, r-z for z>=0.5 (observed, Legacy Surveys total fluxes)
Sample: DESI DR2 loa main dark+bright, ZCAT_PRIMARY, ZWARN=0, COADD_FIBERSTATUS=0,
        SPECTYPE=GALAXY, DELTACHI2>25, 0<z<1.6.
Input:  ~/data/sdssdesi_cache/zall_loa.parquet (built from zall-pix-loa.fits)
"""
import os
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as cosmo

CACHE = os.path.expanduser('~/data/sdssdesi_cache/zall_loa.parquet')
DZ, DM, DC = 0.05, 1.0, 0.25
ZMAX = 1.6

z = pd.read_parquet(CACHE, columns=['SURVEY', 'PROGRAM', 'Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2',
                                    'COADD_FIBERSTATUS', 'ZCAT_PRIMARY', 'FLUX_G', 'FLUX_R', 'FLUX_Z'])
g = z[(z.SURVEY == 'main') & z.PROGRAM.isin(['dark', 'bright']) & z.ZCAT_PRIMARY & (z.ZWARN == 0)
      & (z.COADD_FIBERSTATUS == 0) & (z.SPECTYPE == 'GALAXY') & (z.DELTACHI2 > 25)
      & (z.Z > 0.001) & (z.Z < ZMAX) & (z.FLUX_G > 0) & (z.FLUX_R > 0) & (z.FLUX_Z > 0)].copy()

zgrid = np.linspace(0, ZMAX, 161)
dm = np.interp(g.Z, zgrid, cosmo.distmod(np.clip(zgrid, 1e-3, None)).value)
g['Mr'] = 22.5 - 2.5 * np.log10(g.FLUX_R) - dm                      # total luminosity, no K-correction
g['col'] = np.where(g.Z < 0.5, -2.5 * np.log10(g.FLUX_G / g.FLUX_R), -2.5 * np.log10(g.FLUX_R / g.FLUX_Z))
g['zb'] = np.floor(g.Z / DZ).astype(int)
g['Mb'] = np.floor(g.Mr / DM).astype(int)
g['cb'] = np.floor(g.col / DC).astype(int)

print(f'galaxies used: {len(g):,}')
print('M_r 1-99%:', np.percentile(g.Mr, [1, 99]).round(1), ' color 1-99%:', np.percentile(g.col, [1, 99]).round(2))
nz, nm, nc = int(ZMAX / DZ), 8, 8
print(f'nominal 3-axis cells: {nz} x {nm} x {nc} = {nz*nm*nc}')
cells = g.groupby(['zb', 'Mb', 'cb']).size()
for n in [100, 300, 1000, 3000, 10000]:
    print(f'cells with >= {n:>5} galaxies: {(cells >= n).sum():>5}   galaxies in them: {cells[cells >= n].sum()/len(g)*100:5.1f}%')
per = cells[cells >= 300].groupby(level=0).size()
print('populated (>=300) per z slice:', ' '.join(f'{k*DZ:.2f}:{v}' for k, v in per.items()))
cells.rename('N').reset_index().to_csv(os.path.expanduser('~/github/projects_desi/data/grid_cells_z_Mr_color.csv'), index=False)
