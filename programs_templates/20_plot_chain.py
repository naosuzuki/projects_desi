#!/usr/bin/env python3
"""Plot the chained galaxy stacks (19_chain_stacks.py), like for like: one column per effective Dn4000 bin (DB_EFF, re-typed by
[OII] EW above z = 0.8); per z slice the galaxy-weighted mean of ALL its stacks of that type, coloured by z.
Top = each stack in its own window (5400-5500 below z = 0.30, 3000-3100 above), bottom = chained to 5400-5500.
Per pixel a stack enters only where NSPEC >= 0.5 N_STACKED and 100 A inside the observed range (3600-9824 A, DESI) at its
z; the mean is smoothed with a 9-pixel (~600 km/s) boxcar. Linear flux.
Usage: python3 20_plot_chain.py --chain DIR [DIR ...] --out PNG
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import fitsio
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DBLAB = {0: 'Dn4000 < 1.15', 1: '1.15-1.40', 2: '1.40-1.65', 3: '1.65-1.90', 4: '> 1.90'}
OBS = (3600 + 100, 9824 - 100); SMOOTH = 9; YMAX = {0: 3.0, 1: 2.0, 2: 1.6, 3: 1.5, 4: 1.5}


def smooth(y):
    ok = np.isfinite(y); num = uniform_filter1d(np.where(ok, y, 0.0), SMOOTH); den = uniform_filter1d(ok.astype(float), SMOOTH)
    return np.where(ok & (den > 0.5), num / np.where(den > 0, den, 1), np.nan)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--chain', nargs='+', required=True); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    fig, ax = plt.subplots(2, 5, figsize=(27, 11), sharex=True); cm = plt.get_cmap('turbo')
    ngal = {}
    for d in a.chain:
        ls = '-' if 'passC' in d else '--'
        for fn in sorted(glob.glob(os.path.join(d, 'stacks_z*.fits'))):
            F = fitsio.FITS(fn); c = pd.DataFrame(F['CELLS'].read()); f = F['FLUX'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); F.close()
            for db, g in c.groupby('DB_EFF'):
                idx = g.index.to_numpy(); N = g.N_STACKED.to_numpy(float); z = g.Z_MED.to_numpy(float)
                good = (ns[idx] >= 0.5 * N[:, None]) & (w[None, :] > OBS[0] / (1 + z[:, None])) & (w[None, :] < OBS[1] / (1 + z[:, None])) \
                       & np.isfinite(f[idx])
                wt = np.where(good, N[:, None], 0.0); sw = wt.sum(0)
                for r_, scale in ((1, np.ones(len(idx))), (0, 1.0 / g.CHAIN.to_numpy(float))):
                    y = np.where(sw > 0, (wt * np.where(good, f[idx] * scale[:, None], 0)).sum(0) / np.where(sw > 0, sw, 1), np.nan)
                    ax[r_, int(db)].plot(w, smooth(y), ls, color=cm(min(np.average(z, weights=N) / 1.6, 1)), lw=0.7)
                ngal[(d, int(db))] = ngal.get((d, int(db)), 0) + int(N.sum())
    for db in range(5):
        for r_ in range(2):
            ax[r_, db].set_ylim(0, YMAX[db]); ax[r_, db].set_xlim(1400, 9900); ax[r_, db].axhline(1, color='0.6', lw=0.5)
            for lo, hi, c_ in ((3000, 3100, 'C0'), (5400, 5500, 'C3')):
                ax[r_, db].axvspan(lo, hi, color=c_, alpha=0.15)
        nd = sum(v for (d, k), v in ngal.items() if k == db and 'passC' in d); nsd = sum(v for (d, k), v in ngal.items() if k == db and 'passC' not in d)
        ax[0, db].set_title(f'DB_EFF {db}  {DBLAB[db]}\nDESI {nd:,}, SDSS {nsd:,} galaxies', fontsize=11); ax[1, db].set_xlabel('rest wavelength [A]')
    ax[0, 0].set_ylabel('own window: F = 1 in 5400-5500 (z < 0.30) or 3000-3100 (z >= 0.30)')
    ax[1, 0].set_ylabel('chained: F = 1 in 5400-5500')
    sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, 1.6)); fig.colorbar(sm, ax=ax, label='z (galaxy-weighted)', fraction=0.015)
    fig.suptitle('Galaxy stacks before / after chaining: per (type, z slice) the galaxy-weighted mean of all stacks; type = Dn4000 bin, '
                 're-typed by [OII] EW above z = 0.8; solid DESI, dashed SDSS; shaded = the two normalization windows', fontsize=12)
    fig.savefig(a.out, dpi=75); plt.close(fig); print(a.out)


if __name__ == '__main__':
    main()
