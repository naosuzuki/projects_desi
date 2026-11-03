#!/usr/bin/env python3
"""ln(SDSS / DESI) of the same stellar cell in the telluric bands, v2 and corrected, from the residuals 06_star_compare.py saved. The star cancels, so what is left is the atmosphere of the two surveys and the noise.
Top row: median over cells with the band nuisance curve removed. Prints the rms per band, median over cells, and the same after a noise-free comparison: rms(v2)^2 - rms(corrected)^2.
Output: plots_luminosity/plots_telluric_starcompare.png, csv_files/plots_telluric_starcompare_crossband.csv
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import tellapply
Z = np.load('/Volumes/exdisk3/sdssdesi_stars/v2_tell/starcompare_cross.npz'); W, V, T = Z['wave'], Z['v2'], Z['tell']; BLUE, GRAY, INK, MUTED = '#2a78d6', '#a8a79f', '#0b0b0b', '#52514e'
NODES = [n for n in tellapply.basis('desi').bands if n[0] < 9200]


def rms(r):
    return np.sqrt(np.nanmean(r**2, axis=1))


def boot(x, n=500):
    x = x[np.isfinite(x)]; g = np.random.default_rng(1); return np.median(x), np.std([np.median(g.choice(x, len(x))) for _ in range(n)])


def main():
    good = Z['noise'] < 0.02; rows = []
    fig, ax = plt.subplots(2, len(NODES), figsize=(14, 5.2), sharey='row', gridspec_kw=dict(height_ratios=[3, 1.2]))
    for j, n in enumerate(NODES):
        i = (W >= n[0]) & (W <= n[-1]); a, b = V[good][:, i], T[good][:, i]; ra, rb = rms(a), rms(b); d = ra**2 - rb**2; x = W[i]; ma, mb = np.nanmedian(a, axis=0), np.nanmedian(b, axis=0)
        rows.append(dict(band=f'{n[0]:.0f}-{n[-1]:.0f}', cells=int(good.sum()), rms_v2=np.nanmedian(ra), rms_tell=np.nanmedian(rb), removed=np.sign(np.nanmedian(d)) * np.sqrt(abs(np.nanmedian(d))), removed_err=boot(d)[1] / (2 * np.sqrt(abs(np.nanmedian(d)))), frac_better=np.nanmean(d > 0), median_v2=np.nanstd(ma), median_tell=np.nanstd(mb)))
        q = ax[0, j]; q.axhline(0, color=GRAY, lw=0.5); q.plot(x, ma, color=GRAY, lw=1.4, label='v2'); q.plot(x, mb, color=BLUE, lw=1.0, label='telluric corrected'); q.set_title(f'{n[0]:.0f}–{n[-1]:.0f} Å', fontsize=9, loc='left', color=INK)
        p = ax[1, j]; p.axhline(0, color=GRAY, lw=0.5); p.plot(x, mb - ma, color=INK, lw=0.7); p.set_xlabel('wavelength [Å]', fontsize=8, color=MUTED)
        for u in (q, p):
            u.tick_params(labelsize=7, colors=MUTED, length=2); u.locator_params(axis='x', nbins=3); [u.spines[s].set_visible(False) for s in ('top', 'right')]; [u.spines[s].set_color(GRAY) for s in ('left', 'bottom')]
    ax[0, 0].set_ylabel('ln(SDSS / DESI), median of cells', fontsize=8, color=MUTED); ax[1, 0].set_ylabel('corrected − v2', fontsize=8, color=MUTED); h, l = ax[0, -1].get_legend_handles_labels(); fig.legend(h, l, loc='upper right', bbox_to_anchor=(0.99, 0.995), fontsize=9, frameon=False, ncol=2)
    fig.suptitle(f'Same stellar cells in SDSS and DESI, {good.sum()} cells with noise below 2 %', fontsize=11, x=0.045, ha='left', color=INK); fig.tight_layout()
    fig.savefig(os.path.join(ROOT, 'plots_luminosity', 'plots_telluric_starcompare.png'), dpi=120); t = pd.DataFrame(rows); t.to_csv(os.path.join(ROOT, 'csv_files', 'plots_telluric_starcompare_crossband.csv'), index=False); print(t.round(5).to_string())


if __name__ == '__main__':
    main()
