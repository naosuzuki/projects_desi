#!/usr/bin/env python3
"""Observed-frame check of the stacks for atmospheric / instrumental features (Nao 2026-11-01).
Each stack's small-scale structure, f / running median (WIN pixels) - 1, is moved to the observed frame at the stack's
median redshift (lambda_obs = lambda_rest (1 + Z_MED)) and averaged over all stacks (median per observed pixel). Galaxy lines are
fixed in the rest frame and smear out over the redshifts; telluric bands, sky lines, camera seams and calibration features
are fixed in the observed frame and stay sharp. Also the same in the rest frame for comparison.
Output: PNG with both panels per survey, and a CSV of the observed-frame median (survey, wave_obs, med, n).
Usage: python3 25_obs_residuals.py --chain DIR [DIR ...] --out PNG
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import fitsio
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WIN = 51
OBS = np.arange(3550.0, 10000.0, 1.0)
MARK = {'O I 5577': 5577.3, 'O I 6300': 6300.3, 'Na D 5890': 5892.0, 'B band': 6870.0, 'A band': 7605.0,
        'H2O 7200': 7200.0, 'H2O 8200': 8200.0, 'H2O 9350': 9350.0, 'DESI b/r 5800': 5800.0, 'DESI r/z 7600': 7600.0,
        'SDSS b/r 6000': 6000.0}


def residuals(chains, max_stacks=None):
    rows_o, rows_r, w = [], [], None
    for d in chains:
        for fn in sorted(glob.glob(os.path.join(d, 'stacks_z*.fits'))):
            F = fitsio.FITS(fn); f = F['FLUX'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); c = F['CELLS'].read(); F.close()
            for i in range(len(c)):
                N = max(int(c['N_STACKED'][i]), 1); z = float(c['Z_MED'][i])
                good = np.isfinite(f[i]) & (ns[i] >= 0.5 * N) & (f[i] > 0)
                if good.sum() < 500:
                    continue
                y = np.where(good, f[i], np.nan); j = np.where(good)[0]
                sm = np.full(len(w), np.nan); sm[j] = median_filter(y[j], WIN, mode='nearest')
                r = np.where(good, y / sm - 1, np.nan)
                rows_r.append(r)
                rows_o.append(np.interp(OBS, w * (1 + z), r, left=np.nan, right=np.nan))
    return np.array(rows_o), np.array(rows_r), w


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--chain', nargs='+', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--label', nargs='+', default=None)
    a = ap.parse_args(); labels = a.label or [os.path.basename(d.rstrip('/')) for d in a.chain]
    fig, ax = plt.subplots(2 * len(a.chain), 1, figsize=(16, 3.4 * 2 * len(a.chain)))
    out = []
    for k, (d, lab) in enumerate(zip(a.chain, labels)):
        ro, rr, w = residuals([d])
        with np.errstate(all='ignore'):
            mo = np.nanmedian(ro, axis=0); no = np.isfinite(ro).sum(0); mr = np.nanmedian(rr, axis=0)
        mo[no < 20] = np.nan
        out.append(pd.DataFrame(dict(SURVEY=lab, WAVE_OBS=OBS, MED=mo, N=no)))
        a0, a1 = ax[2 * k], ax[2 * k + 1]
        a0.plot(OBS, mo, lw=0.5, color='k'); a0.set_xlim(OBS[0], OBS[-1]); a0.set_ylim(-0.03, 0.03); a0.grid(alpha=0.3)
        a0.set_title(f'{lab}: {len(ro)} stacks, median of f / running median - 1, OBSERVED frame (lambda_rest (1 + Z_MED))', fontsize=10)
        for name, x in MARK.items():
            a0.axvline(x, color='C3', lw=0.6, alpha=0.6); a0.text(x, 0.027, name, rotation=90, fontsize=7, va='top', color='C3')
        a1.plot(w, mr, lw=0.5, color='C0'); a1.set_xlim(1400, 9900); a1.set_ylim(-0.1, 0.1); a1.grid(alpha=0.3)
        a1.set_title(f'{lab}: same, REST frame', fontsize=10)
        top = pd.Series(np.abs(mo), index=OBS).dropna().nlargest(200)
        print(f'{lab}: {len(ro)} stacks; observed-frame rms {np.nanstd(mo):.4f}; strongest |median| at obs ' +
              ', '.join(f'{x:.0f}' for x in sorted(set((top.index // 20 * 20).astype(int)))[:40]))
    ax[-1].set_xlabel('wavelength [A]')
    fig.tight_layout(); fig.savefig(a.out, dpi=80); plt.close(fig)
    pd.concat(out).to_csv(os.path.splitext(a.out)[0] + '.csv', index=False); print(a.out)


if __name__ == '__main__':
    main()
