#!/usr/bin/env python3
"""Luminosity-binned galaxy spectra from the chained stacks (19_chain_stacks.py; Nao 2026-10-31: stack by absolute magnitude).
One column per type (DB_EFF: Dn4000 bin, re-typed by [OII] EW above z = 0.8). Per 1-mag bin of the stack's median MABS_G
(the luminosity at the reference, MABS_G_MED) the galaxy-weighted mean over ALL chained stacks of that type and luminosity,
every redshift together, so one spectrum per luminosity bin covers the UV (high-z stacks) to the optical (low-z stacks).
All on one scale (F = 1 in 5400-5500 A), then shown as ABSOLUTE AB magnitude: F_nu = F_lambda lambda^2, anchored in rest g
(4600-4950 A) to the galaxy-weighted mean MABS_G of the bin:  M_AB(lambda) = <MABS_G> - 2.5 log10(F_nu(lambda) / F_nu(g)).
Per pixel a stack enters only where NSPEC >= 0.5 N_STACKED and 100 A inside the observed range at its z, with its weight
N_STACKED tapered linearly to 0 over the last TAPER pixels at each end of that range (groups of stacks fade in and out instead
of switching); 9-pixel boxcar. Bottom row: galaxies contributing per pixel.
--by ew: columns are emission-line classes from the EW axis instead (passive < 2, weak 2-8, moderate 8-32, strong > 32 A).
Usage: python3 21_plot_lum.py --chain DIR [DIR ...] --out PNG [--dm 1.0] [--by db|ew]
The PNG goes where --out says (Nao 2026-11-02: plots_luminosity/ in the repo root); the table of the plotted curves goes to csv_files/ in the repo root.
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
GROUPS = {'db': ('DB_EFF', lambda v: int(v), {k: f'DB_EFF {k}  {v}' for k, v in DBLAB.items()}),
          # emission-line class from the EW axis (8 log2 bins: <1, 1-2, 2-4, ..., 32-64, >64 A; Halpha z < 0.46, [OII] above), pairs merged
          'ew': ('EB', lambda v: int(v) // 2, {0: 'passive  EW < 2 A', 1: 'weak  EW 2-8 A', 2: 'moderate  EW 8-32 A', 3: 'strong  EW > 32 A'})}
OBS = (3600 + 100, 9824 - 100); SMOOTH = 9; MRANGE = (-24.0, -16.0); TAPER = 200; GWIN = (4600, 4950)
FLIM = (5e-4, 30.0)                                        # F_lambda at 10 pc [erg/s/cm^2/A]


CSV_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'csv_files')   # Nao 2026-11-02

def smooth(y):
    ok = np.isfinite(y); num = uniform_filter1d(np.where(ok, y, 0.0), SMOOTH); den = uniform_filter1d(ok.astype(float), SMOOTH)
    return np.where(ok & (den > 0.5), num / np.where(den > 0, den, 1), np.nan)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--chain', nargs='+', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--dm', type=float, default=1.0); ap.add_argument('--by', default='db', choices=list(GROUPS))
    ap.add_argument('--axis', default='G', choices=['G', '5450'], help='luminosity axis: MABS_G (anchor rest g) or M_AB(5450) (anchor 5400-5500)')
    ap.add_argument('--title', default='')
    ap.add_argument('--unit', default='mag', choices=['mag', 'flam'], help='y axis: absolute AB mag, or F_lambda at 10 pc (log)')
    a = ap.parse_args(); gcol, gfun, glab = GROUPS[a.by]; ncol = len(glab)
    MCOL, GW, MLAB = ('MABS_G_MED', GWIN, 'M_g') if a.axis == 'G' else ('M5450_MED', (5400, 5500), 'M_5450')
    MR = MRANGE if a.axis == 'G' else (-24.0, -14.0)
    edges = np.arange(MR[0], MR[1] + 1e-9, a.dm); nb = len(edges) - 1
    acc = {}                                               # (db, mbin) -> [sum N f, sum N]
    w = None
    for d in a.chain:
        for fn in sorted(glob.glob(os.path.join(d, 'stacks_z*.fits'))):
            F = fitsio.FITS(fn); c = pd.DataFrame(F['CELLS'].read()); f = F['FLUX'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); F.close()
            mb = np.digitize(c[MCOL], edges) - 1
            for i in range(len(c)):
                if not (0 <= mb[i] < nb):
                    continue
                N = float(c.N_STACKED.iat[i]); z = float(c.Z_MED.iat[i])
                good = (ns[i] >= 0.5 * N) & (w > OBS[0] / (1 + z)) & (w < OBS[1] / (1 + z)) & np.isfinite(f[i])
                if not good.any():
                    continue
                j = np.arange(len(w)); g0, g1 = j[good].min(), j[good].max()
                wt = N * np.clip(np.minimum(j - g0, g1 - j) / TAPER, 0, 1) * good
                k = (gfun(c[gcol].iat[i]), int(mb[i]))
                if k not in acc:
                    acc[k] = [np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0.0, 0.0]
                acc[k][0] += np.where(good, wt * f[i], 0.0); acc[k][1] += wt; acc[k][2] += np.where(good, N, 0.0)
                acc[k][3] += N * float(c[MCOL].iat[i]); acc[k][4] += N
    fig, ax = plt.subplots(2, ncol, figsize=(5.4 * ncol, 10.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    cm = plt.get_cmap('rainbow_r'); norm = plt.Normalize(*MR)          # faint purple -> bright red
    rows = []
    gw = (w >= GW[0]) & (w <= GW[1])
    for (db, b), (sf, swt, sn, sm_, sN) in sorted(acc.items()):
        mc = 0.5 * (edges[b] + edges[b + 1]); mg_ = sm_ / sN
        y = smooth(np.where(swt > 0, sf / np.where(swt > 0, swt, 1), np.nan)); fnu = y * w**2
        ref = np.nanmean(fnu[gw])
        if not np.isfinite(ref) or ref <= 0:
            continue
        with np.errstate(all='ignore'):
            mab = mg_ - 2.5 * np.log10(fnu / ref)
        yv = mab if a.unit == 'mag' else 10**(-0.4 * (mab + 48.6)) * 2.998e18 / w**2      # F_lambda at 10 pc, erg/s/cm^2/A
        ax[0, db].plot(w, yv, color=cm(norm(mc)), lw=0.9, label=f'{MLAB} {edges[b]:.0f} to {edges[b+1]:.0f}  (<{MLAB}> {mg_:.2f}, max {int(sn.max()):,})')
        ax[1, db].plot(w, sn, color=cm(norm(mc)), lw=0.9)
        rows.append(dict(GROUP=db, M_lo=edges[b], M_hi=edges[b + 1], MABS_G_mean=round(mg_, 3), max_galaxies_per_pixel=int(sn.max()),
                         wave_lo=float(w[sn > 0].min()), wave_hi=float(w[sn > 0].max())))
    for db in range(ncol):
        a0 = ax[0, db]; a0.set_xlim(1400, 9900); a0.grid(alpha=0.3)
        if a.unit == 'mag':
            a0.set_ylim(MR[1] + 2, -25.5)
        else:
            a0.set_yscale('log'); a0.set_ylim(*FLIM)
        a0.axvspan(5400, 5500, color='C3', alpha=0.15); a0.axvspan(3000, 3100, color='C0', alpha=0.15)
        a0.set_title(glab[db], fontsize=11); a0.legend(fontsize=7, loc='lower right'); a0.axvspan(*GW, color='0.8', alpha=0.3)
        ax[1, db].set_yscale('log'); ax[1, db].set_xlabel('rest wavelength [A]')
    ax[0, 0].set_ylabel(f'absolute AB magnitude M_AB(lambda)  (F_nu, anchored to <{MLAB}>)' if a.unit == 'mag' else
                        f'F_lambda at 10 pc [erg s^-1 cm^-2 A^-1]  (anchored to <{MLAB}>)'); ax[1, 0].set_ylabel('galaxies per pixel')
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm); fig.colorbar(sm, ax=ax, label=f'{MLAB} (stack median)', fraction=0.015)
    fig.suptitle(a.title + f'Galaxy spectra by luminosity ({a.dm:.0f}-mag bins of {MLAB}) and ' + ('Dn4000 type' if a.by == 'db' else
                 'emission-line class (Halpha EW z < 0.46, [OII] above)') + '\nchained stacks over all z, '
                 f'galaxy-weighted, edges tapered; ' + ('absolute AB mag' if a.unit == 'mag' else 'F_lambda at 10 pc') + f' anchored to <{MLAB}> in the grey band; shaded = normalization windows', fontsize=12)
    fig.savefig(a.out, dpi=75); plt.close(fig)
    csv = os.path.join(CSV_DIR, os.path.splitext(os.path.basename(a.out))[0] + '.csv'); os.makedirs(CSV_DIR, exist_ok=True)
    pd.DataFrame(rows).to_csv(csv, index=False); print(a.out, csv)


if __name__ == '__main__':
    main()
