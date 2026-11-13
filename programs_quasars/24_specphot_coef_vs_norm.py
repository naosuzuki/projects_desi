"""Eigencoefficients against the normalization in a redshift slice (Nao 2026-11-12: 'plot c1 vs norm, c2 vs norm and do all
of the eigen coefficients'). Reads the non-BAL coefficient table of 21_specphot_fit.py, keeps the slice with a photometric
correction (PCORR_MODE 0-2) and a fit, and draws c_1..c_9 against M_AB(3025), the absolute AB magnitude of NORM (model flux in the 3000-3050 A pivot window,
F_lam(1+z) in erg/s/cm2/A, converted per quasar with its own distance modulus; top axis log NORM at the slice centre), as a 2D density of all
quasars with the running median and 16-84 % band, and the per-survey medians.
Outputs plots_luminosity/plots_qso_coef_vs_norm_z<zlo>-<zhi>_<tag>.png and csv_files/... (the running medians)."""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import qsolib

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70)
ap.add_argument('--tag', default='v1'); ap.add_argument('--nmin', type=int, default=200, help='quasars per bin for a median point'); a = ap.parse_args()
cols = ['SURVEY', 'Z', 'NORM', 'PCORR_MODE', 'CHI2', 'NPIX'] + [f'C{k}' for k in range(1, 10)]
d = pd.read_parquet(f'{PROD}/specphot/qso_specphot_coeff_{a.tag}.parquet', columns=cols)
s = d[(d.Z >= a.zlo) & (d.Z < a.zhi) & (d.NORM > 0) & d.PCORR_MODE.isin([0, 1, 2])].copy(); s['LOGN'] = np.log10(s.NORM)
zc = 0.5 * (a.zlo + a.zhi); dmc = qsolib.passB.COSMO.distmod(zc).value
K = -48.6 - 2.5 * np.log10(3025.0**2 / 2.998e18)                       # F_lam at 3025 A -> AB (F_nu) magnitude
dm = qsolib.passB.COSMO.distmod(s.Z.to_numpy()).value
s['MAB'] = -2.5 * (s.LOGN + 0.4 * dm) + K                                # absolute AB magnitude at 3025 A, per quasar
logn_of_M = lambda M: (K - M) / 2.5 - 0.4 * dmc; M_of_logn = lambda logn: -2.5 * (logn + 0.4 * dmc) + K
edges = np.arange(np.floor(s.MAB.quantile(0.001) * 4) / 4, s.MAB.quantile(0.999) + 0.25, 0.25); mid = 0.5 * (edges[:-1] + edges[1:])
ib = np.digitize(s.MAB, edges) - 1
surveys = [('desi', 'DESI DR2', 'tab:blue'), ('dr17', 'SDSS DR17', 'tab:orange'), ('dr8', 'SDSS DR8', 'tab:green')]
fig, axes = plt.subplots(3, 3, figsize=(14, 12), sharey=True); rows = []
for k, ax in zip(range(1, 10), axes.ravel()):
    c = s[f'C{k}'].to_numpy(); lo, hi = np.nanpercentile(c, [0.5, 99.5])
    ax.hist2d(c, s.MAB, bins=[np.linspace(lo, hi, 80), edges], cmap='Greys', norm=matplotlib.colors.LogNorm())
    for key, lab, col in [('all', 'all', 'k')] + surveys:
        sel = np.ones(len(s), bool) if key == 'all' else (s.SURVEY == key).to_numpy()
        med = np.full(len(mid), np.nan); p16 = med.copy(); p84 = med.copy(); n = np.zeros(len(mid), int)
        for j in range(len(mid)):
            m = sel & (ib == j); n[j] = m.sum()
            if n[j] >= a.nmin:
                med[j], p16[j], p84[j] = np.nanpercentile(c[m], [50, 16, 84])
        if key == 'all':
            ax.fill_betweenx(mid, p16, p84, color='tab:red', alpha=0.15); ax.plot(med, mid, color='tab:red', lw=2, label='all: median, 16-84 %')
        else:
            ax.plot(med, mid, color=col, lw=1.2, label=lab)
        rows += [dict(coef=f'C{k}', survey=key, mab_lo=edges[j], mab_hi=edges[j + 1], n=n[j], median=med[j], p16=p16[j], p84=p84[j]) for j in range(len(mid)) if n[j] > 0]
    ax.set_xlim(lo, hi); ax.set_xlabel(f'$c_{k}$'); ax.grid(alpha=0.3)
    if k == 1:
        ax.legend(fontsize=8, loc='upper left')
    if k % 3 == 0:
        right = ax.secondary_yaxis('right', functions=(logn_of_M, M_of_logn)); right.set_ylabel(f'log NORM at $z={zc:.3f}$ [erg s$^{{-1}}$ cm$^{{-2}}$ $\\AA^{{-1}}$]', fontsize=9); right.tick_params(labelsize=8)
for ax in axes[:, 0]:
    ax.set_ylabel('$M_{\\rm AB}(3025)$')
axes[0, 0].invert_yaxis()
fig.suptitle(f'Non-BAL eigencoefficients vs NORM as absolute AB magnitude at 3025 $\\AA$, $z={a.zlo}$-{a.zhi}: {len(s):,} quasars ' + ', '.join(f'{lab} {int((s.SURVEY == key).sum()):,}' for key, lab, _ in surveys), y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
stem = f'plots_qso_coef_vs_mag_z{a.zlo}-{a.zhi}_{a.tag}'
fig.savefig(f'{REPO}/plots_luminosity/{stem}.png', dpi=120); pd.DataFrame(rows).to_csv(f'{REPO}/csv_files/{stem}.csv', index=False); print(stem, len(s))
