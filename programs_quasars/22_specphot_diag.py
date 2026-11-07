"""Diagnostics of the spectrophotometric correction (21_specphot_fit.py): per survey the distribution of the photometric
minus synthetic magnitude offsets dm per band (point sources / PSF objects with the band covered), the tilt a1 of the
corrected spectra, and the normalization NORM (pivot flux of the non-BAL fit) against the eigenspectra normalization of
the training quasars (LOGFABS_3025 of qso_cells_v3bal_all, brought back to the observed frame with the distance modulus).
Outputs plots_luminosity/plots_qso_specphot_diag_<tag>.png and csv_files/plots_qso_specphot_diag_<tag>.csv (histograms)."""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import qsolib

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BANDS = {'desi': 'GRZ', 'dr17': 'GRI', 'dr8': 'GRI'}; LABEL = {'desi': 'DESI DR2', 'dr17': 'SDSS DR17', 'dr8': 'SDSS DR8'}

ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='v1'); a = ap.parse_args()
cols = ['TARGETID', 'PLATE', 'MJD', 'FIBERID', 'SURVEY', 'Z', 'NORM', 'PCORR_MODE', 'PCORR_A0', 'PCORR_A1', 'BAL'] + [f'PCORR_DM_{b}' for b in 'GRIZ']
d = pd.concat([pd.read_parquet(f'{PROD}/specphot/qso_specphot_{k}_{a.tag}.parquet', columns=cols) for k in ('coeff', 'bal')], ignore_index=True)
t = pd.read_parquet(f'{PROD}/qso_cells_v3bal_all.parquet', columns=['TARGETID', 'PLATE', 'MJD', 'FIBERID', 'SURVEY', 'TRAIN', 'LOGFABS_3025']); t = t[t.TRAIN]

fig, ax = plt.subplots(3, 3, figsize=(13, 10)); rows = []
edges_dm = np.linspace(-1, 1, 81); edges_a1 = np.linspace(-3, 3, 61); edges_n = np.linspace(-0.5, 0.5, 81)
for i, sv in enumerate(['desi', 'dr17', 'dr8']):
    g = d[(d.SURVEY == sv) & d.PCORR_MODE.isin([0, 1, 2])]
    for b in BANDS[sv]:
        x = g[f'PCORR_DM_{b}'].dropna().to_numpy(); h, _ = np.histogram(x, edges_dm)
        ax[i, 0].step(edges_dm[:-1], h / h.sum(), where='post', label=f'{b.lower()}: med {np.median(x):+.3f}, 16/84 {np.percentile(x, 16):+.2f}/{np.percentile(x, 84):+.2f}')
        rows += [dict(survey=sv, panel='dm', band=b.lower(), lo=lo, hi=hi, n=n) for lo, hi, n in zip(edges_dm[:-1], edges_dm[1:], h)]
    ax[i, 0].set_title(f'{LABEL[sv]}: $\\Delta m$, {len(g):,} point sources', fontsize=10); ax[i, 0].legend(fontsize=8); ax[i, 0].set_xlabel('$\\Delta m$ [mag]')
    x = g[g.PCORR_MODE == 2].PCORR_A1.to_numpy(); h, _ = np.histogram(x, edges_a1); ax[i, 1].step(edges_a1[:-1], h / h.sum(), where='post')
    rows += [dict(survey=sv, panel='a1', band='', lo=lo, hi=hi, n=n) for lo, hi, n in zip(edges_a1[:-1], edges_a1[1:], h)]
    ax[i, 1].set_title(f'tilt $a_1$: {len(x):,} spectra ({len(x) / len(g):.0%}), med {np.median(x):+.2f}', fontsize=10); ax[i, 1].set_xlabel('$a_1$ [mag/dex]')
    key = ['TARGETID'] if sv == 'desi' else ['PLATE', 'MJD', 'FIBERID']
    m = t[t.SURVEY == sv].merge(d[d.SURVEY == sv], on=key, suffixes=('_t', ''))
    dm = qsolib.passB.COSMO.distmod(np.clip(m.Z.to_numpy(float), 1e-4, None)).value
    x = (np.log10(m.NORM) + 0.4 * dm - m.LOGFABS_3025).to_numpy(); xu = x + 0.4 * m.PCORR_A0.fillna(0).to_numpy()
    for v, lab in ((xu, 'before correction'), (x, 'after correction')):
        v = v[np.isfinite(v)]; h, _ = np.histogram(v, edges_n)
        ax[i, 2].step(edges_n[:-1], h / h.sum(), where='post', label=f'{lab}: med {np.median(v):+.3f}, 16/84 {np.percentile(v, 16):+.3f}/{np.percentile(v, 84):+.3f}')
        rows += [dict(survey=sv, panel='norm_' + lab.split()[0], band='', lo=lo, hi=hi, n=n) for lo, hi, n in zip(edges_n[:-1], edges_n[1:], h)]
    ax[i, 2].set_title(f'log NORM $-$ log norm(eigenspectra), {len(m):,} training', fontsize=10); ax[i, 2].legend(fontsize=8); ax[i, 2].set_xlabel('dex')
    for k in range(3):
        ax[i, k].set_ylabel('fraction')
fig.tight_layout(); fig.savefig(f'{REPO}/plots_luminosity/plots_qso_specphot_diag_{a.tag}.png', dpi=130)
pd.DataFrame(rows).to_csv(f'{REPO}/csv_files/plots_qso_specphot_diag_{a.tag}.csv', index=False); print('done')
