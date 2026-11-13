"""Held-out rms of the brightness prediction against redshift (Nao 2026-11-12, rounds 6 and 7):
   the pooled cosmology-free model of round 6w (one model, z 0.3-2.2, target M_rel = apparent magnitude minus its 0.05-bin
   median, no redshift input; per-bin table written by 26_specphot_ml.py --save-pred) against the per-bin models of round 7
   (one model per 0.1 bin, csv_files/qso_ml_bench_zbins_v1_r7.csv).  Output plots_luminosity/plots_qso_ml_zbins_v1.png."""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='v1'); ap.add_argument('--target', default='M_rel'); a = ap.parse_args()
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
col = {'ridge': 'tab:gray', 'HGB tuned': 'tab:blue', 'torch MLP': 'tab:red'}
for ax, tgt in zip(axs, [a.target, a.target + '_raw']):
    first = True
    for f in sorted(glob.glob(f'{REPO}/csv_files/qso_ml_zbins_*_A2_{tgt}_r6w.csv')):
        m = os.path.basename(f).split('_A2_')[0].replace('qso_ml_zbins_', '').replace('_', ' '); d = pd.read_csv(f); d = d[d.n > 2000]
        if first:
            ax.plot(d.zlo + 0.025, d.sigma, 'k:', lw=1.5, label='$\\sigma$ of the target (no prediction)'); first = False
        ax.plot(d.zlo + 0.025, d.rms, '-', color=col.get(m, 'k'), lw=2, label=f'pooled model, {m}')
        ax.plot(d.zlo + 0.025, d.rms_bright, '--', color=col.get(m, 'k'), lw=1, label=f'pooled, {m}, bright half')
    p7 = f'{REPO}/csv_files/qso_ml_bench_zbins_{a.tag}_r7.csv'
    if os.path.exists(p7):
        b = pd.read_csv(p7); b = b[(b.target == tgt) & (b.error.isna() | (b.error == ''))]; b['zlo'] = b.slice.str.split('-').str[0].astype(float)
        for m, g in b.groupby('method'):
            g = g.sort_values('zlo'); ax.plot(g.zlo + 0.05, g.rms, 'o', color=col.get(m, 'k'), ms=5, mfc='none', label=f'per-bin models, {m}')
    ax.set_xlabel('redshift'); ax.set_title(f'target {tgt.replace("_", " ")}', fontsize=10); ax.grid(alpha=0.3)
axs[0].set_ylabel('held-out rms  [mag]'); axs[0].legend(fontsize=7.5, loc='upper left', ncol=2)
fig.suptitle('Brightness from the coefficients against redshift: one cosmology-free model (round 6w) against one model per 0.1 bin (round 7)', fontsize=10)
fig.tight_layout(); out = f'{REPO}/plots_luminosity/plots_qso_ml_zbins_{a.tag}.png'; fig.savefig(out, dpi=120); print(out)
