#!/usr/bin/env python3
"""Eigenvalue histogram: the variance fraction of each component, lambda_k / sum lambda (Nao 2026-11-04).
The bases store the cumulative weighted R^2 after e0 (the mean) and after each added vector (VARIANCE extension);
the increments are the eigenvalues as fractions of the total weighted variance of the stacks about zero, and the
bars show them (a) as fractions of the variance left after the mean e0 and (b) of the variance the 10 vectors explain.
Usage: python3 20_eigval_hist.py [--basis name=path ...]; writes plots_luminosity/plots_qso_eigval_hist_<tag>.png
and csv_files/plots_qso_eigval_hist_<tag>.csv."""
import os, argparse, numpy as np, pandas as pd, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, '..'); D = os.path.join(R, 'data_eigenspectra')
DEF = [f'v3 (single basis)={D}/quasar_eigenspectra_v3.fits', f'v4 non-BAL={D}/quasar_eigenspectra_v4_nonbal.fits', f'v4 BAL={D}/quasar_eigenspectra_v4_bal.fits']
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--basis', action='append', default=None, help='label=path'); ap.add_argument('--tag', default='v4')
    a = ap.parse_args(); bases = [s.split('=', 1) for s in (a.basis or DEF)]
    fig, ax = plt.subplots(1, len(bases), figsize=(5 * len(bases), 4.2), sharey=True); ax = np.atleast_1d(ax); rows = []
    for (lab, path), x in zip(bases, ax):
        r2 = fitsio.read(path, 'VARIANCE'); lam = np.diff(r2)                      # lambda_k / total variance, k = 1..n
        rest = 1.0 - r2[0]                                                           # variance about the mean e0
        f_rest, f_sum = lam / rest, lam / lam.sum(); k = np.arange(1, len(lam) + 1)
        x.bar(k - 0.2, 100 * f_rest, 0.4, label='of the variance about e0'); x.bar(k + 0.2, 100 * f_sum, 0.4, label='of the variance in e1-e%d' % len(lam))
        for i, v in zip(k, f_sum): x.text(i + 0.2, 100 * v * 1.15, f'{100 * v:.1f}', ha='center', fontsize=7)
        x.set_yscale('log'); x.set_xticks(k); x.set_xticklabels([f'e{i}' for i in k]); x.set_xlabel('component')
        x.set_title(f'{lab}\ne0 (mean) explains {100 * r2[0]:.2f} %, e1-e{len(lam)} {100 * (r2[-1] - r2[0]):.2f} % of the total', fontsize=9); x.grid(alpha=0.3, axis='y')
        for i in k: rows.append(dict(basis=lab, k=i, lambda_over_total=lam[i - 1], frac_about_e0=f_rest[i - 1], frac_of_vectors=f_sum[i - 1], R2_cum=r2[i]))
    ax[0].set_ylabel('variance fraction  lambda_k / sum  [%]'); ax[0].legend(fontsize=8)
    fig.suptitle('Eigenvalue histogram: variance fraction per component'); fig.tight_layout()
    out = f'{R}/plots_luminosity/plots_qso_eigval_hist_{a.tag}.png'; fig.savefig(out, dpi=110); print(out)
    pd.DataFrame(rows).to_csv(f'{R}/csv_files/plots_qso_eigval_hist_{a.tag}.csv', index=False)
