#!/usr/bin/env python3
"""Eigenvalue histogram (Nao 2026-11-04): the variance fraction of each component, lambda_k / sum lambda, for one or more galaxy
bases. Two estimates per vector: (a) the increment of the cumulative R^2 (VARIANCE extension: R^2 of the residuals about e0 for
the passive and centered bases) times the total weighted variance, i.e. the variance the vector explains in the weighted EM-PCA
sense; (b) lambda_k = sum_i w_i c_ik^2 from the orthonormal coefficients (COEFF_ORTHO), w_i the stack's mean pixel weight.
The bar chart shows (a); (b) is printed. Output plots_luminosity/gal_eigenvalues_<tags>.png, csv_files/gal_eigenvalues.csv.
Usage: python3 35_eigenvalues.py FITS [FITS ...]
"""
import os, sys
import numpy as np, pandas as pd, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def eig(fn):
    F = fitsio.FITS(fn); R2 = F['VARIANCE'].read(); hdr = F['EIGVEC'].read_header()
    frac_a = np.diff(np.concatenate([[0.0], R2]))            # increments of R^2 about e0: fraction of the residual variance
    out = dict(tag=os.path.basename(fn).replace('galaxy_empca_', '').replace('.fits', ''), R2=R2, frac=frac_a / R2[-1], unexpl=1 - R2[-1])
    try:
        C = F['COEFF_ORTHO'].read().astype(float); k0 = 1 if hdr.get('CENTER', 0) or hdr.get('E0', '') in ('passive', 'select') else 0
        lam = (C[:, k0:]**2).sum(0); out['lam'] = lam / lam.sum()
    except Exception:
        out['lam'] = None
    F.close(); return out


def main(files):
    res = [eig(f) for f in files]; rows = []
    fig, ax = plt.subplots(1, 2, figsize=(14, 5)); nb = len(res); wd = 0.8 / nb
    for j, r in enumerate(res):
        k = np.arange(1, len(r['frac']) + 1)
        ax[0].bar(k + (j - (nb - 1) / 2) * wd, r['frac'], wd, label=f"{r['tag']}  (unexplained {r['unexpl']:.2%})")
        ax[1].plot(k, r['R2'], 'o-', label=r['tag'])
        for i in range(len(k)):
            rows.append(dict(basis=r['tag'], k=k[i], var_fraction=r['frac'][i], R2_cum=r['R2'][i], lam_fraction=None if r['lam'] is None else r['lam'][i]))
        print(r['tag'], 'variance fraction per vector:', ' '.join(f'{x:.3f}' for x in r['frac']), f'| unexplained {r["unexpl"]:.4f}')
        if r['lam'] is not None:
            print('   from coefficients:', ' '.join(f'{x:.3f}' for x in r['lam']))
    ax[0].set_xlabel('component k'); ax[0].set_ylabel('variance fraction  $\\lambda_k / \\sum \\lambda$  (of the variance about $e_0$)'); ax[0].set_yscale('log'); ax[0].legend(fontsize=8)
    ax[0].set_xticks(np.arange(1, 11))
    ax[1].set_xlabel('number of vectors'); ax[1].set_ylabel('cumulative $R^2$ about $e_0$'); ax[1].set_ylim(0.5, 1.0); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3); ax[1].set_xticks(np.arange(1, 11))
    fig.suptitle('Galaxy eigenspectra: importance of each component'); fig.tight_layout()
    tags = '_'.join(r['tag'].split('_')[-1] for r in res)
    pn = os.path.join(ROOT, 'plots_luminosity', f'gal_eigenvalues_{tags}.png'); fig.savefig(pn, dpi=110); print('wrote', pn)
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, 'csv_files', 'gal_eigenvalues.csv'), index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
