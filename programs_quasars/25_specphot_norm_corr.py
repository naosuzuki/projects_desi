"""How well do the eigencoefficients predict the brightness? (Nao 2026-11-12). In a redshift slice of the non-BAL
coefficient table (photometric correction, fit, NORM > 0): per coefficient the Pearson and Spearman correlation with
M_AB(3025) (NORM as absolute AB magnitude, each quasar's own distance modulus), the slope of M against c_k alone and the
scatter it leaves; then the joint linear prediction M = a + sum_k b_k c_k with 5-fold cross-validated R^2 and rms, built
up one coefficient at a time in order of gain. Done for the whole slice and for the bright half (M_AB < median), where
the coefficients are not noise-dominated. Outputs csv_files/qso_norm_corr_z<zlo>-<zhi>_<tag>.csv and the LaTeX rows."""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import qsolib

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70); ap.add_argument('--tag', default='v1'); a = ap.parse_args()
cols = ['SURVEY', 'Z', 'NORM', 'PCORR_MODE'] + [f'C{k}' for k in range(1, 10)]
d = pd.read_parquet(f'{PROD}/specphot/qso_specphot_coeff_{a.tag}.parquet', columns=cols)
s = d[(d.Z >= a.zlo) & (d.Z < a.zhi) & (d.NORM > 0) & d.PCORR_MODE.isin([0, 1, 2])].dropna().copy()
K = -48.6 - 2.5 * np.log10(3025.0**2 / 2.998e18); s['MAB'] = -2.5 * (np.log10(s.NORM) + 0.4 * qsolib.passB.COSMO.distmod(s.Z.to_numpy()).value) + K
C = s[[f'C{k}' for k in range(1, 10)]].to_numpy(); M = s.MAB.to_numpy()


def cv_fit(X, y, nfold=5, seed=1):
    """Linear least squares y = a + X b with k-fold cross-validation: returns (R2, rms) on the held-out folds and the full-sample coefficients."""
    rng = np.random.default_rng(seed); fold = rng.integers(0, nfold, len(y)); pred = np.empty_like(y)
    A = np.column_stack([np.ones(len(y)), X])
    for f in range(nfold):
        tr = fold != f; coef = np.linalg.lstsq(A[tr], y[tr], rcond=None)[0]; pred[~tr] = A[~tr] @ coef
    r = y - pred; return 1 - r.var() / y.var(), r.std(), np.linalg.lstsq(A, y, rcond=None)[0]


rows = []; tex = []
for label, sel in (('all', np.ones(len(M), bool)), ('bright half', M < np.median(M))):
    Cs, Ms = C[sel], M[sel]; sig = Ms.std()
    tex.append(f'\\multicolumn{{7}}{{l}}{{\\textit{{{label}}}: {sel.sum():,} quasars, $M_{{\\rm AB}}$ {Ms.min():.1f} to {Ms.max():.1f}, $\\sigma(M)={sig:.2f}$}}\\\\')
    single = []
    for k in range(9):
        pr = stats.pearsonr(Cs[:, k], Ms)[0]; sr = stats.spearmanr(Cs[:, k], Ms)[0]
        r2, rms, coef = cv_fit(Cs[:, [k]], Ms); single.append((k, r2))
        rows.append(dict(slice=f'{a.zlo}-{a.zhi}', sample=label, n=int(sel.sum()), coef=f'C{k + 1}', pearson=pr, spearman=sr, slope_mag_per_unit=coef[1], r2_single=r2, rms_single=rms))
        tex.append(f'$c_{k + 1}$ & {pr:+.3f} & {sr:+.3f} & {1000 * coef[1]:+.1f} & {r2:.3f} & {rms:.2f} & \\\\')
    # greedy build-up of the joint prediction
    left = list(range(9)); order = []; built = []
    while left:
        best = max(left, key=lambda k: cv_fit(Cs[:, order + [k]], Ms)[0]); order.append(best); left.remove(best)
        r2, rms, coef = cv_fit(Cs[:, order], Ms); built.append((best, r2, rms))
        rows.append(dict(slice=f'{a.zlo}-{a.zhi}', sample=label, n=int(sel.sum()), coef='+'.join(f'C{k + 1}' for k in order), pearson=np.nan, spearman=np.nan, slope_mag_per_unit=np.nan, r2_single=r2, rms_single=rms))
    tex.append('joint, added in order & & & & & & ' + ', '.join(f'$c_{k + 1}$ ({r2:.3f}, {rms:.2f})' for k, r2, rms in built) + '\\\\')
    print(label, sel.sum(), 'sigma(M) %.2f' % sig); print('  single R2:', ' '.join(f'c{k + 1}={r2:.3f}' for k, r2 in single)); print('  joint:', ' '.join(f'+c{k + 1}->{r2:.3f}/{rms:.2f}' for k, r2, rms in built))
stem = f'qso_norm_corr_z{a.zlo}-{a.zhi}_{a.tag}'; pd.DataFrame(rows).to_csv(f'{REPO}/csv_files/{stem}.csv', index=False)
open(f'{REPO}/csv_files/{stem}.tex', 'w').write('\n'.join(tex) + '\n'); print('wrote', stem)
