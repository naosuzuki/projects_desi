#!/usr/bin/env python3
"""Stellar EM-PCA with sqrt(N) stack weights against equal weights, both on the telluric-corrected stacks (Nao 2026-11-03; step7_empca.py --stack-weight sqrtn, as programs_quasars/06_empca.py).
summary  per class: pixels, R2, overlap of the first n vectors, pixel noise of the vectors, share of the weight in the largest tenth of the stacks.
byn      10-vector residual of the stacks in thirds of N: rms of (stack - fit) / mean over the basis pixels, and the same divided by the stack noise.
amp      telluric amplitude held by each basis, heldout residuals in the bands and controls: the tables of 06_star_compare.py with v2 = equal, tell = sqrt(N).
Output: csv_files/plots_star_sqrtn_{summary,byn,amp,heldout}.csv
"""
import os, sys, importlib.util
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tellapply
sp = importlib.util.spec_from_file_location('compare', os.path.join(HERE, '06_star_compare.py')); C = importlib.util.module_from_spec(sp); sp.loader.exec_module(C)
SFX = '_sqrtn'


def summary(s, c, st, a, b):
    k = a['ok'] & b['ok']; E, F = b['EIGVEC'][:, k], a['EIGVEC'][:, k]; noise = lambda e: 1.4826 * np.median(np.abs(np.diff(e, n=2, axis=1)), axis=1)
    i = np.array([np.where(st['name'] == x)[0][0] for x in a['id'][a['kind'] == 'stack']]); n = st['n'][i].max(axis=1); w = np.sqrt(n / np.median(n)); top = np.sort(w)[-max(1, len(w) // 10):].sum() / w.sum(); r = C.space(s, c, a, b)
    return dict(survey=s, cls=c, stacks=len(i), n_min=int(n.min()), n_med=int(np.median(n)), n_max=int(n.max()), top10pct_share_equal=max(1, len(w) // 10) / len(w), top10pct_share_sqrtn=top, pix_equal=int(b['ok'].sum()), pix_sqrtn=int(a['ok'].sum()), **{f'r2{m}_equal': b['VARIANCE'][m - 1] for m in (1, 3, 10)}, **{f'r2{m}_sqrtn': a['VARIANCE'][m - 1] for m in (1, 3, 10)}, **{f'mincos{m}': r[f'mincos{m}'] for m in (1, 3, 5, 10)}, noise_ratio_med=np.median(noise(F) / noise(E)), mean_change_rms=np.sqrt(np.mean((a['MEAN'][k] / b['MEAN'][k] - 1)**2)))


def by_n(s, c, st, a, b):
    i = np.array([np.where(st['name'] == x)[0][0] for x in a['id'][a['kind'] == 'stack']]); N = st['n'][i].max(axis=1); third = np.searchsorted(np.quantile(N, [1 / 3, 2 / 3]), N); rows = []
    for lab, q in (('equal', b), ('sqrtn', a)):
        nm = q['norm'][q['kind'] == 'stack']; ok = q['ok']; res = []; chi = []
        for j, m in zip(i, nm):
            g = ok & (st['n'][j] > 0) & np.isfinite(st['f'][j]) & (st['e'][j] > 0); x = st['f'][j] / m; w = (m / np.where(g, st['e'][j], 1))**2; r0 = x - q['MEAN']; E = q['EIGVEC']; r = r0 - E.T @ C.lsq(E[:, g], r0[g], w[g]); u = g & (np.abs(q['MEAN']) > 1e-3) & (np.abs(r / np.where(q['MEAN'] != 0, q['MEAN'], 1)) < C.RCLIP)
            res.append(np.sqrt(np.mean((r[u] / q['MEAN'][u])**2))); chi.append(np.sqrt(np.mean(r[u]**2 * w[u])))
        res, chi = np.array(res), np.array(chi)
        for t in range(3): rows.append(dict(survey=s, cls=c, third=t, n_med=np.median(N[third == t]), basis=lab, rms=np.median(res[third == t]), over_noise=np.median(chi[third == t])))
    return rows


def main():
    pd.set_option('display.width', 200); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200)
    B8, Bd = tellapply.basis('dr8'), tellapply.basis('desi'); Bs = {'sdss': tellapply.Basis('dr8', shift=C.SHIFT), 'desi': tellapply.Basis('desi', shift=C.SHIFT)}; su = []; bn = []; am = []; he = []
    for x in Bs.values():
        i = np.where(x.ok)[0]; x.Ffull[0, i] = np.roll(x.Ffull[0, i], C.SHIFT)
    for s in C.CLASSES:
        st = C.stacks(C.TELL, s); ind = C.individuals(s)[0]; bas = Bd if s == 'desi' else B8
        for c in C.CLASSES[s]:
            a, b = C.base(C.TELL, s, c, SFX), C.base(C.TELL, s, c); su.append(summary(s, c, st, a, b)); bn += by_n(s, c, st, a, b); am += C.amp(s, c, a, b, bas, Bs[s]); he += C.heldout(s, c, a, b, st, st, ind, ind, bas); print(' ', s, c, flush=True)
    su, bn, am, he = pd.DataFrame(su), pd.DataFrame(bn), pd.DataFrame(am), pd.DataFrame(he)
    for n, t in (('summary', su), ('byn', bn), ('amp', am), ('heldout', he)): t.to_csv(f'{C.CSV}/plots_star_sqrtn_{n}.csv', index=False)
    print(su.round(4).T.to_string()); print(bn.pivot_table(index=['survey', 'cls', 'third', 'n_med'], columns='basis', values=['rms', 'over_noise']).round(4))
    print(am[(am.kind == 'star') & (am.control == 'real')].pivot_table(index=['survey', 'T'], columns='basis', values='rms', aggfunc='median').round(3)); print(he[he.kind == 'stack'].pivot_table(index=['survey', 'cls', 'window'], columns='basis', values='rms').round(5).loc[lambda t: t.index.get_level_values('window').isin(['all bands', 'control'])])


if __name__ == '__main__':
    main()
