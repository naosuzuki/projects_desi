"""Figure of the brightness-prediction benchmark (26/28_specphot_ml*.py): left, held-out rms per method for the
coefficient-only inputs (A: c_1..c_9; A2: + general-basis g_k; D: direction only) against the controls (N noise proxies,
S shape measures, B all, C + calibration) and the pixel-level models; right, predicted against true M_AB for the best
coefficient-only model (from the saved held-out predictions, if present). Outputs plots_luminosity/plots_qso_ml_bench_<slice>_<tag>.png."""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--slice', default='1.65-1.7'); ap.add_argument('--tag', default='v1'); ap.add_argument('--target', default='M'); a = ap.parse_args()
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f'{REPO}/csv_files/qso_ml_bench_z{a.slice}_{a.tag}*.csv'))], ignore_index=True)   # main csv + per-round files (--out)
pix = f'{REPO}/csv_files/qso_ml_pix_z{a.slice}_{a.tag}.csv'
if os.path.exists(pix):
    d = pd.concat([d[~d.method.str.startswith('pix:')], pd.read_csv(pix)], ignore_index=True)
d = d[d.error.isna() | (d.error == '')]
d = d[(d.target == a.target) & (d['round'].astype(str) != '3')] if 'round' in d else d[d.target == a.target]
sigma = d.sigma.iloc[0]
order = {'A': 0, 'A2': 1, 'D': 2, 'N': 3, 'S': 4, 'B': 5, 'C': 6}; col = {'A': 'tab:red', 'A2': 'tab:orange', 'D': 'tab:pink', 'N': 'grey', 'S': 'tab:green', 'B': 'tab:blue', 'C': 'tab:purple'}
lab = {'A': 'A: $c_1..c_9$', 'A2': 'A2: $c_k$ + general-basis $g_k$', 'D': 'D: direction $c/|c|$ only', 'N': 'N: noise proxies only (control)', 'S': 'S: shape measures', 'B': 'B: all fit outputs', 'C': 'C: + calibration record'}
coef = d[d.features.isin(order)].sort_values('rms').drop_duplicates(['method', 'features']).copy()     # the same method/input may appear in several rounds: keep the best
pix = d[d.method.str.startswith('pix')].sort_values('rms').drop_duplicates(['method', 'features']).copy()
best = coef[coef.features == 'A'].sort_values('rms').iloc[0] if (coef.features == 'A').any() else None
meths = coef.groupby('method').rms.min().sort_values(ascending=False).index.tolist()
fig = plt.figure(figsize=(15, 0.42 * len(meths) + 4)); gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])
ax = fig.add_subplot(gs[0]); y = np.arange(len(meths)); h = 0.11
for f, k in order.items():
    sub = coef[coef.features == f].set_index('method').reindex(meths)
    ax.barh(y + (k - 3) * h, sub.rms, h, color=col[f], alpha=0.9 if f in ('A', 'A2', 'D') else 0.45, label=lab[f])
ax.axvline(sigma, color='k', ls='--', lw=1, label=f'no prediction: $\\sigma(M)$ = {sigma:.2f}')
ax.set_yticks(y); ax.set_yticklabels(meths); ax.set_xlabel(f'held-out rms of $M_{{\\rm AB}}(3025)$ [mag], target {a.target}'); ax.set_xlim(0.3, 1.15); ax.grid(axis='x', alpha=0.3)
ax.set_title(f'Brightness from the eigencoefficients, $z={a.slice}$: methods $\\times$ inputs'); ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=3)
if len(pix):
    for i, r in pix[pix.method.str.contains('cnn|pca_hgb')].sort_values('rms', ascending=False).reset_index().iterrows():
        ax.text(1.14, -1.2 - 0.45 * i, f"{r.method} [{r.features}]: rms {r.rms:.2f}, $R^2$ {r.r2:.2f}", ha='right', va='top', fontsize=7.5, color='tab:brown')
    ax.set_ylim(-1.5 - 0.45 * pix.method.str.contains('cnn|pca_hgb').sum(), len(meths) - 0.5)
ax2 = fig.add_subplot(gs[1])
import re
def parse(f):                                                             # qso_ml_pred_<method>_<fset>_<target>_r<round>.npz
    m = re.match(r'qso_ml_pred_(.+)_(A2?|AC|XC?)_(M(?:_raw|_rel|_rel_raw)?)_r(.+)\.npz$', os.path.basename(f)); return m.groups() if m else None
pf = [f for f in sorted(glob.glob(f'{REPO}/csv_files/qso_ml_pred_*.npz')) if parse(f) and parse(f)[2] == a.target and parse(f)[1] in ('A', 'A2')]
pf = sorted(pf, key=lambda f: (parse(f)[1] == 'A2', parse(f)[3] == '3p'))  # prefer the 18 coefficients, then the round-3 rerun
if pf and best is not None:
    p = pf[-1]; z = np.load(p); M, P = z['M'], z['pred']
    ax2.hist2d(M, P, bins=[np.linspace(-27.5, -20.5, 90)] * 2, cmap='Greys', norm=matplotlib.colors.LogNorm()); ax2.plot([-27.5, -20.5], [-27.5, -20.5], 'r-', lw=1)
    e = np.arange(-27.5, -20.4, 0.25); ib = np.digitize(M, e) - 1; med = [np.median(P[ib == j]) if (ib == j).sum() > 50 else np.nan for j in range(len(e) - 1)]
    ax2.plot(0.5 * (e[:-1] + e[1:]), med, color='tab:red', lw=2, label='median prediction'); ax2.invert_xaxis(); ax2.invert_yaxis()
    ax2.set_xlabel('true $M_{\\rm AB}(3025)$'); ax2.set_ylabel('predicted $M_{\\rm AB}(3025)$'); ax2.legend(fontsize=8)
    r2p = 1 - np.mean((M - P) ** 2) / np.var(M); ax2.set_title(f'{parse(p)[0].replace("_", " ")}, inputs {parse(p)[1]}, round {parse(p)[3]}: rms {np.sqrt(np.mean((M - P) ** 2)):.2f}, $R^2$ {r2p:.2f}\n({os.path.basename(p)})', fontsize=10)
else:
    ax2.text(0.5, 0.5, 'held-out predictions not saved yet', ha='center', va='center', transform=ax2.transAxes); ax2.set_axis_off()
fig.tight_layout(rect=(0, 0.04, 1, 1)); out = f'{REPO}/plots_luminosity/plots_qso_ml_bench_z{a.slice}_{a.tag}.png'; fig.savefig(out, dpi=120); print(out)

# ---- LaTeX table: best held-out result per (input set, method family) for the paper ----
lab = {'A': 'nine $c_k$', 'A2': '$c_k$ + $g_k$ (18)', 'AC': '$c_k$ + colour', 'X': 'colour-free $c_k$', 'XC': 'colour-free $c_k$ + colour', 'D': 'direction $c_k/|c|$', 'N': 'noise proxies only', 'S': 'coefficients + survey/BAL/range',
       'B': 'S + noise proxies', 'C': 'B + calibration record'}
trn = lambda r: (str(r.train) if isinstance(r.train, str) and r.train != 'nan' else a.slice).replace('-', '--')  # training redshift range
dd = d[d.target == a.target].copy(); dd['r8'] = dd['round'].astype(str).str.startswith('8')   # round 8 ran on the DESI quasars of the slice only: own block
rows = []
def block(sub, fs, title):
    best = sub.sort_values('rms').drop_duplicates(['method', 'features'])
    rows.append(f"\\multicolumn{{7}}{{l}}{{\\textit{{{title}}}}}\\\\")
    for f in fs:
        sb = best[best.features == f].sort_values('rms')
        if not len(sb):
            continue
        keep = pd.concat([sb[sb.method.str.startswith(('OLS', 'ridge'))].head(1), sb[sb.method.str.match('HGB')].head(1), sb[sb.method.str.match('kNN')].head(1), sb[sb.method.str.contains('MLP|stack')].head(1)]).drop_duplicates()
        for i, r in enumerate(keep.sort_values('rms', ascending=False).itertuples()):
            rows.append(f"{lab.get(f, f) if i == 0 else ''} & {r.method} & {trn(r)} & {r.r2:.3f} & {r.rms:.3f} & {r.rms_bright:.3f} & {r.mad68:.3f}\\\\")
    return best
main = dd[~dd.r8 & ~dd.method.str.startswith('pix:')]; n0 = int(main.n.max()) if len(main) else 0
block(main, ['A', 'A2', 'D', 'N', 'S', 'B', 'C'], f"all surveys, {n0:,} quasars, $\\sigma(M)={main.sigma.max():.2f}$")
pix = dd[dd.method.str.startswith('pix:') & dd.method.str.contains('cnn|pca_hgb')].sort_values('rms').drop_duplicates(['method', 'features'])
for r in pix.itertuples():
    rows.append(f"pixels ({r.features}) & {r.method[5:]} & {trn(r)} & {r.r2:.3f} & {r.rms:.3f} & {r.rms_bright:.3f} & {r.mad68:.3f}\\\\")
r8 = dd[dd.r8]
if len(r8):
    block(r8, ['A', 'AC', 'X', 'XC'], f"round 8, DESI only, {int(r8.n.max()):,} quasars, $\\sigma(M)={r8.sigma.max():.2f}$")
tex = f'{REPO}/csv_files/qso_ml_bench_z{a.slice}_{a.tag}_{a.target}.tex'; open(tex, 'w').write('\n'.join(rows) + '\n'); print(tex)
