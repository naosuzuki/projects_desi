#!/usr/bin/env python3
"""BAL classification from the eigenspectrum fits (16_eigenfit.py) against the pass-Q classes: tables and a figure
(Nao 2026-11-04: 'identify all of the BALs' before the BAL / non-BAL eigenspectra).
Usage: python3 18_bal_classes.py [--fit v3/qso_eigenfit_v6.parquet]
Writes csv_files/qso_bal_classes_<tag>.csv (crosstab and fractions) and plots_luminosity/plots_qso_bal_classes_<tag>.png
(AI and BI distributions, fraction vs S/N and redshift, class agreement)."""
import os, sys, argparse, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, '..'); PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
NAMES = {-1: 'no coverage', 0: 'non-BAL', 1: 'AI > 0', 2: 'BI > 0', 3: 'LoBAL'}

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', default=f'{PROD}/v3/qso_eigenfit_v6.parquet'); ap.add_argument('--tag', default='v6')
    a = ap.parse_args()
    d = pd.read_parquet(a.fit)
    # the two-class flag for the v7 split: BI(C IV) or BI(Si IV) > 0, or a Mg II / Al III trough with AI >= 1000 km/s (LoBAL);
    # the AI > 0 class against the eigenspectrum continuum is noise-driven (half of the z > 2 quasars, tracking the fit chi2)
    g = lambda c, t=0: d[c].fillna(0) > t
    d['BAL_TWO'] = (g('BAL_BI_CIV') | g('BAL_BI_SIIV') | g('BAL_AI_MGII', 999) | g('BAL_AI_ALIII', 999)).astype(np.int8)
    d.to_parquet(a.fit, index=False); print(f'BAL_TWO written to {a.fit}: {int(d.BAL_TWO.sum()):,} BAL of {len(d):,}')
    d = d[d.Z.notna()]
    d['Q'] = d.BAL_CLASS_Q.fillna(-1).astype(int); d['E'] = d.BAL_CLASS_E.astype(int)
    ct = pd.crosstab(d.Q, d.E, rownames=['pass Q'], colnames=['eigen'], margins=True); print(ct)
    cov = d[d.E >= 0]
    civ = np.isfinite(d.BAL_AI_CIV)
    print('BAL_TWO fraction: C IV covered %.4f, not covered %.4f; agreement with pass-Q BI/LoBAL (Q>=2): %d of %d' % (d.BAL_TWO[civ].mean(), d.BAL_TWO[~civ].mean(), int(((d.BAL_TWO == 1) & (d.Q >= 2)).sum()), int((d.Q >= 2).sum())))
    rows = [dict(sample='all', n=len(cov), **{f'E{k}': int((cov.E == k).sum()) for k in range(4)}, bi_frac=float((cov.E == 2).mean()), ai_frac=float((cov.E >= 1).mean()))]
    for sv in ['desi', 'dr17', 'dr8']:
        s = cov[cov.SURVEY == sv]; rows.append(dict(sample=sv, n=len(s), **{f'E{k}': int((s.E == k).sum()) for k in range(4)}, bi_frac=float((s.E == 2).mean()), ai_frac=float((s.E >= 1).mean())))
    hi = cov[np.isfinite(cov.BAL_AI_CIV)]          # C IV covered
    for lo, hi_ in [(1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 5.0)]:
        s = hi[(hi.Z >= lo) & (hi.Z < hi_)]; rows.append(dict(sample=f'z {lo}-{hi_}', n=len(s), **{f'E{k}': int((s.E == k).sum()) for k in range(4)}, bi_frac=float((s.E == 2).mean()) if len(s) else np.nan, ai_frac=float((s.E >= 1).mean()) if len(s) else np.nan))
    t = pd.DataFrame(rows); print(t.round(4).to_string(index=False))
    os.makedirs(f'{R}/csv_files', exist_ok=True); t.to_csv(f'{R}/csv_files/qso_bal_classes_{a.tag}.csv', index=False); ct.to_csv(f'{R}/csv_files/qso_bal_crosstab_{a.tag}.csv')
    # figure
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for k, (col, lab) in enumerate([('BAL_AI_CIV', 'AI(C IV) [km/s]'), ('BAL_BI_CIV', 'BI(C IV) [km/s]'), ('BAL_AI_MGII', 'AI(Mg II) [km/s]')]):
        v = d[col].to_numpy(float); v = v[np.isfinite(v) & (v > 0)]
        ax[0, k].hist(np.log10(v), bins=60, color='C0', alpha=0.7); ax[0, k].set_xlabel(f'log10 {lab}'); ax[0, k].set_ylabel('quasars'); ax[0, k].set_title(f'{len(v):,} with {lab.split()[0]} > 0')
    hi = d[np.isfinite(d.BAL_AI_CIV)]
    zb = pd.cut(hi.Z, [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0]); g = hi.groupby(zb, observed=True)
    ax[1, 0].plot([iv.mid for iv in g.E.mean().index], g.apply(lambda s: (s.E >= 1).mean()), 'o-', label='AI > 0 (eigen)')
    ax[1, 0].plot([iv.mid for iv in g.E.mean().index], g.apply(lambda s: (s.E == 2).mean()), 's-', label='BI > 0 (eigen)')
    ax[1, 0].plot([iv.mid for iv in g.E.mean().index], g.apply(lambda s: (s.Q >= 1).mean()), 'o--', color='0.5', label='AI > 0 (pass Q)')
    ax[1, 0].plot([iv.mid for iv in g.E.mean().index], g.apply(lambda s: (s.Q == 2).mean()), 's--', color='0.5', label='BI > 0 (pass Q)')
    ax[1, 0].set_xlabel('z'); ax[1, 0].set_ylabel('fraction, C IV covered'); ax[1, 0].legend(fontsize=8)
    cb = pd.cut(hi.CHI2 / hi.NPIX.clip(lower=1), [0, 0.8, 1.0, 1.2, 1.5, 2, 3, 5, 100]); g = hi.groupby(cb, observed=True)
    ax[1, 1].plot(range(len(g)), g.apply(lambda s: (s.E >= 1).mean()), 'o-', label='AI > 0'); ax[1, 1].plot(range(len(g)), g.apply(lambda s: (s.E == 2).mean()), 's-', label='BI > 0')
    ax[1, 1].set_xticks(range(len(g))); ax[1, 1].set_xticklabels([str(i) for i in g.size().index], rotation=45, fontsize=7); ax[1, 1].set_xlabel('reduced chi2 of the fit'); ax[1, 1].set_ylabel('fraction'); ax[1, 1].legend(fontsize=8)
    m = ct.drop('All', axis=0).drop('All', axis=1); im = ax[1, 2].imshow(np.log10(m.to_numpy().astype(float) + 1), cmap='Blues')
    ax[1, 2].set_xticks(range(m.shape[1])); ax[1, 2].set_xticklabels([NAMES.get(c, c) for c in m.columns], rotation=45, fontsize=7); ax[1, 2].set_yticks(range(m.shape[0])); ax[1, 2].set_yticklabels([NAMES.get(c, c) for c in m.index], fontsize=7)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]): ax[1, 2].text(j, i, f'{m.iloc[i, j]:,}', ha='center', va='center', fontsize=6)
    ax[1, 2].set_xlabel('eigenspectrum continuum'); ax[1, 2].set_ylabel('pass Q (power law)'); ax[1, 2].set_title('class agreement')
    fig.suptitle(f'BAL troughs against the {a.tag} eigenspectrum fit: {len(d):,} training quasars'); fig.tight_layout()
    out = f'{R}/plots_luminosity/plots_qso_bal_classes_{a.tag}.png'; fig.savefig(out, dpi=100); print(out)
