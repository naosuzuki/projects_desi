"""Redshift distribution of the specphot sample (21_specphot_fit.py, non-BAL + BAL tables), stacked by survey.
Outputs plots_luminosity/plots_qso_specphot_zhist_<tag>.png and csv_files/plots_qso_specphot_zhist_<tag>.csv."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='v1'); a = ap.parse_args()
d = pd.concat([pd.read_parquet(f'{PROD}/specphot/qso_specphot_{k}_{a.tag}.parquet', columns=['SURVEY', 'Z']) for k in ('coeff', 'bal')], ignore_index=True)
order = [('dr8', 'SDSS DR8', 'tab:green'), ('dr17', 'SDSS DR17', 'tab:orange'), ('desi', 'DESI DR2', 'tab:blue')]
edges = np.arange(0, 5.01, 0.05); mid = 0.5 * (edges[:-1] + edges[1:])
hist = {k: np.histogram(d.Z[d.SURVEY == k], edges)[0] for k, _, _ in order}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for axi, log in zip(ax, (False, True)):
    axi.hist([d.Z[d.SURVEY == k] for k, _, _ in order], edges, stacked=True, color=[c for _, _, c in order],
             label=[f'{lab}: {len(d[d.SURVEY == k]):,}, median z {d.Z[d.SURVEY == k].median():.2f}' for k, lab, _ in order])
    axi.step(edges[:-1], sum(hist.values()), where='post', color='k', lw=0.8, label=f'all: {len(d):,}, median z {d.Z.median():.2f}')
    axi.set_xlabel('redshift'); axi.set_ylabel('quasars per 0.05'); axi.set_xlim(0, 5)
    if log:
        axi.set_yscale('log'); axi.set_ylim(1, None)
    axi.legend(fontsize=9)
ax[0].set_title('Spectrophotometric quasar sample, stacked by survey'); ax[1].set_title('same, log scale')
fig.tight_layout(); fig.savefig(f'{REPO}/plots_luminosity/plots_qso_specphot_zhist_{a.tag}.png', dpi=130)
pd.DataFrame({'z_lo': edges[:-1], 'z_hi': edges[1:], **{f'n_{k}': hist[k] for k, _, _ in order}, 'n_all': sum(hist.values())}).to_csv(f'{REPO}/csv_files/plots_qso_specphot_zhist_{a.tag}.csv', index=False)
print('done')
