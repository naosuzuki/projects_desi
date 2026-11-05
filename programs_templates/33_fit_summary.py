#!/usr/bin/env python3
"""Summary of the per-galaxy fit tables of 32_fit_galaxies.py (Nao 2026-11-04): chi2/NPIX of the full basis and of e0 alone by
program, redshift, S/N, Dn4000 bin (DB_EFF) and line-EW bin (EB), and the coefficient distributions C_k/C0 by EB (the
coefficient-based selection of non-emission galaxies). Joins the pass tables for Z, DB_EFF, EB, MABS_G, SNR_SDSSPIX.
Output: csv_files/gal_fit_<tag>_percentiles.csv (one row per group and quantity), csv_files/gal_fit_<tag>_sample.parquet
(1 in 20 galaxies, all columns, on exdisk3 via ~/data), plots_luminosity/gal_fit_<tag>_chi2.png.
Usage: python3 33_fit_summary.py --fit ~/data/sdssdesi_galfit_p1 --tag p1 [--workers 6]
"""
import os, glob, argparse, time, json
import numpy as np, pandas as pd
from multiprocessing import Pool

TABLES = os.path.expanduser('~/data/sdssdesi_passC_v3t/tables'); STATUS = os.path.expanduser('~/data/sdssdesi_status')
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
KEEP = ['Z', 'DB_EFF', 'EB', 'CB', 'MABS_G', 'SNR_SDSSPIX', 'DN4000', 'EW', 'ZWARN', 'DELTACHI2']


def load(args):
    prog, fn = args
    g = pd.read_parquet(fn).drop(columns=['Z'])
    t = pd.read_parquet(os.path.join(TABLES, prog, os.path.basename(fn)), columns=['TARGETID'] + KEEP).drop_duplicates('TARGETID').set_index('TARGETID')
    g = g.join(t, on='TARGETID'); g['PROG'] = prog
    g['RCHI2'] = g.CHI2 / np.maximum(g.NPIX, 1); g['RCHI2_E0'] = g.CHI2_E0 / np.maximum(g.NPIX, 1)
    for k in range(1, 11):
        g[f'R{k}'] = g[f'C{k}'] / g.C0
    g['ZBIN'] = np.digitize(g.Z, [0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4]); g['SBIN'] = np.digitize(g.SNR, [1, 2, 3, 5, 10, 20])
    cols = ['PROG', 'ZBIN', 'SBIN', 'DB_EFF', 'EB', 'RCHI2', 'RCHI2_E0', 'C0'] + [f'R{k}' for k in range(1, 11)]
    agg = g[cols].copy(); agg['N'] = 1
    return g.iloc[::20], agg


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', required=True); ap.add_argument('--tag', default='p1'); ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args(); t0 = time.time()
    jobs = [(p, f) for p in ['bright', 'dark'] for f in sorted(glob.glob(os.path.join(a.fit, p, 'hpx*.parquet')))]
    sj = os.path.join(STATUS, f'step_gal_fitsum_{a.tag}.json')
    json.dump(dict(done=0, total=len(jobs), status='running', started=t0, note='reading fit tables'), open(sj, 'w'))
    samp, aggs = [], []
    with Pool(a.workers) as pool:
        for i, (s, g) in enumerate(pool.imap_unordered(load, jobs, chunksize=4)):
            samp.append(s); aggs.append(g)
            if i % 500 == 0:
                json.dump(dict(done=i, total=len(jobs), status='running', started=t0, note=f'{time.time()-t0:.0f}s'), open(sj, 'w')); print(i, flush=True)
    S = pd.concat(samp, ignore_index=True); A = pd.concat(aggs, ignore_index=True); del aggs
    print(f'{len(A):,} galaxies read in {time.time()-t0:.0f}s; sample {len(S):,}', flush=True)
    S.to_parquet(os.path.expanduser(f'~/data/sdssdesi_galfit_{a.tag}/gal_fit_{a.tag}_sample.parquet'))
    q = [0.05, 0.16, 0.5, 0.84, 0.95]; rows = []
    def add(name, key, d):
        for col in ['RCHI2', 'RCHI2_E0', 'C0'] + [f'R{k}' for k in range(1, 11)]:
            v = d[col].to_numpy(float); v = v[np.isfinite(v)]
            if len(v):
                rows.append(dict(group=name, key=str(key), n=len(v), quantity=col, **{f'p{int(x*100)}': np.quantile(v, x) for x in q}))
    add('all', 'all', A)
    for grp in ['PROG', 'ZBIN', 'SBIN', 'DB_EFF', 'EB']:
        for key, d in A.groupby(grp):
            add(grp, key, d)
    for key, d in A.groupby(['DB_EFF', 'EB']):
        add('DB_EFF_EB', key, d)
    R = pd.DataFrame(rows); fn = os.path.join(ROOT, 'csv_files', f'gal_fit_{a.tag}_percentiles.csv'); R.to_csv(fn, index=False)
    print('wrote', fn, flush=True)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for axx, grp, lab in zip(ax.ravel(), ['ZBIN', 'SBIN', 'EB', 'DB_EFF', 'PROG'], ['z bin (edges 0.1 .. 1.4)', 'S/N bin (1,2,3,5,10,20)', 'line EW bin', 'Dn4000 bin', 'program']):
        d = R[(R.group == grp)]
        for col, c in [('RCHI2', 'k'), ('RCHI2_E0', 'r')]:
            e = d[d.quantity == col]; x = np.arange(len(e))
            axx.errorbar(x, e.p50, yerr=[e.p50 - e.p16, e.p84 - e.p50], fmt='o-', color=c, label=col, capsize=3)
        axx.set_xticks(np.arange(len(e))); axx.set_xticklabels(e.key); axx.set_xlabel(lab); axx.set_ylabel('chi2 / NPIX (median, 16-84 %)'); axx.set_yscale('log'); axx.legend(fontsize=8)
    axx = ax[1, 2]
    for eb, d in S.groupby('EB'):
        axx.hist(np.clip(d.R1, -200, 200), bins=100, range=(-200, 200), histtype='step', label=f'EB {eb}', density=True)
    axx.set_xlabel('C1 / C0 (star-formation vector)'); axx.legend(fontsize=7); axx.set_yscale('log')
    fig.suptitle(f'Fit of basis {a.tag} to every DESI galaxy: {len(A):,} galaxies'); fig.tight_layout()
    pn = os.path.join(ROOT, 'plots_luminosity', f'gal_fit_{a.tag}_chi2.png'); fig.savefig(pn, dpi=110); print('wrote', pn, flush=True)
    json.dump(dict(done=len(jobs), total=len(jobs), status='done', started=t0, note=f'{len(A):,} galaxies, {time.time()-t0:.0f}s'), open(sj, 'w'))


if __name__ == '__main__':
    main()
