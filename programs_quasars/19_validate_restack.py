#!/usr/bin/env python3
"""Compare the store-based restack of the run-5 cells (17_restack_store.py --cells qso_cells_v3_all) with the run-5 stacks
(v3/qso_stacks.npz): per common stack the median |WMEAN_store / WMEAN_5 - 1| over pixels with N >= 50 and |N_store - N_5|.
Usage: python3 19_validate_restack.py [--a PROD/qso_cells_v3_all/qso_stacks.npz] [--b PROD/v3/qso_stacks.npz]"""
import os, argparse, numpy as np, pandas as pd
PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
ap = argparse.ArgumentParser(); ap.add_argument('--a', default=f'{PROD}/qso_cells_v3_all/qso_stacks.npz'); ap.add_argument('--b', default=f'{PROD}/v3/qso_stacks.npz'); a = ap.parse_args()
A, B = np.load(a.a, allow_pickle=True), np.load(a.b, allow_pickle=True)
ca, cb = list(A['STACK_CELL']), list(B['STACK_CELL']); common = [c for c in cb if c in set(ca)]
print(f'{len(ca)} stacks in A, {len(cb)} in B, {len(common)} common')
rows = []
for c in common:
    i, j = ca.index(c), cb.index(c); na, nb = A['N'][i], B['N'][j]; ok = (na >= 50) & (nb >= 50)
    r = A['WMEAN'][i][ok] / B['WMEAN'][j][ok] - 1
    rows.append(dict(cell=c, nmax_a=int(na.max()), nmax_b=int(nb.max()), dn_max=int(np.abs(na - nb).max()), med_abs_ratio=float(np.nanmedian(np.abs(r))) if ok.any() else np.nan, p99_abs_ratio=float(np.nanpercentile(np.abs(r), 99)) if ok.any() else np.nan))
t = pd.DataFrame(rows); print(t.describe().round(5).to_string())
print('worst 8 by median ratio:'); print(t.sort_values('med_abs_ratio', ascending=False).head(8).to_string(index=False))
print('stacks with |dN| > 0:', int((t.dn_max > 0).sum()))
