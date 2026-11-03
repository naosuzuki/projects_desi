#!/usr/bin/env python3
"""Per-quasar normalization record and window fluxes (Nao 2026-10-31), from the v2 rest-frame store (one pass).

For every stacked quasar: NORM (the factor its spectrum was divided by = ivar-weighted mean of F_lam*(1+z), cgs,
in its own window), NORM_SRC (0 window mean, 1 power-law fallback, 2 none), and the mean of the NORMALIZED spectrum in each of
the three windows (W5115 = 5100-5130, W3025 = 3000-3050, W1282 = 1275-1290 A) with the good-pixel fraction COV*.
The window it was normalized in reads 1; another covered window gives the flux ratio directly, e.g. a z<0.8 quasar with
COV3025 >= 0.8 has R1 = F(3025)/F(5115) = W3025.
Joined with qso_cells_all (Z, SURVEY, LOGFABS_NORM, LBIN, SBIN, ALPHA_LAM ...) -> quasar/v2/qso_norm_table.parquet
Usage: python3 09_norm_table.py [--store DIR] [--nproc 4]
"""
import os, glob, argparse, numpy as np, pandas as pd
from multiprocessing import Pool

P = os.path.expanduser('~/data/sdssdesi_products/quasar')
WINDOWS = {'5115': (5100.0, 5130.0), '3025': (3000.0, 3050.0), '1282': (1275.0, 1290.0)}
DLOG = 1e-4


def one(path):
    d = np.load(path)
    q0 = int(d['QRID_MIN']); start = d['START']; off = d['OFFSET']; F = d['FLUX'].astype(np.float32); I = d['IVAR'].astype(np.float32)
    n = len(d['ID']); out = {'ID': d['ID'], 'CELL_I': d['CELL_I'], 'NORM': d['NORM'], 'NORM_SRC': d['NORM_SRC']}
    for k, (lo, hi) in WINDOWS.items():
        r0 = int(np.ceil(np.log10(lo) / DLOG)); r1 = int(np.floor(np.log10(hi) / DLOG)); npx = r1 - r0 + 1
        m = np.full(n, np.nan, np.float32); cov = np.zeros(n, np.float32)
        for i in range(n):
            a = off[i] + (r0 - q0 - start[i]); b = a + npx; L = off[i + 1] - off[i]
            ja, jb = max(a, off[i]), min(b, off[i] + L)
            if jb <= ja:
                continue
            f = F[ja:jb]; iv = I[ja:jb]; g = (iv > 0) & np.isfinite(f) & np.isfinite(iv)
            cov[i] = g.sum() / npx
            if g.any():
                m[i] = (f[g] * iv[g]).sum() / iv[g].sum()
        out['W' + k] = m; out['COV' + k] = cov
    out['STORE_SURVEY'] = np.full(n, os.path.basename(os.path.dirname(path)))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--store', default='/Volumes/exdisk3/sdssdesi_qso_restframe_v2')
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    files = sorted(glob.glob(f'{a.store}/*/*.npz'))
    print(len(files), 'store files', flush=True)
    with Pool(a.nproc) as p:
        parts = []
        for k, df in enumerate(p.imap_unordered(one, files, chunksize=20)):
            parts.append(df)
            if k % 2000 == 0:
                print(k, flush=True)
    t = pd.concat(parts, ignore_index=True)
    # join: SDSS ID = (PLATE*100000 + MJD)*10000 + FIBERID; DESI ID = TARGETID
    c = pd.read_parquet(f'{P}/qso_cells_all.parquet')
    c = c[c.TRAIN & (c.STACK_CELL != '')]
    sd = c.SURVEY.isin(['dr8', 'dr17'])
    c['ID'] = np.where(sd, (c.PLATE.fillna(0).astype(np.int64) * 100000 + c.MJD.fillna(0).astype(np.int64)) * 10000 + c.FIBERID.fillna(0).astype(np.int64),
                       c.TARGETID.fillna(0).astype(np.int64))
    keep = [x for x in ['ID', 'SURVEY', 'PLATE', 'MJD', 'FIBERID', 'TARGETID', 'PROGRAM', 'Z', 'LOGFABS_NORM', 'LOGL_NORM', 'ALPHA_LAM',
                        'LBIN', 'SBIN', 'EBIN', 'ZBIN', 'STACK_CELL', 'SNR_ALL'] if x in c.columns]
    c['_S'] = np.where(sd, c.SURVEY, 'desi')
    t = t.merge(c[keep + ['_S']].rename(columns={'_S': 'STORE_SURVEY'}), on=['ID', 'STORE_SURVEY'], how='left')
    out = f'{P}/v2/qso_norm_table.parquet'
    t.to_parquet(out, index=False)
    print(f'{len(t):,} quasars ({t.Z.notna().sum():,} matched to qso_cells_all) -> {out}')


if __name__ == '__main__':
    main()
