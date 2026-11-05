#!/usr/bin/env python3
"""Release v3: chain the stacks of one version (G or 5450; v3lib.py) to unit F_lambda in 5400-5500 A, as 19_chain_stacks.py.
Per stack (cell of that version's grid), from the per-galaxy v3 tables of the members (ACC_<ver> >= 0, NORM_OK):
MABS_G_MED (MABS_G_FULL), M5450_MED (M_AB5450). The stack's type (DB_EFF) is its Dn4000 bin, re-typed by its median line EW
above z = 0.8 and for Dn4000 bins 2-4 fainter than MABS_G = -19.5 (v3lib.db_eff). Stacks of the 3000-3100 slices (z >= 0.30)
are multiplied by R(DB_EFF, M) from the lookup on the version's own axis: version G the MABS_G lookup at MABS_G_MED, version 5450
the M_AB(5450) lookup at M5450_MED. R_DIRECT (where the stack covers both windows) is kept as a check.
Output <stacks>_chain/ (as 19): stacks_z*.fits + cell_chain.csv.
Usage: python3 23_chain_v3.py --version G|5450 --stacks DIR --tables GLOB
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import fitsio
import v3lib

W5450, W3050 = (5400, 5500), (3000, 3100)
SWITCH_Z, DZ, COVER = 0.30, 0.05, 0.8
KEY = ['ZB', 'MB', 'CB', 'EB', 'DB']
QA = os.path.expanduser('~/data/sdssdesi_passC_v2_resp/qa_2win')
LOOKUP = {'G': (f'{QA}/gal_ratio_lookup_uniform2000_program_z_sdssfill.csv', 'MABS_G', 'MABS_G_MED'),
          '5450': (f'{QA}/gal_ratio_lookup_uniform2000_program_z_sdssfill_MSPEC.csv', 'M_AB5450', 'M5450_MED')}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['G', '5450']); ap.add_argument('--stacks', required=True)
    ap.add_argument('--tables', required=True); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or a.stacks.rstrip('/') + '_chain'; os.makedirs(out, exist_ok=True)
    cols = ['Z', 'ZB', 'CB', 'EB', 'DB', f'MB_{a.version}', f'ACC_{a.version}', 'NORM_OK', 'MABS_G_FULL', 'M_AB5450']
    import pyarrow.parquet as pq                       # tables of coadds without galaxies lack the v3 columns: skip them
    fs = [f for f in sorted(glob.glob(a.tables)) if set(cols) <= set(pq.read_schema(f).names)]
    t = pd.concat([pd.read_parquet(f, columns=cols) for f in fs], ignore_index=True)
    t = t[(t[f'ACC_{a.version}'] >= 0) & t.NORM_OK].rename(columns={f'MB_{a.version}': 'MB'})
    med = t.groupby(KEY).agg(MABS_G_MED=('MABS_G_FULL', 'median'), M5450_MED=('M_AB5450', 'median')).reset_index()
    lp, lx, mcol = LOOKUP[a.version]
    L = pd.read_csv(lp).sort_values(lx); look = {db: (g[lx].to_numpy(), g.logR.to_numpy()) for db, g in L.groupby('DB')}
    rows = []
    for fn in sorted(glob.glob(os.path.join(a.stacks, 'stacks_z*.fits'))):
        F = fitsio.FITS(fn); f = F['FLUX'].read(); iv = F['IVAR'].read(); sc = F['SCATTER'].read(); ns = F['NSPEC'].read()
        w = F['WAVE'].read(); c = F['CELLS'].read(); c = pd.DataFrame(c.astype(c.dtype.newbyteorder('='))); F.close()
        c = c.drop(columns=[x for x in ('MABS_G_MED', 'M5450_MED') if x in c]).merge(med, on=KEY, how='left')
        zb = int(c.ZB.iloc[0]); upper = (zb + 1) * DZ > SWITCH_Z + 1e-9
        s5 = (w >= W5450[0]) & (w <= W5450[1]); s3 = (w >= W3050[0]) & (w <= W3050[1]); own = s3 if upper else s5
        N = np.maximum(c.N_STACKED.to_numpy(), 1)
        with np.errstate(all='ignore'):
            m_own = np.nanmean(np.where(ns[:, own] > 0, f[:, own], np.nan), axis=1)
            m5 = np.nanmean(np.where(ns[:, s5] > 0, f[:, s5], np.nan), axis=1); m3 = np.nanmean(np.where(ns[:, s3] > 0, f[:, s3], np.nan), axis=1)
        cov5 = ns[:, s5].min(axis=1) / N; cov3 = ns[:, s3].min(axis=1) / N
        c['DB_EFF'] = v3lib.db_eff(c.Z_MED, c.DB, c.EW_MED, c.MABS_G_MED)
        rlook = np.array([10**np.interp(m, *look[int(d)]) if np.isfinite(m) else np.nan for m, d in zip(c[mcol], c.DB_EFF)])
        renorm = 1.0 / m_own; chain = rlook if upper else np.ones(len(c)); k = renorm * chain
        c['NORM_WIN'] = '3000-3100' if upper else '5400-5500'; c['RENORM'] = renorm; c['R_LOOKUP'] = rlook
        c['R_DIRECT'] = np.where((cov5 >= COVER) & (cov3 >= COVER) & (m5 > 0), m3 / m5, np.nan); c['CHAIN'] = chain
        tmp = os.path.join(out, os.path.basename(fn)) + '.tmp'
        with fitsio.FITS(tmp, 'rw', clobber=True) as fo:
            fo.write((f * k[:, None]).astype(np.float32), extname='FLUX'); fo.write((iv / k[:, None]**2).astype(np.float32), extname='IVAR')
            fo.write((sc * k[:, None]).astype(np.float32), extname='SCATTER'); fo.write(ns, extname='NSPEC'); fo.write(w, extname='WAVE')
            fo.write(c.to_records(index=False), extname='CELLS')
        os.replace(tmp, os.path.join(out, os.path.basename(fn))); rows.append(c.assign(FILE=os.path.basename(fn)))
    C = pd.concat(rows, ignore_index=True); C.to_csv(os.path.join(out, 'cell_chain.csv'), index=False)
    up = C[C.NORM_WIN == '3000-3100']; both = up[(up.R_DIRECT > 0) & np.isfinite(up.R_LOOKUP)]; d = np.log10(both.R_DIRECT / both.R_LOOKUP)
    print(f'{out}: version {a.version}, {len(C)} stacks ({int(C.N_STACKED.sum()):,} galaxies), {len(up)} chained from 3000-3100 '
          f'(no magnitude {up[mcol].isna().sum()}), re-typed {(C.DB_EFF != C.DB).sum()}; renorm median {np.nanmedian(C.RENORM):.4f}')
    if len(both):
        print(f'  check ({len(both)} stacks covering both windows): log R_direct - log R_lookup median {np.median(d):+.3f}, '
              f'16-84 % {np.percentile(d, 16):+.3f}..{np.percentile(d, 84):+.3f} dex')


if __name__ == '__main__':
    main()
