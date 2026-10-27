#!/usr/bin/env python3
"""Scan all pass-B tables once (internal disk) and write exact per-class S/N counts, star calibration statistics
and a compact per-spectrum column subset for the cell assignment.

Outputs: ~/data/sdssdesi_products/passB_summary_columns.parquet  (54 M rows x ~25 columns, ~5 GB)
         ~/data/sdssdesi_products/passB_summary.json
"""
import os, glob, json, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

OUT = os.path.expanduser('~/data/sdssdesi_products')
COLS = ['TARGETID', 'PROGRAM', 'SPECTYPE', 'Z', 'ZWARN', 'ZCAT_PRIMARY', 'DELTACHI2', 'SNR_ALL', 'SNR_SDSSPIX',
        'CAL_DG', 'CAL_DR', 'CAL_DZ_CHECK', 'CAL_OK', 'CAL_MODE', 'CAL_A', 'CAL_B_EST', 'DN4000', 'DN4000_ERR',
        'OII_SUM_EW', 'OII_SUM_EW_ERR', 'OII_SUM_FLUX', 'OII_SUM_FLUX_ERR', 'HALPHA_EW', 'HALPHA_EW_ERR',
        'HALPHA_FLUX', 'HALPHA_FLUX_ERR', 'HBETA_FLUX', 'HBETA_FLUX_ERR', 'MABS_R', 'MABS_G', 'MABS_U', 'MABS_NUV',
        'K_R', 'K_G', 'K_U', 'K_NUV', 'M_R_TOTAL_DERED', 'FIBER_FRACTION_R', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'EBV', 'PHOTSYS']

t0 = time.time()
files = sorted(glob.glob(os.path.expanduser('~/data/sdssdesi_passB/*/*.parquet')))
parts = []
for i, f in enumerate(files):
    names = pq.read_schema(f).names
    want = [c for c in COLS if c in names] + [c for c in ('COADD_FIBERSTATUS', 'COADD_FIBERSTATUS_x') if c in names]
    d = pq.read_table(f, columns=want).to_pandas()
    if 'COADD_FIBERSTATUS_x' in d:
        d = d.rename(columns={'COADD_FIBERSTATUS_x': 'COADD_FIBERSTATUS'})
    d['HEALPIX'] = int(os.path.basename(f)[3:-8])
    parts.append(d)
    if i % 2000 == 0:
        print(i, f'{time.time()-t0:.0f}s', flush=True)
d = pd.concat(parts, ignore_index=True); del parts
print('rows', len(d), f'{time.time()-t0:.0f}s', flush=True)
d.to_parquet(f'{OUT}/passB_summary_columns.parquet')
print('written columns parquet', f'{time.time()-t0:.0f}s', flush=True)

good = (d.ZWARN == 0) & (d.COADD_FIBERSTATUS == 0) & d.ZCAT_PRIMARY.fillna(False).astype(bool)
out = {'rows': int(len(d)), 'good': int(good.sum()), 'files': len(files)}
for t in ['GALAXY', 'QSO', 'STAR']:
    g = d[good & (d.SPECTYPE == t)]
    out[t] = {'n': int(len(g)), 'median_snr_sdsspix': float(g.SNR_SDSSPIX.median()),
              **{f'snr_gt_{th}': int((g.SNR_SDSSPIX > th).sum()) for th in (5, 10, 20)},
              'by_program_snr_gt_20': {p: int(((g.PROGRAM == p) & (g.SNR_SDSSPIX > 20)).sum()) for p in ('bright', 'dark')}}
s = d[good & (d.SPECTYPE == 'STAR') & (d.SNR_ALL > 3)]
q = lambda x: [float(v) for v in np.nanpercentile(x, [16, 50, 84])]
out['star_cal'] = {'n': int(len(s)), 'dg_p16_50_84': q(s.CAL_DG), 'dr_p16_50_84': q(s.CAL_DR),
                   'dz_after_p16_50_84': q(s[s.CAL_OK].CAL_DZ_CHECK), 'cal_ok_frac': float(s.CAL_OK.mean())}
gal = d[good & (d.SPECTYPE == 'GALAXY') & (d.Z > 0) & (d.Z < 1.6)]
with np.errstate(all='ignore'):
    out['galaxy_lines'] = {'n': int(len(gal)),
                           'halpha_det_frac_covered': float(((gal.HALPHA_FLUX / gal.HALPHA_FLUX_ERR) > 3).sum() / max(gal.HALPHA_FLUX.notna().sum(), 1)),
                           'oii_det_frac_covered': float(((gal.OII_SUM_FLUX / gal.OII_SUM_FLUX_ERR) > 3).sum() / max(gal.OII_SUM_FLUX.notna().sum(), 1)),
                           'dn4000_median': float(gal.DN4000.median()), 'dn4000_gt_1.6_frac': float((gal.DN4000 > 1.6).mean())}
json.dump(out, open(f'{OUT}/passB_summary.json', 'w'), indent=1)
print(json.dumps(out, indent=1)); print('done', f'{time.time()-t0:.0f}s')
