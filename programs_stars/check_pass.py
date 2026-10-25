#!/usr/bin/env python3
"""Checks on pass-S output: Gaia-G/BP/RP check of the pipeline flux scale, Step-1 fiber-flux check for stars,
RV distribution, spectral indices and band ratios per Redrock SUBTYPE.
Usage: python3 check_pass.py [outdir]   (default ~/data/sdssdesi_passS)
"""
import os, sys, glob
import numpy as np, pandas as pd

OUT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/data/sdssdesi_passS')
files = sorted(glob.glob(os.path.join(OUT, '*', '*.parquet')))
df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
print(f'files {len(files)}  stars {len(df):,}  mean time/file {df.groupby("FILE").PROC_SEC.first().mean():.1f}s  columns {df.shape[1]}')
good = (df.ZWARN == 0) & ((df.COADD_FIBERSTATUS == 0) if 'COADD_FIBERSTATUS' in df else True)
if 'SUBTYPE' not in df and 'SUBCLASS1' in df: df['SUBTYPE'] = df.SUBCLASS1
print('good (ZWARN=0, fiberstatus=0):', good.sum(), ' SUBTYPE:', df[good].SUBTYPE.value_counts().to_dict())

def stats(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0: return 'n=0'
    q = np.percentile(x, [16, 50, 84])
    return f'n={len(x):>5,} median={q[1]:+.3f} sigma68={(q[2]-q[0])/2:.3f}'

print('\n== Gaia check: catalog (Vega) minus synthetic from the pipeline spectrum (point-source flux, no fiber factor)')
for snr_lo, snr_hi in ((3, 10), (10, 30), (30, 1e9)):
    s = df[good & (df.SNR_ALL >= snr_lo) & (df.SNR_ALL < snr_hi)]
    print(f'  S/N {snr_lo:>3}-{snr_hi if snr_hi < 1e9 else "":<4}  G: {stats(s.GAIA_DG)}   BP: {stats(s.GAIA_DBP)}   RP: {stats(s.GAIA_DRP)}   BP-RP: {stats(s.GAIA_DBPRP)}')
print('  band weight in the padding (extrapolated): G %.3f  BP %.3f  RP %.3f' % tuple(df[f'GAIA_EXTRAP_{b}'].iloc[0] for b in ('G', 'BP', 'RP')))
print('  by SUBTYPE (S/N>10), G and BP-RP:')
for t, s in df[good & (df.SNR_ALL > 10)].groupby('SUBTYPE'):
    print(f'    {t:3s} G: {stats(s.GAIA_DG)}   BP-RP: {stats(s.GAIA_DBPRP)}')

s = df[good & (df.SNR_ALL > 3)]
if 'CAL_DZ_CHECK' in df:
    print('\n== Step 2 (DESI, pass B): catalog minus synthetic FIBER mag, spectrum x MEAN_PSF_TO_FIBER_SPECFLUX (expect 0 +- 0.05)')
    print(f'  g: {stats(s.CAL_DG)}   r: {stats(s.CAL_DR)}   z after scale: {stats(s[s.CAL_OK].CAL_DZ_CHECK)}   CAL_OK {s.CAL_OK.mean():.3f}')
else:
    print('\n== Step 2 (SDSS): PSF mag minus synthetic (expect 0 +- 0.05); fiber mag minus synthetic for the record')
    for sv, g2 in s.groupby('SURVEY'):
        print(f'  {sv}: g: {stats(g2.CAL_DG)}   r: {stats(g2.CAL_DR)}   CAL_OK {g2.CAL_OK.mean():.3f}   fiber g: {stats(g2.CAL_DG_FIBER)}   Gaia G after the PSF scale: {stats(g2.GAIA_DG_CAL)}')
m = np.isfinite(s.GAIA_DG.to_numpy(float)) & np.isfinite(s.CAL_DR.to_numpy(float))
print(f'  Gaia G vs catalog r residuals: corr(GAIA_DG, CAL_DR) = {np.corrcoef(s.GAIA_DG[m], s.CAL_DR[m])[0,1]:.2f}  (n={m.sum()})')

print('\n== RV (km/s) from Redrock z and the integer-shift residual')
if 'DIST_PC' in df: print('  distance pc 5/50/95:', np.nanpercentile(df[good].DIST_PC, [5,50,95]).round(0), ' M_G_SYN0 - M_G0 (abs. synthetic G vs Gaia):', stats(df[good].M_G_SYN0 - (df[good].GAIA_PHOT_G_MEAN_MAG - df[good].DM - df[good].A_G_USED)))
print('  RV:', stats(df[good].RV_KMS), '  |RV|>300:', f'{(np.abs(df[good].RV_KMS) > 300).mean():.3f}', '  residual:', stats(df[good].RV_RESID_KMS))

print('\n== indices by SUBTYPE (median, good, S/N>10): EW in A (positive = absorption), ratios dimensionless')
cols = ['CAK_EW', 'HDELTA_A_EW', 'G4300_EW', 'HBETA_LK_EW', 'MGB_EW', 'FE5270_EW', 'NAD_EW', 'HALPHA_ST_EW', 'CAT_EW', 'HEI_5876_EW', 'HEII_4686_EW', 'TIO5', 'CAH2', 'CAH3', 'C_BLUE_DERED', 'C_RED_DERED']
s = df[good & (df.SNR_ALL > 10)]
print('  ' + ' '.join(f'{c[:9]:>9s}' for c in ['SUBTYPE', 'n'] + cols))
for t, g in s.groupby('SUBTYPE'):
    print('  ' + f'{t:>9s} {len(g):>9d} ' + ' '.join(f'{g[c].median():9.2f}' for c in cols))
print('  coverage (finite fraction, all good):', {c: round(float(df[good][c].notna().mean()), 2) for c in ['CAK_EW', 'HBETA_LK_EW', 'HALPHA_ST_EW', 'CAT_EW', 'TIO5']})
print('  Halpha emission (EW < -1 A, >3 sigma):', int(((df.HALPHA_ST_EW < -1) & (df.HALPHA_ST_EW / df.HALPHA_ST_EW_ERR < -3) & good).sum()), 'of', int(good.sum()))
