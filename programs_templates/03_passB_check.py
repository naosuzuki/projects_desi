#!/usr/bin/env python3
"""Checks on pass-B output: Step 1 star validation, timing, line sanity, rebin equivalence with astrolib.

Usage: python3 03_passB_check.py [outdir]
"""
import os, sys, glob
import numpy as np
import pandas as pd

OUT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/data/sdssdesi_passB')
files = sorted(glob.glob(os.path.join(OUT, '*', '*.parquet')))
df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
print(f'files {len(files)}  spectra {len(df):,}  mean time/file {df.groupby("FILE").PROC_SEC.first().mean():.1f}s')
fs = df['COADD_FIBERSTATUS'] if 'COADD_FIBERSTATUS' in df else df['COADD_FIBERSTATUS_x']
good = (df.ZWARN == 0) & (fs == 0) & df.ZCAT_PRIMARY.fillna(False).astype(bool)
print('with catalog redshift:', df.Z.notna().sum(), ' good:', good.sum())

def stats(x):
    x = x[np.isfinite(x)]
    if len(x) == 0: return 'n=0'
    q = np.percentile(x, [16, 50, 84])
    return f'n={len(x):>6,} median={q[1]:+.3f} sigma68={(q[2]-q[0])/2:.3f}'

print('\n== Step 1: catalog minus synthetic fiber magnitude BEFORE correction (spectrum x MEAN_PSF_TO_FIBER_SPECFLUX)')
for t in ['STAR', 'GALAXY', 'QSO']:
    s = df[good & (df.SPECTYPE == t) & (df.SNR_ALL > 3)]
    print(f'  {t:6s} g: {stats(s.CAL_DG.to_numpy())}   r: {stats(s.CAL_DR.to_numpy())}')
print('   (stars should be 0 +- 0.05; a common offset means the convention is wrong)')
print('\n== Step 1: z-band check AFTER the g,r correction (independent band)')
for t in ['STAR', 'GALAXY', 'QSO']:
    s = df[good & (df.SPECTYPE == t) & (df.SNR_ALL > 3) & df.CAL_OK]
    print(f'  {t:6s} z: {stats(s.CAL_DZ_CHECK.to_numpy())}   CAL_OK fraction {df[good & (df.SPECTYPE==t)].CAL_OK.mean():.3f}')
print('\n== Step 1: correction parameters (a = scale at 6000 A, b = tilt per dex)')
s = df[good & df.CAL_OK]
print('  all   a:', stats(s.CAL_A.to_numpy()), '  b:', stats(s.CAL_B.to_numpy()))
if 'CAL_MODE' in df:
    for t in ['STAR', 'GALAXY', 'QSO']:
        m = df[good & (df.SPECTYPE == t)].CAL_MODE.value_counts().to_dict()
        s2 = df[good & (df.SPECTYPE == t) & (df.CAL_MODE == 2)]
        print(f'  {t:6s} modes {m}   tilt (mode 2): {stats(s2.CAL_B.to_numpy())}   syn err g: {stats(df[good & (df.SPECTYPE==t)].SYN_FIBER_G_ERR.to_numpy())}')

print('\n== Step 5: S/N by type (median per-pixel S/N, good spectra)')
for t in ['STAR', 'GALAXY', 'QSO']:
    s = df[good & (df.SPECTYPE == t)]
    print(f'  {t:6s} n={len(s):>6,} median SNR_ALL={s.SNR_ALL.median():.2f} SDSS-pixel scale={s.SNR_SDSSPIX.median():.2f}')

print('\n== Step 5: lines in galaxies (fraction detected at >3 sigma, median EW of detections)')
g = df[good & (df.SPECTYPE == 'GALAXY')]
for name in ['OII_SUM', 'OII_3726', 'OII_3729', 'HBETA', 'OIII_5007', 'HALPHA', 'NII_6583']:
    if f'{name}_FLUX' not in g: continue
    fl, er = g[f'{name}_FLUX'], g[f'{name}_FLUX_ERR']
    cov = fl.notna(); det = cov & (fl / er > 3)
    print(f'  {name:10s} covered {cov.sum():>6,}  detected {det.sum():>6,} ({det.sum()/max(cov.sum(),1)*100:4.1f}%)  median EW(det) {g.loc[det, name+"_EW"].median():6.1f} A')
print('  Dn4000 galaxies:', stats(g.DN4000.to_numpy()), '  LRG-like (Dn4000>1.6):', (g.DN4000 > 1.6).mean().round(3))
print('  Halpha sigma choice (km/s):', g.ha_SIGMA.value_counts(dropna=True).to_dict())

print('\n== Step 4: absolute magnitudes (galaxies)')
for b, (lo, hi) in {'R': (0, 0.35), 'G': (0.35, 0.72), 'U': (0.72, 1.45), 'NUV': (1.45, 1.6)}.items():
    s = g[(g.Z >= lo) & (g.Z < hi)]
    print(f'  M_{b:3s} z {lo}-{hi}: n={len(s):>6,} covered {s[f"MABS_{b}"].notna().mean()*100:5.1f}%  median {s[f"MABS_{b}"].median():.2f}  K median {s[f"K_{b}"].median():+.2f}')

# ---- rebin equivalence with sdsslib astrolib on one spectrum
try:
    sys.path.append(os.environ.get('PYLIB', os.path.expanduser('~/github/sdsslib/pylib')))
    import numpy, math
    numpy.int = int   # astrolib uses the removed alias
    import astrolib, fitsio
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import importlib; pb = importlib.import_module('02_passB_measure')
    p = [l.split()[1] for l in open(pb.MANIFEST)][0]
    p = pb.local_path(p)
    f = fitsio.FITS(p); w = f['B_WAVELENGTH'].read(); fl = f['B_FLUX'][0:1, :][0].astype(float); iv = f['B_IVAR'][0:1, :][0].astype(float); mk = f['B_MASK'][0:1, :][0]
    rw, rf, ri, rid, rm = astrolib.exec_logrebinning_fluxivarmask(1e-4, w, fl, iv, mk)
    rid2, M, M2 = pb.rebin_matrix(w)
    rf2 = (M @ fl)
    common = np.intersect1d(rid, rid2); i1 = np.searchsorted(rid, common); i2 = np.searchsorted(rid2, common)
    d = rf2[i2] - rf[i1]
    print(f'\n== rebin vs astrolib (B camera, 1 spectrum): common pixels {len(common)}, max |diff| {np.nanmax(np.abs(d)):.3e}, rms {np.sqrt(np.nanmean(d**2)):.3e}, flux rms {np.sqrt(np.nanmean(rf[i1]**2)):.3e}')
except Exception as e:
    print('\nrebin comparison skipped:', repr(e))
