#!/usr/bin/env python3
"""Stars, step 3b (Nao, 2026-10-24 night): Milky Way E(B-V) in front of each star from 3D dust maps at the parallax distance.

Maps (dustmaps, data in ~/data/GAIA; the global dustmaps config is not changed):
  Edenhofer+24 (dustmaps 'edenhofer2023', mean of the integrated extinction; 69 pc - 1.25 kpc, 14') for d < 1.25 kpc,
  Bayestar19 (Green+19; Dec > -30 deg, 3.4'-13.7') beyond, and as the fallback where Edenhofer has no value.
Units (calibrated empirically on the first run, see the constants): Bayestar19's E is on the raw SFD scale -> E(B-V)_SF11 =
  0.86 E_B19 (matches 0.86 x SFD at d > 3 kpc to 0.1 %); Edenhofer's E (Zhang, Green & Rix 2023 unit) is already on the SF11
  E(B-V) scale -> factor 1.00 (matches Bayestar in the 0.3-1.25 kpc overlap to 0.6 %).
Stars closer than 69 pc (inside Edenhofer's inner edge) get the 69 pc value scaled by d/69 pc (local dust ~ uniform; SRC = 3).
Columns: EBV_B19, EBV_E24 (converted to E(B-V)), EBV_3D (the combination), EBV_3D_SRC (2 Edenhofer, 4 Bayestar, 3 inner scaling,
  0 none), EBV_3D_ERR (NaN with the mean-only Edenhofer file).
Inputs/outputs:
  SDSS: ~/data/sdssdesi_cache/sdss_stars_gaia.parquet (RA, DEC, PLX_CORR) -> same file, columns added
  DESI: ~/data/sdssdesi_cache/zall_stars_main.parquet (TARGET_RA/DEC, PARALLAX) -> ~/data/sdssdesi_cache/desi_stars_ebv3d.parquet
Usage: python3 step3b_ebv3d.py [--sdss] [--desi]   (default: both)
"""
import os, sys, time, argparse
import numpy as np, pandas as pd
from dustmaps.config import config
config['data_dir'] = os.path.expanduser('~/data/GAIA')
import astropy.units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from dustmaps.edenhofer2023 import Edenhofer2023Query

CACHE = os.path.expanduser('~/data/sdssdesi_cache')
# Unit calibration (first run, 2026-10-24 night, 10.3 M DESI stars): with B19_TO_EBV = 1 the Bayestar19 value at d > 3 kpc was
# 1.161 x (0.86 SFD), i.e. Bayestar19's unit is the raw SFD E(B-V) scale, so the SF11 factor applies to it too; and with
# E24_TO_EBV = 2.8/3.1 the Edenhofer/Bayestar ratio at 0.3-1.25 kpc was 0.772 = 0.903 x 0.855, i.e. Edenhofer's E is already
# on the SF11 E(B-V) scale. Both ratios become 1.00 with the factors below.
B19_TO_EBV = 0.86
E24_TO_EBV = 1.00
E24_DMIN, E24_DMAX = 69.0, 1250.0
CHUNK = 500_000


STATUS = os.path.expanduser('~/data/sdssdesi_status/step_star_ebv3d.json')
_T0 = time.time()


def _progress(done, total, status='running'):
    """progress protocol of the monitors (monitor_stars.py reads this file)."""
    try:
        os.makedirs(os.path.dirname(STATUS), exist_ok=True)
        tmp = STATUS + '.tmp'
        import json; json.dump({'done': int(done), 'total': int(total), 'status': status, 'started': _T0}, open(tmp, 'w')); os.replace(tmp, STATUS)
    except Exception:
        pass


class Dust3D:
    def __init__(self):
        t0 = time.time()
        self.b19 = BayestarQuery(version='bayestar2019', max_samples=1)
        self.e24 = Edenhofer2023Query(integrated=True, load_samples=False)
        print(f'maps loaded {time.time()-t0:.0f}s', flush=True)

    def query(self, ra, dec, d_pc):
        n = len(ra); out = {k: np.full(n, np.nan) for k in ('EBV_B19', 'EBV_E24', 'EBV_3D', 'EBV_3D_ERR')}
        src = np.zeros(n, np.int8)
        ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(d_pc) & (d_pc > 0)
        for i0 in range(0, n, CHUNK):
            sl = slice(i0, min(i0 + CHUNK, n)); m = ok[sl]
            if not m.any():
                continue
            idx = np.where(m)[0] + i0
            d = d_pc[idx]
            c_b = SkyCoord(ra[idx] * u.deg, dec[idx] * u.deg, distance=np.clip(d, 1.0, 6e4) * u.pc, frame='icrs')
            out['EBV_B19'][idx] = B19_TO_EBV * np.asarray(self.b19(c_b, mode='median'), float)
            d_e = np.clip(d, E24_DMIN, E24_DMAX)
            c_e = SkyCoord(ra[idx] * u.deg, dec[idx] * u.deg, distance=d_e * u.pc, frame='icrs')
            e = E24_TO_EBV * np.asarray(self.e24(c_e, mode='mean'), float)
            es = np.full(len(e), np.nan)          # the mean-only map carries no samples: no per-star error (needs samples_healpix.fits, 19.5 GB)
            inner = d < E24_DMIN
            e = np.where(inner, e * d / E24_DMIN, e); es = np.where(inner, es * d / E24_DMIN, es)
            e = np.where(d <= E24_DMAX, e, np.nan); es = np.where(d <= E24_DMAX, es, np.nan)
            out['EBV_E24'][idx] = e
            use_e = np.isfinite(e)
            out['EBV_3D'][idx] = np.where(use_e, e, out['EBV_B19'][idx]); out['EBV_3D_ERR'][idx] = np.where(use_e, es, np.nan)
            src[idx] = np.where(use_e, np.where(inner, 3, 2), np.where(np.isfinite(out['EBV_B19'][idx]), 4, 0))
            print(f'  {min(i0 + CHUNK, n):>10,}/{n:,}', flush=True)
            _progress(min(i0 + CHUNK, n), n)
        out['EBV_3D_SRC'] = src
        return out


def report(df, ebv_sfd, d):
    e3 = df['EBV_3D']; src = df['EBV_3D_SRC']
    print('  source counts:', pd.Series(src).value_counts().to_dict(), ' (2 Edenhofer, 3 inner scaling, 4 Bayestar, 0 none)')
    far = (d > 3000) & np.isfinite(df['EBV_B19']) & (ebv_sfd > 0.02)
    if far.sum() > 100:
        print(f'  d > 3 kpc, n={far.sum():,}: median EBV_B19 / (0.86 SFD) = {np.nanmedian(df["EBV_B19"][far] / (0.86 * ebv_sfd[far])):.3f}  (expect ~1 if the Bayestar unit is E(B-V)_SF11)')
    ov = (d > 300) & (d < 1250) & np.isfinite(df['EBV_E24']) & np.isfinite(df['EBV_B19']) & (df['EBV_B19'] > 0.02)
    if ov.sum() > 100:
        print(f'  0.3-1.25 kpc, n={ov.sum():,}: median EBV_E24 / EBV_B19 = {np.nanmedian(df["EBV_E24"][ov] / df["EBV_B19"][ov]):.3f}  (expect ~1 after the unit conversions)')
    near = (d < 300) & np.isfinite(e3)
    if near.sum() > 100:
        print(f'  d < 300 pc, n={near.sum():,}: median EBV_3D / (0.86 SFD) = {np.nanmedian(e3[near] / np.where(ebv_sfd[near] > 0, 0.86 * ebv_sfd[near], np.nan)):.3f}  ->  the fraction of the column in front of the nearby stars')
    for lo, hi in ((0, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 1e9)):
        s = (d >= lo) & (d < hi) & np.isfinite(e3)
        if s.sum():
            print(f'  {lo:>5.0f}-{hi if hi < 1e9 else 99999:>5.0f} pc n={s.sum():>9,}  median E(B-V)_3D {np.nanmedian(e3[s]):.4f}  0.86 SFD {np.nanmedian(0.86 * ebv_sfd[s]):.4f}')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sdss', action='store_true'); ap.add_argument('--desi', action='store_true')
    a = ap.parse_args(); do_sdss = a.sdss or not (a.sdss or a.desi); do_desi = a.desi or not (a.sdss or a.desi)
    dust = Dust3D()
    if do_sdss:
        p = f'{CACHE}/sdss_stars_gaia.parquet'; df = pd.read_parquet(p)
        d = (1000.0 / df.PLX_CORR.where(df.PLX_CORR > 0)).to_numpy(float)
        print(f'SDSS: {len(df):,} rows, {np.isfinite(d).sum():,} with a positive parallax', flush=True)
        out = dust.query(df.RA.to_numpy(float), df.DEC.to_numpy(float), d)
        for k, v in out.items():
            df[k] = v
        df.to_parquet(p + '.tmp'); os.replace(p + '.tmp', p)
        # SFD for the SDSS stars is not in the catalog (it comes from the plugmap in the pass); compare with the DESI-like far limit later
        print('  source counts:', pd.Series(out['EBV_3D_SRC']).value_counts().to_dict())
        for lo, hi in ((0, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 1e9)):
            s = (d >= lo) & (d < hi) & np.isfinite(out['EBV_3D'])
            if s.sum(): print(f'  {lo:>5.0f}-{hi if hi < 1e9 else 99999:>5.0f} pc n={s.sum():>8,}  median E(B-V)_3D {np.nanmedian(out["EBV_3D"][s]):.4f}')
    if do_desi:
        z = pd.read_parquet(f'{CACHE}/zall_stars_main.parquet', columns=['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'PARALLAX', 'EBV', 'ZCAT_PRIMARY', 'REF_CAT'])
        z = z[z.ZCAT_PRIMARY & (z.REF_CAT == 'G2') & (z.PARALLAX > 0)].drop_duplicates('TARGETID').reset_index(drop=True)
        d = 1000.0 / (z.PARALLAX.to_numpy(float) + 0.017)
        print(f'DESI: {len(z):,} primary stars with a Gaia parallax', flush=True)
        out = dust.query(z.TARGET_RA.to_numpy(float), z.TARGET_DEC.to_numpy(float), d)
        res = pd.DataFrame({'TARGETID': z.TARGETID, 'DIST_PC': d, 'EBV_SFD': z.EBV.to_numpy(float), **out})
        res.to_parquet(f'{CACHE}/desi_stars_ebv3d.parquet')
        print(f'  written {CACHE}/desi_stars_ebv3d.parquet')
        report(res, res.EBV_SFD.to_numpy(float), d)
    _progress(1, 1, 'done')


if __name__ == '__main__':
    main()
