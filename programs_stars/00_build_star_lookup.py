#!/usr/bin/env python3
"""Step 0a of the stellar thread: all main-survey STAR rows of zall-pix-loa.fits with the stellar and Gaia columns
that the pass-B compact lookup does not carry (SUBTYPE, REF_ID = Gaia DR3 source_id, Gaia G/BP/RP, PARALLAX, PM,
MWS_TARGET, ...) -> ~/data/sdssdesi_cache/zall_stars_main.parquet  (12.18 M rows, 2 min, chunked: < 2 GB RAM).
Reads the internal disk only (safe while pass B runs). Prints a census at the end.
Usage: python3 00_build_star_lookup.py
"""
import os, time
import numpy as np, pandas as pd, fitsio

ZALL = os.path.expanduser('~/desiredux/zall-pix-loa.fits')
OUT = os.path.expanduser('~/data/sdssdesi_cache/zall_stars_main.parquet')
COLS = ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'DELTACHI2', 'COADD_FIBERSTATUS',
        'ZCAT_PRIMARY', 'EBV', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
        'REF_ID', 'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'PARALLAX', 'PMRA', 'PMDEC',
        'MWS_TARGET', 'SCND_TARGET', 'DESI_TARGET', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_EXPTIME', 'COADD_NUMEXP', 'MEAN_MJD',
        'TSNR2_BGS', 'TSNR2_LRG', 'TSNR2_QSO', 'TARGET_RA', 'TARGET_DEC']
CKMS = 299792.458


def build(step=2_000_000):
    t0 = time.time()
    f = fitsio.FITS(ZALL)['ZCATALOG']
    n = f.get_nrows(); parts = []
    for s in range(0, n, step):
        r = f.read(rows=np.arange(s, min(s + step, n)), columns=COLS)
        r = r[(r['SPECTYPE'] == 'STAR') & (r['SURVEY'] == 'main')]
        parts.append(pd.DataFrame({c: (np.char.strip(r[c].astype(str)) if r[c].dtype.kind in 'SU'
                                       else r[c].astype(r[c].dtype.newbyteorder('='))) for c in COLS}))
        print(f'{min(s + step, n):>10,}/{n:,} kept {sum(len(p) for p in parts):,} {time.time() - t0:.0f}s', flush=True)
    df = pd.concat(parts, ignore_index=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_parquet(OUT + '.tmp'); os.replace(OUT + '.tmp', OUT)
    print('written', OUT, f'{len(df):,} rows, {time.time() - t0:.0f}s')
    return df


def census(d):
    print('main STAR rows', f'{len(d):,}', ' primary', f'{d.ZCAT_PRIMARY.sum():,}')
    g = d[d.ZCAT_PRIMARY & (d.ZWARN == 0) & (d.COADD_FIBERSTATUS == 0)]
    print('good (primary, ZWARN=0, fiberstatus=0):', f'{len(g):,}', ' by program:', g.PROGRAM.value_counts().to_dict())
    gb = g[g.PROGRAM.isin(['dark', 'bright'])]
    print('dark+bright good:', f'{len(gb):,}', ' by SUBTYPE:', gb.SUBTYPE.value_counts().to_dict())
    gaia = gb.REF_CAT == 'G2'
    print('with Gaia (REF_CAT G2):', f'{gaia.sum():,}', ' G<17.65 (XP):', f'{(gaia & (gb.GAIA_PHOT_G_MEAN_MAG < 17.65)).sum():,}',
          ' G percentiles 5/50/95:', np.nanpercentile(gb.GAIA_PHOT_G_MEAN_MAG[gaia], [5, 50, 95]).round(2))
    rv = gb.Z.to_numpy() * CKMS
    print('RV = cz km/s percentiles 1/50/99:', np.percentile(rv, [1, 50, 99]).round(1), ' |RV|>300 km/s fraction:',
          round(float((np.abs(rv) > 300).mean()), 4), ' median c*ZERR km/s:', round(float(gb.ZERR.median() * CKMS), 2))
    mws = gb.MWS_TARGET.astype(np.int64) != 0
    print('MWS_TARGET != 0:', f'{mws.sum():,}', ' secondary only:', f'{((~mws) & (gb.SCND_TARGET.astype(np.int64) != 0)).sum():,}',
          ' neither:', f'{((~mws) & (gb.SCND_TARGET.astype(np.int64) == 0)).sum():,}')
    print('backup-program good stars (coadds NOT on disk):', f'{(g.PROGRAM == "backup").sum():,}')


if __name__ == '__main__':
    df = pd.read_parquet(OUT) if os.path.exists(OUT) else build()
    census(df)
