#!/usr/bin/env python3
"""Stars, step 1 (Nao, 2026-10-24 night): generate the catalog of SDSS stars that have a Gaia DR3 parallax.

Sources
  spAll-v5_13_2.fits (DR17 = BOSS/eBOSS/SEGUE on BOSS) and specObj-dr8.fits (SDSS-I/II + SEGUE): CLASS == 'STAR' and
  SPECPRIMARY == 1 rows with SUBCLASS, Z (= radial velocity / c), ZWARNING, S/N, PSF/fiber magnitudes (DR17), ELODIE
  Teff/logg/[Fe/H]/RV (the pipeline's stellar fit; truth columns for the SDSS stars).
  Gaia DR3 x SDSS star crossmatches of projects_gaia (csvfiles/gaiadr3_sdssdr{17,8}_star.csv): joined on (plate, mjd,
  fiber); when several Gaia sources match one fiber the nearest (smallest separation) is kept. The CSVs also carry the
  SDSS psfmag_u..z, used for DR8 where specObj has no PSF magnitudes.
Selection: a Gaia match with a finite parallax (Nao: "we collect the ones which have a Gaia parallax"). The other stars
are kept in the file with GAIA_PLX = False so the pass can include them with --all-stars.
Derived: PLX_CORR = parallax + 0.017 mas (global Lindegren+21 offset; the full recipe later), QUALITY = plx S/N > 10 &
RUWE < 1.4 & single & not variable (ms.tex), M_G, M_G0 and (BP-RP)_0 with GSP-Phot A_G / E(BP-RP) where available.
Output: ~/data/sdssdesi_cache/sdss_stars_gaia.parquet (one row per primary star spectrum).
Usage: python3 step1_sdss_catalog.py
"""
import os, time
import numpy as np, pandas as pd, fitsio

SPECTRO = os.path.expanduser('~/sdssredux')
SPALL = f'{SPECTRO}/spAll-v5_13_2.fits'; SPECOBJ = f'{SPECTRO}/specObj-dr8.fits'
GAIA_CSV = {'dr17': os.path.expanduser('~/github/projects_gaia/csvfiles/gaiadr3_sdssdr17_star.csv'),
            'dr8': os.path.expanduser('~/github/projects_gaia/csvfiles/gaiadr3_sdssdr8_star.csv')}
OUT = os.path.expanduser('~/data/sdssdesi_cache/sdss_stars_gaia.parquet')
GAIA_COLS = ['source_id', 'separation', 'parallax', 'parallax_error', 'parallax_over_error', 'pmra', 'pmdec', 'ruwe',
             'astrometric_params_solved', 'nu_eff_used_in_astrometry', 'pseudocolour', 'ecl_lat',
             'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'phot_bp_rp_excess_factor',
             'radial_velocity', 'radial_velocity_error', 'phot_variable_flag', 'non_single_star', 'has_xp_continuous', 'has_rvs',
             'teff_gspphot', 'logg_gspphot', 'mh_gspphot', 'ag_gspphot', 'ebpminrp_gspphot',
             'plate', 'mjd', 'fiber', 'psfmag_u', 'psfmag_g', 'psfmag_r', 'psfmag_i', 'psfmag_z']
ELODIE = ['ELODIE_TEFF', 'ELODIE_LOGG', 'ELODIE_FEH', 'ELODIE_Z', 'ELODIE_Z_ERR', 'ELODIE_BV', 'ELODIE_SPTYPE']
CKMS = 299792.458


def sdss_stars(survey):
    if survey == 'dr17':
        cols = ['PLATE', 'MJD', 'FIBERID', 'CLASS', 'SUBCLASS', 'Z', 'Z_ERR', 'ZWARNING', 'SPECPRIMARY', 'PLATEQUALITY', 'PROGRAMNAME',
                'SN_MEDIAN_ALL', 'PSFMAG', 'FIBER2MAG', 'FIBERMAG', 'EXTINCTION', 'THING_ID', 'RA', 'DEC'] + ELODIE
        t = fitsio.read(SPALL, columns=cols)
    else:
        cols = ['PLATE', 'MJD', 'FIBERID', 'CLASS', 'SUBCLASS', 'Z', 'Z_ERR', 'ZWARNING', 'SPECPRIMARY', 'PLATEQUALITY', 'PROGRAMNAME',
                'SN_MEDIAN', 'PLUG_RA', 'PLUG_DEC'] + ELODIE
        t = fitsio.read(SPECOBJ, columns=cols)
    cls = np.char.strip(t['CLASS'].astype(str))
    t = t[(cls == 'STAR') & (t['SPECPRIMARY'] == 1)]
    df = pd.DataFrame({'SURVEY': survey, 'PLATE': t['PLATE'].astype(np.int32), 'MJD': t['MJD'].astype(np.int32), 'FIBERID': t['FIBERID'].astype(np.int16),
                       'SUBCLASS': np.char.strip(t['SUBCLASS'].astype(str)), 'Z': t['Z'].astype(float), 'ZERR': t['Z_ERR'].astype(float),
                       'ZWARN': t['ZWARNING'].astype(np.int32), 'PLATEQUALITY': np.char.strip(t['PLATEQUALITY'].astype(str)),
                       'PROGRAMNAME': np.char.strip(t['PROGRAMNAME'].astype(str))})
    for c in ELODIE:
        df[c] = np.char.strip(t[c].astype(str)) if t[c].dtype.kind in 'SU' else t[c].astype(float)
    if survey == 'dr17':
        df['SNR_CAT'] = t['SN_MEDIAN_ALL'].astype(float)
        for i, b in enumerate('ugriz'):
            df[f'PSFMAG_{b.upper()}'] = t['PSFMAG'][:, i].astype(float); df[f'FIBER2MAG_{b.upper()}'] = t['FIBER2MAG'][:, i].astype(float)
            df[f'EXT_{b.upper()}'] = t['EXTINCTION'][:, i].astype(float)
        df['THING_ID'] = t['THING_ID'].astype(np.int64); df['RA'] = t['RA'].astype(float); df['DEC'] = t['DEC'].astype(float)
    else:
        sn = t['SN_MEDIAN'].astype(float); df['SNR_CAT'] = np.nanmedian(np.where(sn > 0, sn, np.nan)[:, 1:4], axis=1) if sn.ndim == 2 else sn
        df['THING_ID'] = -1; df['RA'] = t['PLUG_RA'].astype(float); df['DEC'] = t['PLUG_DEC'].astype(float)
    df['SUBCLASS1'] = df.SUBCLASS.str[:1].where(df.SUBCLASS.str.len() > 0, 'N')   # O B A F G K M W(D) C(V) L
    return df


def gaia(survey):
    g = pd.read_csv(GAIA_CSV[survey], usecols=GAIA_COLS)
    g = g.sort_values('separation').drop_duplicates(['plate', 'mjd', 'fiber'])
    g = g.rename(columns={'plate': 'PLATE', 'mjd': 'MJD', 'fiber': 'FIBERID'})
    g = g.rename(columns={c: 'GAIA_' + c.upper() for c in g.columns if c not in ('PLATE', 'MJD', 'FIBERID') and not c.startswith('psfmag')})
    g = g.rename(columns={f'psfmag_{b}': f'PSFMAG_{b.upper()}_CSV' for b in 'ugriz'})
    g['PLATE'] = g.PLATE.astype(np.int32); g['MJD'] = g.MJD.astype(np.int32); g['FIBERID'] = g.FIBERID.astype(np.int16)
    return g


def main():
    t0 = time.time(); parts = []
    for survey in ('dr17', 'dr8'):
        s = sdss_stars(survey); g = gaia(survey)
        d = s.merge(g, on=['PLATE', 'MJD', 'FIBERID'], how='left')
        assert len(d) == len(s)
        if survey == 'dr8':
            for b in 'ugriz':
                d[f'PSFMAG_{b.upper()}'] = d[f'PSFMAG_{b.upper()}_CSV']       # specObj has no PSF mags: take them from the crossmatch
        d = d.drop(columns=[f'PSFMAG_{b.upper()}_CSV' for b in 'ugriz'])
        parts.append(d)
        print(f'{survey}: primary stars {len(s):,}  Gaia match {d.GAIA_SOURCE_ID.notna().sum():,}  finite parallax {d.GAIA_PARALLAX.notna().sum():,}  {time.time()-t0:.0f}s', flush=True)
    df = pd.concat(parts, ignore_index=True)
    df['GAIA_PLX'] = df.GAIA_PARALLAX.notna()
    df['PLX_CORR'] = df.GAIA_PARALLAX + 0.017
    df['QUALITY'] = df.GAIA_PLX & (df.GAIA_PARALLAX_OVER_ERROR > 10) & (df.GAIA_RUWE < 1.4) & (df.GAIA_NON_SINGLE_STAR == 0) \
                    & (df.GAIA_PHOT_VARIABLE_FLAG != 'VARIABLE') & (df.PLX_CORR > 0)
    with np.errstate(all='ignore'):
        df['M_G'] = df.GAIA_PHOT_G_MEAN_MAG + 5 * np.log10(df.PLX_CORR.where(df.PLX_CORR > 0) / 100.0)
        df['M_G0'] = df.M_G - df.GAIA_AG_GSPPHOT.fillna(0); df['BPRP0'] = df.GAIA_BP_RP - df.GAIA_EBPMINRP_GSPPHOT.fillna(0)
    df['RV_PIPE_KMS'] = df.Z * CKMS
    df.to_parquet(OUT + '.tmp'); os.replace(OUT + '.tmp', OUT)
    print(f'written {OUT}: {len(df):,} rows, {df.shape[1]} columns, {time.time()-t0:.0f}s\n')
    for survey, d in df.groupby('SURVEY'):
        p = d[d.GAIA_PLX]; q = d[d.QUALITY]
        print(f'== {survey}: {len(d):,} primary stars; with Gaia parallax {len(p):,} (ZWARN=0: {(p.ZWARN == 0).sum():,}); '
              f'quality {len(q):,}; plx S/N>50 {(q.GAIA_PARALLAX_OVER_ERROR > 50).sum():,}; XP {(p.GAIA_HAS_XP_CONTINUOUS == True).sum():,}; '
              f'GSP-Phot [M/H] {p.GAIA_MH_GSPPHOT.notna().sum():,}; ELODIE Teff {(p.ELODIE_TEFF > 0).sum():,}')
        print(f'   S/N>20: {(p.SNR_CAT > 20).sum():,}  S/N>10: {(p.SNR_CAT > 10).sum():,}   SUBCLASS: {p.SUBCLASS1.value_counts().head(8).to_dict()}')
        print(f'   programs: {p.PROGRAMNAME.value_counts().head(6).to_dict()}')
        print(f'   RV km/s 1/50/99: {np.nanpercentile(p.RV_PIPE_KMS, [1, 50, 99]).round(0)}   Gaia RV available: {p.GAIA_RADIAL_VELOCITY.notna().sum():,}')
        print(f'   M_G0 range (quality) 5/50/95: {np.nanpercentile(q.M_G0, [5, 50, 95]).round(2)}  (BP-RP)_0: {np.nanpercentile(q.BPRP0, [5, 50, 95]).round(2)}')


if __name__ == '__main__':
    main()
