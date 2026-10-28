#!/usr/bin/env python3
"""Master quasar table for DESI DR2 (loa) from the pass-Q output of the shared walk, cross-matched with the DESI DR1
quasar catalog (qso_cat_dr1_main_dark_healpix_zlya-v0), plus checks (docs/HANDOVER_quasars.md sec. 8.15 step C):
  1. concatenate ~/data/sdssdesi_passQ/{bright,dark}/*.parquet (selected columns), quality cuts
     (ZCAT_PRIMARY, ZWARN == 0, COADD_FIBERSTATUS == 0, 0 < z < 5), the rest kept and flagged;
  2. Z_SYS recomputed without O I 1304 (qsolib.systemic_z: [O III] peak -> Mg II peak -> pipeline z);
  3. luminosity: absolute flux density at the normalization window (qsolib.abs_norm_flux, fiber scale divided out);
  4. match DR1 QSO catalog by TARGETID: Z, BI_CIV, AI_CIV;
  5. checks: Z_SYS vs DR1 Z; our BAL_CLASS vs DR1 BI_CIV/AI_CIV; slope by z; per-slice EW quantiles;
  6. write ~/data/sdssdesi_products/quasar/qso_master_desi.parquet, qso_checks_desi.txt, qso_ew_quantiles_desi.csv.
Usage: python3 03_qso_table_desi.py
"""
import os, glob, time, json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import fitsio
import qsolib as Q

HOME = os.path.expanduser('~')
PASSQ = f'{HOME}/data/sdssdesi_passQ'
OUT = f'{HOME}/data/sdssdesi_products/quasar'
DR1Q = f'{HOME}/sdssredux/DESI_DR1_zlya/qso_cat_dr1_main_dark_healpix_zlya-v0.fits'
STATUS = f'{HOME}/data/sdssdesi_status/step_qso_table_desi.json'
CKMS = 299792.458
os.makedirs(OUT, exist_ok=True)

KEEP = ['TARGETID', 'PROGRAM', 'FILE', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'ZCAT_PRIMARY', 'COADD_FIBERSTATUS',
        'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'SNR_B', 'SNR_R', 'SNR_Z', 'SNR_ALL', 'SNR_SDSSPIX', 'EBV',
        'CAL_A', 'CAL_B', 'CAL_B_EST', 'CAL_MODE', 'CAL_DG', 'CAL_DR', 'CAL_OK', 'CAL_DZ_CHECK', 'SYN_FIBER_G_CAL', 'SYN_FIBER_R_CAL', 'SYN_FIBER_Z_CAL',
        'DM', 'DN4000', 'ALPHA_LAM', 'ALPHA_LAM_ERR', 'CONT_LOGF3000', 'CONT_NWIN', 'CONT_LAMMIN', 'CONT_LAMMAX', 'CONT_CHI2NU',
        'REST_LAMMIN', 'REST_LAMMAX', 'FOREST_NPIX', 'FOREST_TRANS',
        'BAL_AI_CIV', 'BAL_BI_CIV', 'BAL_NTROUGH_CIV', 'BAL_VMIN_CIV', 'BAL_VMAX_CIV', 'BAL_AI_SIIV', 'BAL_BI_SIIV', 'BAL_AI_MGII', 'BAL_AI_ALIII', 'BAL_CLASS',
        'Z_MGII', 'Z_MGII_PK', 'Z_OIII', 'Z_OIII_PK', 'Z_CIV', 'Z_CIV_PK', 'Z_HBETA', 'Z_HBETA_PK', 'Z_OI_PK',
        'OII_SUM_EW', 'OII_SUM_EW_ERR', 'OIII_5007_EW', 'OIII_5007_EW_ERR', 'HBETA_EW', 'HALPHA_EW']
LINES = ['OVI', 'LYA', 'NV', 'OI_1304', 'SIIV', 'CIV', 'OIII_1663', 'ALIII', 'CIII', 'MGII', 'HBETA_B', 'HALPHA_B', 'OIII_N', 'OI_6300']
for l in LINES:
    KEEP += [f'{l}_EW', f'{l}_EW_ERR', f'{l}_SNR', f'{l}_SIGMA', f'{l}_VOFF', f'{l}_VPEAK']

# primary / oxygen line per redshift (grid v3, handover 8.7-8.9)
PRIMARY = [(0.0, 0.4, 'HALPHA_B'), (0.4, 0.9, 'HBETA_B'), (0.9, 2.3, 'MGII'), (2.3, 5.0, 'CIV')]
OXYGEN = [(0.0, 0.9, 'OIII_5007'), (0.9, 1.6, 'OII_SUM'), (1.6, 2.6, 'OI_1304'), (2.6, 5.0, 'OVI')]
ZEDGES = np.r_[np.arange(0, 1.01, 0.1), np.arange(1.2, 3.01, 0.2), [3.5, 4.0, 4.5, 5.0]]


def status(done, total, st, note=''):
    json.dump(dict(done=done, total=total, status=st, note=note, started=T0), open(STATUS + '.tmp', 'w')); os.replace(STATUS + '.tmp', STATUS)


def main():
    global T0
    T0 = time.time(); out = []
    files = sorted(glob.glob(f'{PASSQ}/*/*.parquet'))
    status(0, len(files) + 3, 'running', 'reading pass-Q files')
    parts = []; cols = None
    for k, f in enumerate(files):
        if pq.read_metadata(f).num_rows == 0:
            continue
        if cols is None:
            names = pq.read_schema(f).names; cols = [c for c in KEEP if c in names]
        parts.append(pd.read_parquet(f, columns=cols))
        if k % 1000 == 0:
            status(k, len(files) + 3, 'running', f'{k} files read')
    d = pd.concat(parts, ignore_index=True); del parts
    d['SURVEY'] = 'desi'
    d['QCUT'] = d.ZCAT_PRIMARY.astype(bool) & (d.ZWARN == 0) & (d.COADD_FIBERSTATUS == 0) & (d.Z > 0) & (d.Z < 5)
    out.append(f'pass-Q DESI rows: {len(d):,} (dark {(d.PROGRAM=="dark").sum():,}, bright {(d.PROGRAM=="bright").sum():,}); '
               f'ZCAT_PRIMARY {d.ZCAT_PRIMARY.sum():,}; after cuts {d.QCUT.sum():,}')
    d['Z_SYS'], d['Z_SYS_SRC'] = Q.systemic_z(d)
    for k, v in Q.abs_norm_flux(d, fiber_scale=True).items():
        d[k] = v
    status(len(files), len(files) + 3, 'running', 'matching DR1 QSO catalog')
    # --- DR1 quasar catalog (main dark)
    c = fitsio.read(DR1Q, columns=['TARGETID', 'Z', 'BI_CIV', 'AI_CIV'])
    cd = pd.DataFrame({'TARGETID': c['TARGETID'].astype(np.int64), 'DR1Q_Z': c['Z'].astype(float),
                       'DR1Q_BI_CIV': c['BI_CIV'].astype(float), 'DR1Q_AI_CIV': c['AI_CIV'].astype(float)}).drop_duplicates('TARGETID')
    d = d.merge(cd, on='TARGETID', how='left')
    d['IN_DR1Q'] = d.DR1Q_Z.notna()
    out.append(f'matched to DR1 QSO catalog: {d.IN_DR1Q.sum():,} of {len(d):,} (catalog {len(cd):,}; dark after cuts '
               f'{d[d.QCUT & (d.PROGRAM=="dark")].IN_DR1Q.mean()*100:.1f} %)')
    d.to_parquet(f'{OUT}/qso_master_desi.parquet')
    status(len(files) + 1, len(files) + 3, 'running', 'checks')
    g = d[d.QCUT & (d.SNR_ALL > 1)]
    out.append(f'QCUT & SNR_ALL>1 (per 0.8 A pixel; ~S/N 5 per SDSS pixel is SNR_SDSSPIX): {len(g):,}; SNR_SDSSPIX>5: {(d.QCUT & (d.SNR_SDSSPIX>5)).sum():,}')
    g = d[d.QCUT & (d.SNR_SDSSPIX > 5)]
    # --- check 1: Z_SYS vs DR1 Z
    out.append('\n== Z_SYS vs DR1Q Z (km/s), SNR_SDSSPIX>5, |dv|<3000 for the statistics')
    for src, name in ((1, '[O III] peak'), (2, 'Mg II peak'), (0, 'pipeline')):
        s = g[(g.Z_SYS_SRC == src) & g.IN_DR1Q]
        dv = (s.Z_SYS - s.DR1Q_Z) / (1 + s.DR1Q_Z) * CKMS; dv = dv[np.isfinite(dv) & (np.abs(dv) < 3000)]
        if len(dv) > 50:
            out.append(f'  {name:13s}: n={len(dv):8,} median {dv.median():7.0f}  MAD {1.48*np.abs(dv-dv.median()).median():6.0f}  |dv|>1000: {(np.abs(dv)>1000).mean()*100:.1f} %')
    out.append('  Z_SYS_SRC counts (all after cuts): ' + ', '.join(f'{k}: {v:,}' for k, v in d[d.QCUT].Z_SYS_SRC.value_counts().sort_index().items()))
    # --- check 2: BAL vs DR1 catalog
    out.append('\n== BAL: our BAL_CLASS (rows) vs DR1Q BI_CIV>0 (cols), SNR_SDSSPIX>5, 1.6<z<5, in DR1Q')
    b = g[g.IN_DR1Q & (g.Z > 1.6) & (g.Z < 5) & (g.BAL_CLASS >= 0)]
    out.append(pd.crosstab(b.BAL_CLASS, b.DR1Q_BI_CIV > 0, margins=True).to_string())
    bi = b.DR1Q_BI_CIV > 0; ai0 = (b.DR1Q_AI_CIV <= 0) & (b.DR1Q_BI_CIV <= 0)
    out.append(f'  DR1Q BI_CIV>0: {bi.sum():,} ({bi.mean()*100:.1f} %); of those our class>0: {(b.BAL_CLASS[bi]>0).mean()*100:.1f} %, class 2: {(b.BAL_CLASS[bi]==2).mean()*100:.1f} %')
    out.append(f'  DR1Q AI=BI=0: {ai0.sum():,}; of those our class>0 (false positives): {(b.BAL_CLASS[ai0]>0).mean()*100:.1f} %')
    # --- check 3: luminosity and slope by z
    out.append('\n== by z slice: n, median log L_lam at the normalization window (erg/s/A), window inside fit, median ALPHA_LAM (SNR_SDSSPIX>5)')
    rows = []
    for a, bz in zip(ZEDGES[:-1], ZEDGES[1:]):
        s = g[(g.Z >= a) & (g.Z < bz)]
        rows.append({'z': f'{a:.1f}-{bz:.1f}', 'n': len(s), 'lam_norm': s.LAM_NORM.median() if len(s) else np.nan,
                     'logL_med': round(s.LOGL_NORM.median(), 2) if len(s) > 20 else np.nan,
                     'in_fit': round(s.NORM_IN_FIT.mean(), 3) if len(s) else np.nan,
                     'alpha_med': round(s.ALPHA_LAM.median(), 2) if len(s) > 20 else np.nan})
    out.append(pd.DataFrame(rows).to_string(index=False))
    # --- check 4: EW quantiles per slice
    out.append('\n== per-slice quantiles (16/50/84 %) of the primary-line EW, oxygen-line EW and slope, SNR_SDSSPIX>5, line detected >=3 sigma')
    qrows = []
    for a, bz in zip(ZEDGES[:-1], ZEDGES[1:]):
        s = g[(g.Z >= a) & (g.Z < bz)]
        zc = 0.5 * (a + bz)
        pl = [l for lo, hi, l in PRIMARY if lo <= zc < hi][0]; ol = [l for lo, hi, l in OXYGEN if lo <= zc < hi][0]
        r = {'z_lo': a, 'z_hi': bz, 'n': len(s), 'primary': pl, 'oxygen': ol}
        for tag, l in (('p', pl), ('o', ol)):
            ew = s[f'{l}_EW']; err = s[f'{l}_EW_ERR']; det = ew[(ew > 0) & (err > 0) & (ew / err >= 3)]
            r[f'{tag}_ndet'] = len(det)
            for qn in (16, 50, 84):
                r[f'{tag}_q{qn}'] = round(np.percentile(det, qn), 2) if len(det) > 50 else np.nan
        al = s.ALPHA_LAM.dropna()
        for qn in (16, 50, 84):
            r[f'alpha_q{qn}'] = round(np.percentile(al, qn), 2) if len(al) > 50 else np.nan
        bal = s[s.BAL_CLASS >= 0]
        r['bal_frac'] = round((bal.BAL_CLASS > 0).mean(), 3) if len(bal) > 50 else np.nan
        qrows.append(r)
    qdf = pd.DataFrame(qrows); qdf.to_csv(f'{OUT}/qso_ew_quantiles_desi.csv', index=False)
    out.append(qdf.to_string(index=False))
    txt = '\n'.join(out)
    open(f'{OUT}/qso_checks_desi.txt', 'w').write(txt); print(txt)
    status(len(files) + 3, len(files) + 3, 'done', f'{len(d):,} rows')
    print(f'\nwrote {OUT}/qso_master_desi.parquet ({len(d):,} rows, {len(d.columns)} cols), {time.time()-T0:.0f}s')


if __name__ == '__main__':
    main()
