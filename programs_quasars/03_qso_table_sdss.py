#!/usr/bin/env python3
"""Master quasar table for SDSS (DR17 + DR8) from the pass-Q output, cross-matched with DR16Q, plus the checks that
needed real measurements (docs/HANDOVER_quasars.md sec. 8.15 step C):
  1. concatenate ~/data/sdssdesi_passQ_sdss/{dr17,dr8}/*.parquet (selected columns), apply the quality cuts
     (ZWARN == 0, PLATEQUALITY == good, 0 < z < 5), keep the rest flagged;
  2. match DR16Q_v4 by (PLATE, MJD, FIBERID): BAL_PROB, BI_CIV, AI_CIV, BI_SIIV, AI_SIIV, Z_PCA, Z_MGII, Z_CIV, Z_HBETA, Z_VI;
  3. checks: Z_SYS vs DR16Q Z_MGII / Z_PCA; our BAL_CLASS vs BAL_PROB / BI_CIV; slope by survey and z; EW distributions
     per redshift slice -> proposed bin centres (median) for the six-axis grid;
  4. write ~/data/sdssdesi_products/quasar/qso_master_sdss.parquet and ~/data/sdssdesi_products/quasar/qso_checks_sdss.txt
     (+ per-slice EW/slope quantiles csv for the grid: qso_ew_quantiles_sdss.csv).
Usage: python3 03_qso_table_sdss.py
"""
import os, glob, time, json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import fitsio

HOME = os.path.expanduser('~')
PASSQ = f'{HOME}/data/sdssdesi_passQ_sdss'
OUT = f'{HOME}/data/sdssdesi_products/quasar'
DR16Q = f'{HOME}/sdssredux/DR16Q/DR16Q_v4.fits'
STATUS = f'{HOME}/data/sdssdesi_status/step_qso_table_sdss.json'
CKMS = 299792.458
os.makedirs(OUT, exist_ok=True)

KEEP = ['PLATE', 'MJD', 'FIBERID', 'SURVEY', 'THING_ID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZERR', 'ZWARN', 'PLATEQUALITY', 'PROGRAMNAME',
        'SNR_ALL', 'SNR_CAT', 'SNR_B', 'SNR_R', 'EBV', 'CAL_A', 'CAL_B_EST', 'CAL_DG', 'CAL_DR', 'CAL_OK', 'CAL_MODE', 'CAL_DG_FIBER', 'CAL_DR_FIBER',
        'PSFMAG_I', 'FIBERMAG_I', 'M_I_OBS_DERED', 'M_I_SRC', 'DM', 'MABS_I_Z2', 'DN4000',
        'ALPHA_LAM', 'ALPHA_LAM_ERR', 'CONT_LOGF3000', 'CONT_NWIN', 'CONT_LAMMIN', 'CONT_LAMMAX', 'CONT_CHI2NU',
        'REST_LAMMIN', 'REST_LAMMAX', 'FOREST_NPIX', 'FOREST_TRANS',
        'BAL_AI_CIV', 'BAL_BI_CIV', 'BAL_NTROUGH_CIV', 'BAL_VMIN_CIV', 'BAL_VMAX_CIV', 'BAL_AI_SIIV', 'BAL_BI_SIIV', 'BAL_AI_MGII', 'BAL_AI_ALIII', 'BAL_CLASS',
        'Z_SYS', 'Z_SYS_SRC', 'Z_MGII', 'Z_MGII_PK', 'Z_OIII', 'Z_OIII_PK', 'Z_CIV', 'Z_CIV_PK', 'Z_HBETA', 'Z_HBETA_PK', 'Z_OI_PK',
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
    parts = []
    for k, f in enumerate(files):
        md = pq.read_metadata(f)
        if md.num_rows == 0:
            continue
        cols = [c for c in KEEP if c in pq.read_schema(f).names]
        parts.append(pd.read_parquet(f, columns=cols))
        if k % 500 == 0:
            status(k, len(files) + 3, 'running', f'{k} files read')
    d = pd.concat(parts, ignore_index=True); del parts
    d['QCUT'] = (d.ZWARN == 0) & (d.PLATEQUALITY == 'good') & (d.Z > 0) & (d.Z < 5)
    out.append(f'pass-Q SDSS rows: {len(d):,} (DR17 {(d.SURVEY=="dr17").sum():,}, DR8 {(d.SURVEY=="dr8").sum():,}); after cuts {d.QCUT.sum():,}')
    status(len(files), len(files) + 3, 'running', 'matching DR16Q')
    # --- DR16Q
    q = fitsio.read(DR16Q, columns=['PLATE', 'MJD', 'FIBERID', 'Z', 'Z_PCA', 'Z_MGII', 'Z_CIV', 'Z_HBETA', 'Z_VI', 'BAL_PROB', 'BI_CIV', 'AI_CIV', 'BI_SIIV', 'AI_SIIV', 'IS_QSO_FINAL'])
    qd = pd.DataFrame({'PLATE': q['PLATE'].astype(np.int32), 'MJD': q['MJD'].astype(np.int32), 'FIBERID': q['FIBERID'].astype(np.int16)})
    for c in ['Z', 'Z_PCA', 'Z_MGII', 'Z_CIV', 'Z_HBETA', 'Z_VI', 'BAL_PROB', 'BI_CIV', 'AI_CIV', 'BI_SIIV', 'AI_SIIV']:
        qd['DR16Q_' + c] = q[c].astype(float)
    qd['DR16Q_IS_QSO_FINAL'] = q['IS_QSO_FINAL'].astype(np.int8)
    qd = qd.drop_duplicates(['PLATE', 'MJD', 'FIBERID'])
    d = d.merge(qd, on=['PLATE', 'MJD', 'FIBERID'], how='left')
    d['IN_DR16Q'] = d.DR16Q_Z.notna()
    out.append(f'matched to DR16Q: {d.IN_DR16Q.sum():,} of {len(d):,} (DR17 rows: {d[d.SURVEY=="dr17"].IN_DR16Q.mean()*100:.1f} %)')
    d.to_parquet(f'{OUT}/qso_master_sdss.parquet')
    status(len(files) + 1, len(files) + 3, 'running', 'checks')
    g = d[d.QCUT & (d.SNR_ALL > 5)]
    out.append(f'S/N>5 after cuts: {len(g):,} (DR17 {(g.SURVEY=="dr17").sum():,}, DR8 {(g.SURVEY=="dr8").sum():,})')
    # --- check 1: Z_SYS vs DR16Q
    out.append('\n== Z_SYS vs DR16Q (km/s; (z_ours - z_ref)/(1+z) c), S/N>5, |dv|<3000 for the statistics')
    for src, name in ((1, '[O III] peak'), (2, 'Mg II peak'), (3, 'O I 1304 peak')):
        s = g[g.Z_SYS_SRC == src]
        for ref in ['DR16Q_Z_MGII', 'DR16Q_Z_PCA', 'DR16Q_Z_VI']:
            dv = (s.Z_SYS - s[ref]) / (1 + s[ref]) * CKMS; dv = dv[(s[ref] > 0) & np.isfinite(dv)]; dv = dv[np.abs(dv) < 3000]
            if len(dv) > 50:
                out.append(f'  {name:15s} vs {ref:13s}: n={len(dv):7,} median {dv.median():7.0f}  MAD {1.48*np.abs(dv-dv.median()).median():6.0f}  |dv|>1000: {(np.abs(dv)>1000).mean()*100:.1f} %')
    dv = (g.Z - g.DR16Q_Z_PCA) / (1 + g.DR16Q_Z_PCA) * CKMS; dv = dv[(g.DR16Q_Z_PCA > 0)]; dv = dv[np.abs(dv) < 3000]
    out.append(f'  pipeline Z     vs DR16Q_Z_PCA : n={len(dv):7,} median {dv.median():7.0f}  MAD {1.48*np.abs(dv-dv.median()).median():6.0f}')
    # --- check 2: BAL class vs DR16Q
    out.append('\n== BAL: our BAL_CLASS (rows) vs DR16Q BAL_PROB (cols: <0.5, >=0.5) and BI_CIV>0, S/N>5, 1.57<z<5 (C IV in BOSS window), DR17 only')
    b = g[(g.SURVEY == 'dr17') & (g.Z > 1.57) & (g.Z < 5) & g.IN_DR16Q & (g.BAL_CLASS >= 0)]
    tab = pd.crosstab(b.BAL_CLASS, b.DR16Q_BAL_PROB >= 0.5, margins=True); out.append(tab.to_string())
    bi = b.DR16Q_BI_CIV > 0
    out.append(f'  DR16Q BI_CIV>0: {bi.sum():,} ({bi.mean()*100:.1f} %); of those our class>0: {(b.BAL_CLASS[bi]>0).mean()*100:.1f} %, class 2 (BI>0): {(b.BAL_CLASS[bi]==2).mean()*100:.1f} %')
    nb = b.DR16Q_BAL_PROB < 0.1
    out.append(f'  DR16Q BAL_PROB<0.1: {nb.sum():,}; of those our class>0 (false positives): {(b.BAL_CLASS[nb]>0).mean()*100:.1f} %')
    for lo, hi in [(3, 5), (5, 10), (10, 100)]:
        bb = b[(b.SNR_ALL >= lo) & (b.SNR_ALL < hi)]
        if len(bb):
            out.append(f'  S/N {lo}-{hi}: n={len(bb):,}  our BAL fraction {(bb.BAL_CLASS>0).mean()*100:.1f} %  DR16Q BAL_PROB>=0.5 {(bb.DR16Q_BAL_PROB>=0.5).mean()*100:.1f} %  agreement {((bb.BAL_CLASS>0)==(bb.DR16Q_BAL_PROB>=0.5)).mean()*100:.1f} %')
    # --- check 3: slope by survey and z
    out.append('\n== continuum slope ALPHA_LAM (median, S/N>5) by survey and z slice')
    rows = []
    for a, bz in zip(ZEDGES[:-1], ZEDGES[1:]):
        s = g[(g.Z >= a) & (g.Z < bz)]
        r = {'z': f'{a:.1f}-{bz:.1f}'}
        for sv in ('dr17', 'dr8'):
            ss = s[s.SURVEY == sv]; r[f'n_{sv}'] = len(ss); r[f'alpha_{sv}'] = round(ss.ALPHA_LAM.median(), 2) if len(ss) > 20 else np.nan
        rows.append(r)
    out.append(pd.DataFrame(rows).to_string(index=False))
    # --- check 4: EW quantiles per slice for the grid (primary + oxygen lines, slope)
    out.append('\n== per-slice quantiles (16/50/84 %) of the primary-line EW, oxygen-line EW and slope, S/N>5, line detected >=3 sigma')
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
    qdf = pd.DataFrame(qrows); qdf.to_csv(f'{OUT}/qso_ew_quantiles_sdss.csv', index=False)
    out.append(qdf.to_string(index=False))
    txt = '\n'.join(out)
    open(f'{OUT}/qso_checks_sdss.txt', 'w').write(txt); print(txt)
    status(len(files) + 3, len(files) + 3, 'done', f'{len(d):,} rows')
    print(f'\nwrote {OUT}/qso_master_sdss.parquet ({len(d):,} rows, {len(d.columns)} cols), {time.time()-T0:.0f}s')


if __name__ == '__main__':
    main()
