#!/usr/bin/env python3
"""Step 5 (Nao, 2026-10-24 night): formulate the HR grid and test whether the stars of a cell have similar spectra;
report the cells that different stellar groups share.

Input: pass tables of steps 2-3 (SDSS: ~/data/sdssdesi_passS_sdss/*/*.parquet; DESI: ~/data/sdssdesi_passS/*/*.parquet,
joined with the Gaia quality table gaia_desi_stars_hr.parquet for [M/H]).
Selection: quality parallax (QUALITY for SDSS; membership of the DESI quality table), ZWARN = 0, S/N > --snmin.
HR coordinates from the pass columns: M_G0 = G - DM - A_G_USED, (BP-RP)_0 = BP-RP - E(BP-RP)_USED (3D extinction when
step3b has run, else 0.86 x SFD, flagged in EBV_SRC).
Grid (default): M_G0 0.5 mag x (BP-RP)_0 0.1 x [M/H] bins (< -1, -1..-0.5, -0.5..0, > 0; -1 = unknown); [M/H] from
GAIA_MH_GSPPHOT (both) or ELODIE_FEH (SDSS, --feh elodie).
Per cell with >= --nmin stars, "alike" is tested on RESIDUALS: each quantity is first fitted with a plane in (M_G0, (BP-RP)_0)
inside the cell (the smooth temperature/gravity gradient across a 0.5 x 0.1 cell is not "different groups"), the colors are
first corrected with the step-2b tilt (where TILT_OK) and re-dereddened with the 3D E(B-V); then
  * index spread: robust sigma of the residuals / median measurement error -> R; a cell is flagged when R exceeds 2x the median
    R of all cells for that index (the median absorbs the systematic floor of the index errors)
  * bimodality coefficient BC = (skew^2 + 1) / (kurtosis + 3) of the residuals; BC > 0.555 = two groups
  * continuum: residual scatter of M_G_SYN0 - M_G0 (> 0.15 mag: binaries, reddening, calibration) and of the corrected colors
    (> 2x the median cell)
  * Halpha emitters (EW < -1 A at 3 sigma): fraction > 0.10
  * [Fe/H] spread: robust sigma > 0.3 dex or bimodal; the lowest [M/H] bin is split at -1.5 by default (--feh-edges)
Each failing test names a probable cause and a split axis (FEH, HALPHA, BALMER; INDEX/COLOR/CONT are reported, not split).
Na D carries interstellar absorption (distance and dust dependent) and is reported but not used for the verdict. Output: the cell table
(CSV, one row per cell with the statistics, flags, composition by SUBTYPE/SUBCLASS and [M/H]) and a printed report.
Usage:
  python3 step5_hr_cells_check.py --survey sdss [--passdir DIR] [--nmin 100] [--snmin 5] [--dmg 0.5 --dcol 0.1] [--tag sdss]
"""
import os, sys, glob, argparse, time
import numpy as np, pandas as pd
from scipy import stats as st

CACHE = os.path.expanduser('~/data/sdssdesi_cache/stars'); os.makedirs(CACHE, exist_ok=True)
PASSDIR = {'sdss': os.path.expanduser('~/data/sdssdesi_passS_sdss'), 'desi': os.path.expanduser('~/data/sdssdesi_passS')}
FEH_EDGES = [-9, -1.5, -1.0, -0.5, 0.0, 9]
INDICES = ['CAK_EW', 'HDELTA_A_EW', 'G4300_EW', 'HBETA_LK_EW', 'MGB_EW', 'FE5270_EW', 'HALPHA_ST_EW', 'CAT_EW', 'HEI_5876_EW']
REPORT_ONLY = ['NAD_EW']
RATIOS = ['TIO5', 'CAH2', 'CAH3']
R_REL_MAX, BC_MAX, CONT_MAX, COLOR_REL_MAX, HA_MAX, FEH_MAX = 2.0, 0.555, 0.15, 2.0, 0.10, 0.3


def load(survey, passdir):
    files = sorted(glob.glob(os.path.join(passdir, '*', '*.parquet')))
    import pyarrow.parquet as pq
    files = [f for f in files if 'FILE' in pq.ParquetFile(f).schema_arrow.names]          # skip the empty per-file tables
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    if survey == 'desi':
        hr = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet'), columns=['targetid', 'mh_gspphot'])
        df = df.merge(hr.rename(columns={'targetid': 'TARGETID', 'mh_gspphot': 'GAIA_MH_GSPPHOT'}), on='TARGETID', how='inner')
        df['QUALITY'] = True; df['SUBCLASS1'] = df['SUBTYPE'].str[:1]
        df['SID'] = df['TARGETID'].astype(str)
    else:
        df['SID'] = df.PLATE.astype(str) + '-' + df.MJD.astype(str) + '-' + df.FIBERID.astype(str)
    return df, len(files)


def hr_cells(df, dmg, dcol, feh_col):
    """cell keys from the pass columns (shared with step6_stack.py)."""
    with np.errstate(all='ignore'):
        mg0 = df['GAIA_PHOT_G_MEAN_MAG'] - df['DM'] - df['A_G_USED']
        col0 = df['GAIA_BP_RP'] if 'GAIA_BP_RP' in df else df['GAIA_PHOT_BP_MEAN_MAG'] - df['GAIA_PHOT_RP_MEAN_MAG']
        col0 = col0 - df['EBPRP_USED']
    feh = df[feh_col] if feh_col in df else pd.Series(np.nan, index=df.index)
    out = pd.DataFrame({'MG0': mg0, 'COL0': col0, 'FEH': feh,
                        'IA': np.floor(mg0 / dmg), 'IB': np.floor(col0 / dcol),
                        'IF': pd.cut(feh, FEH_EDGES, labels=False).fillna(-1).astype(int)}, index=df.index)
    out['CELL'] = out.IA.astype('Int64').astype(str) + '_' + out.IB.astype('Int64').astype(str) + '_' + out.IF.astype(str)
    return out


def rsig(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 3 else np.nan


def bimod(x):
    """bimodality coefficient; NaN for discrete variables (< 15 distinct values, e.g. ELODIE [Fe/H] on its template grid)."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 20 or len(np.unique(np.round(x, 3))) < 15: return np.nan
    return (st.skew(x)**2 + 1) / (st.kurtosis(x) + 3)


def resid(y, mg, col):
    """residual of y after a plane fit in (M_G0, color) inside the cell (removes the smooth gradient across the cell)."""
    y = np.asarray(y, float); ok = np.isfinite(y) & np.isfinite(mg) & np.isfinite(col)
    out = np.full(len(y), np.nan)
    if ok.sum() < 10:
        out[ok] = y[ok] - np.median(y[ok]); return out
    A = np.stack([np.ones(ok.sum()), mg[ok] - mg[ok].mean(), col[ok] - col[ok].mean()], axis=1)
    coef, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
    out[ok] = y[ok] - A @ coef
    return out


def cell_stats(g):
    r = {'N': len(g)}
    r['MG0_MED'], r['COL0_MED'], r['FEH_MED'] = g.MG0.median(), g.COL0.median(), g.FEH.median()
    r['SNR_MED'] = g.SNR_ALL.median()
    mg = g.MG0.to_numpy(float); col = g.COL0.to_numpy(float)
    for c in INDICES + RATIOS + REPORT_ONLY:
        err = c.replace('_EW', '_EW_ERR') if c.endswith('_EW') else c + '_ERR'
        res = resid(g[c], mg, col)
        s = rsig(res); e = np.nanmedian(g[err]) if err in g else np.nan
        r[f'R_{c}'] = s / e if (np.isfinite(e) and e > 0) else np.nan
        r[f'BC_{c}'] = bimod(res)
    r['SIG_ABSG'] = rsig(resid(g.M_G_SYN0 - g.MG0, mg, col))
    r['SIG_CBLUE'] = rsig(resid(g.CBLUE_C, mg, col)); r['SIG_CRED'] = rsig(resid(g.CRED_C, mg, col))
    r['BC_CBLUE'] = bimod(resid(g.CBLUE_C, mg, col))
    ha = (g.HALPHA_ST_EW < -1) & (g.HALPHA_ST_EW / g.HALPHA_ST_EW_ERR < -3)
    r['FRAC_HA_EM'] = float(ha.mean())
    fres = resid(g.FEH, mg, col)
    r['SIG_FEH'] = rsig(fres); r['BC_FEH'] = bimod(fres); r['FRAC_FEH_LT_M1'] = float((g.FEH < -1).mean()) if g.FEH.notna().any() else np.nan
    comp = g.SUBCLASS1.value_counts(normalize=True)
    r['COMP'] = ' '.join(f'{k}:{v:.2f}' for k, v in comp.head(4).items())
    return r


def verdicts(tab):
    """second pass over all cells: spreads relative to the typical cell, bimodality, causes, split axis."""
    med = {c: np.nanmedian(tab[f'R_{c}']) for c in INDICES + RATIOS}
    med_cb, med_cr = np.nanmedian(tab.SIG_CBLUE), np.nanmedian(tab.SIG_CRED)
    causes_all, split_all, rmax_all, rmax_idx, bimod_all = [], [], [], [], []
    for _, r in tab.iterrows():
        causes = []; rel = {}
        for c in INDICES + RATIOS:
            if np.isfinite(r[f'R_{c}']) and med[c] > 0:
                rel[c] = r[f'R_{c}'] / med[c]
        bim = [c for c in INDICES + RATIOS if np.isfinite(r[f'BC_{c}']) and r[f'BC_{c}'] > BC_MAX and rel.get(c, 0) > 1.3]
        worst = max(rel, key=rel.get) if rel else ''
        rmax_all.append(rel.get(worst, np.nan)); rmax_idx.append(worst); bimod_all.append('+'.join(bim))
        if r.FRAC_HA_EM > HA_MAX: causes.append('HALPHA')
        if (np.isfinite(r.SIG_FEH) and r.SIG_FEH > FEH_MAX) or (np.isfinite(r.BC_FEH) and r.BC_FEH > BC_MAX): causes.append('FEH')
        if ('HBETA_LK_EW' in bim or 'HDELTA_A_EW' in bim) and r.COL0_MED < 0.9: causes.append('BALMER')
        if np.isfinite(r.SIG_ABSG) and r.SIG_ABSG > CONT_MAX: causes.append('CONT')
        if (np.isfinite(r.SIG_CBLUE) and r.SIG_CBLUE > COLOR_REL_MAX * med_cb) or (np.isfinite(r.SIG_CRED) and r.SIG_CRED > COLOR_REL_MAX * med_cr) or \
           (np.isfinite(r.BC_CBLUE) and r.BC_CBLUE > BC_MAX): causes.append('COLOR')
        if bim and not any(x in causes for x in ('FEH', 'HALPHA', 'BALMER')): causes.append('BIMODAL')
        if rel and rel[worst] > R_REL_MAX and not causes: causes.append('INDEX')
        causes_all.append('+'.join(causes))
        split_all.append(next((c for c in causes if c in ('FEH', 'HALPHA', 'BALMER')), 'NONE'))
    tab['R_MAX'] = rmax_all; tab['R_MAX_INDEX'] = rmax_idx; tab['BIMODAL'] = bimod_all
    tab['CAUSES'] = causes_all; tab['HOMOGENEOUS'] = [c == '' for c in causes_all]; tab['SPLIT_AXIS'] = split_all
    return tab


def _progress(name, done, total, status='done', note=''):
    """progress file of the monitors (programs_monitor/monitor.py, monitor_stars.py)."""
    try:
        import json; S = os.path.expanduser('~/data/sdssdesi_status'); os.makedirs(S, exist_ok=True)
        json.dump({'done': int(done), 'total': int(total), 'status': status, 'note': note}, open(f'{S}/step_star_{name}.json', 'w'))
    except Exception:
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='sdss', choices=['sdss', 'desi'])
    ap.add_argument('--passdir'); ap.add_argument('--nmin', type=int, default=100); ap.add_argument('--snmin', type=float, default=5.0)
    ap.add_argument('--dmg', type=float, default=0.5); ap.add_argument('--dcol', type=float, default=0.1)
    ap.add_argument('--feh', default='gspphot', choices=['gspphot', 'elodie']); ap.add_argument('--tag')
    ap.add_argument('--feh-edges', default='-1.5,-1,-0.5,0')
    a = ap.parse_args(); t0 = time.time()
    global FEH_EDGES; FEH_EDGES = [-9] + [float(v) for v in a.feh_edges.split(',')] + [9]
    df, nf = load(a.survey, a.passdir or PASSDIR[a.survey])
    feh_col = 'GAIA_MH_GSPPHOT' if a.feh == 'gspphot' else 'ELODIE_FEH'
    sel = df.QUALITY.fillna(False).astype(bool) & (df.ZWARN == 0) & (df.SNR_ALL > a.snmin) & df.CAL_OK.fillna(False).astype(bool) & np.isfinite(df.DM)
    df = df[sel].reset_index(drop=True)
    cells = hr_cells(df, a.dmg, a.dcol, feh_col); df = pd.concat([df, cells], axis=1)
    # colors: raw -> 3D dereddening (EBV_USED already is the 3D value when step3b ran before the pass; here we recompute from
    # EBV_3D when the table has it) -> step-2b tilt correction where TILT_OK
    import importlib.util, sys as _s; _s.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import starlib
    f99 = {k: float(starlib.passB.extinction.fitzpatrick99(np.array([0.5 * (lo + hi)]), 1.0, 3.1)[0]) for k, (lo, hi) in starlib.COLOR_BANDS.items()}
    ebv = df['EBV_3D'].to_numpy(float) if 'EBV_3D' in df else np.full(len(df), np.nan)
    ebv = np.where(np.isfinite(ebv), ebv, df['EBV_USED'].to_numpy(float)); av = 3.1 * ebv
    df['CBLUE_C'] = df.C_BLUE_RAW - av * (f99['4100'] - f99['5600']); df['CRED_C'] = df.C_RED_RAW - av * (f99['5600'] - f99['8100'])
    sp = os.path.join(CACHE, f'spcal_{a.survey}.parquet')
    if os.path.exists(sp):
        s2 = pd.read_parquet(sp, columns=['SID', 'TILT', 'TILT_OK']); df = df.merge(s2, on='SID', how='left')
        tl = np.where(df.TILT_OK.fillna(False).astype(bool), df.TILT.to_numpy(float), 0.0)
        x = {k: np.log10(0.5 * (lo + hi) / 5500.0) for k, (lo, hi) in starlib.COLOR_BANDS.items()}
        df['CBLUE_C'] = df.CBLUE_C + tl * (x['4100'] - x['5600']); df['CRED_C'] = df.CRED_C + tl * (x['5600'] - x['8100'])
        print(f'step-2b tilt applied to the colors of {int((tl != 0).sum()):,} stars', flush=True)
    print(f'{a.survey}: {nf} files, {len(df):,} quality stars with ZWARN=0, S/N>{a.snmin}, CAL_OK; grid {a.dmg} x {a.dcol} x [M/H]({a.feh}); {time.time()-t0:.0f}s')
    rows = []
    for key, g in df.groupby('CELL'):
        if len(g) < a.nmin: continue
        r = cell_stats(g); r['CELL'] = key; rows.append(r)
    tab = verdicts(pd.DataFrame(rows)).sort_values(['MG0_MED', 'COL0_MED']).reset_index(drop=True)
    tag = a.tag or f'{a.survey}_{a.dmg}x{a.dcol}_{a.feh}'
    out = os.path.join(CACHE, f'cells_{tag}.csv'); tab.to_csv(out, index=False)
    df[['SID', 'CELL', 'MG0', 'COL0', 'FEH', 'SNR_ALL']].to_parquet(os.path.join(CACHE, f'members_{tag}.parquet'))
    n_in = int(tab.N.sum()) if len(tab) else 0
    print(f'cells with >= {a.nmin} stars: {len(tab)} ({n_in:,} stars, {n_in/max(len(df),1)*100:.1f} % of the sample); homogeneous {int(tab.HOMOGENEOUS.sum()) if len(tab) else 0}, mixed {int((~tab.HOMOGENEOUS).sum()) if len(tab) else 0}')
    if len(tab):
        print('mixed cells by cause:', tab[~tab.HOMOGENEOUS].CAUSES.value_counts().to_dict())
        print('median over cells (residual spreads): R_MAX/median %.2f (index %s most often), sigma(absG) %.3f, sigma(C_BLUE) %.3f, Halpha-emitter fraction %.3f, sigma[M/H] %.2f' % (
            tab.R_MAX.median(), tab.R_MAX_INDEX.mode().iloc[0] if tab.R_MAX_INDEX.notna().any() else '', tab.SIG_ABSG.median(), tab.SIG_CBLUE.median(), tab.FRAC_HA_EM.median(), tab.SIG_FEH.median()))
        print('\n== cells shared by different stellar groups (mixed), largest first:')
        print('  M_G0   (BP-RP)0  [M/H]   N   S/N  Rrel(index)       sig_absG  sig_Cblue  f(Ha em)  sig[M/H]  composition            bimodal in          causes')
        for _, r in tab[~tab.HOMOGENEOUS].sort_values('N', ascending=False).head(40).iterrows():
            print(f"  {r.MG0_MED:5.2f}  {r.COL0_MED:6.2f}  {r.FEH_MED if np.isfinite(r.FEH_MED) else float('nan'):5.2f}  {int(r.N):5d} {r.SNR_MED:5.1f}  {r.R_MAX:4.1f} {r.R_MAX_INDEX:12s}  {r.SIG_ABSG:7.3f}  {r.SIG_CBLUE:8.3f}  {r.FRAC_HA_EM:7.3f}  {r.SIG_FEH:7.2f}  {r.COMP:22s} {str(r.BIMODAL)[:18]:18s} {r.CAUSES}")
        print('\n== homogeneous cells (alike stars), largest first:')
        for _, r in tab[tab.HOMOGENEOUS].sort_values('N', ascending=False).head(15).iterrows():
            print(f"  {r.MG0_MED:5.2f}  {r.COL0_MED:6.2f}  {r.FEH_MED if np.isfinite(r.FEH_MED) else float('nan'):5.2f}  {int(r.N):5d} {r.SNR_MED:5.1f}  {r.R_MAX:4.1f} {r.R_MAX_INDEX:12s}  {r.SIG_ABSG:7.3f}  {r.SIG_CBLUE:8.3f}  {r.COMP}")
    print('written', out)
    _progress('cells', 1, 2 if a.survey == 'sdss' else 2, 'running' if a.survey == 'sdss' else 'done', f'{a.survey}: {len(tab)} cells, {int((~tab.HOMOGENEOUS).sum()) if len(tab) else 0} mixed')


if __name__ == '__main__':
    main()
