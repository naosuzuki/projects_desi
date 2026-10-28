#!/usr/bin/env python3
"""Six-axis cell labels for the quasar stacking grid (docs/HANDOVER_quasars.md 8.7-8.10, 8.15 step D; ms.tex tab:qsogrid).

Axes (codes per object; -1 = unknown / not measurable -> the object goes to the collapsed cell of that axis):
  1 ZBIN   redshift slice: dz=0.1 (0-1), 0.2 (1-3), 0.5 (3-5)                                 24 slices
  2 LBIN   log10 F_lam,abs at the normalization window (erg/s/cm2/A at 10 pc), 0.4-dex bins (= 1 mag), floor(logF/0.4)
  3 EBIN   primary-line EW, 6 factor-2 bins centred on the measured per-slice median (qso_ew_quantiles_*.csv):
           edges med x [1/4, 1/2, 1, 2, 4]; line Halpha (z<0.4) / Hbeta (0.4-0.9) / Mg II (0.9-2.3) / C IV (>=2.3);
           fine bin only when EW/EW_ERR >= 3, else -1
  4 SBIN   continuum slope alpha_lambda, 6 bins: < -2.0, -2.0..-1.75, -1.75..-1.5, -1.5..-1.25, -1.25..-1.0, > -1.0
  5 OBIN   oxygen-line EW, 4 bins: <1, 1-4, 4-16, >16 A; line [O III] 5007 (z<0.9) / [O II] (0.9-1.6) / O III] 1663
           (1.6-1.9) / O I 1304 (1.9-2.6) / O VI (>=2.6); fine bin only when EW/EW_ERR >= 3, else -1
  6 BBIN   BAL class 0 non-BAL, 1 weak HiBAL, 2 strong HiBAL, 3 LoBAL; -1 no BAL region covered
Cells: CELL = 'z{ZBIN}_L{LBIN}_E{EBIN}_S{SBIN}_O{OBIN}_B{BBIN}'. Merging to the floor (training set = QCUT & S/N > 5):
a six-axis cell below the floor is collapsed along BAL, then oxygen, then slope (code -> -1) until it reaches the
floor; what is still below the floor ends in the (z, L, EW) pool cell. STACK_CELL is the resulting label.
Usage:
  python3 04_cells.py --survey sdss          # ~/data/sdssdesi_products/quasar/qso_master_sdss.parquet
  python3 04_cells.py --survey all           # sdss + desi master tables
Output: qso_cells_<survey>.parquet (per-object labels), qso_stack_cells_<survey>.csv (one row per stack cell with
counts), and a summary printed / written to qso_cells_<survey>.txt; also data/quasar/qso_stack_cells_<survey>.csv.
"""
import os, sys, json, time, argparse
import numpy as np
import pandas as pd

HOME = os.path.expanduser('~')
PROD = f'{HOME}/data/sdssdesi_products/quasar'
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'quasar')
STATUS = f'{HOME}/data/sdssdesi_status/step_qso_cells.json'
ZEDGES = np.r_[np.arange(0, 1.01, 0.1), np.arange(1.2, 3.01, 0.2), [3.5, 4.0, 4.5, 5.0]]
PRIMARY = [(0.0, 0.4, 'HALPHA_B'), (0.4, 0.9, 'HBETA_B'), (0.9, 2.3, 'MGII'), (2.3, 5.01, 'CIV')]
OXYGEN = [(0.0, 0.9, 'OIII_5007'), (0.9, 1.6, 'OII_SUM'), (1.6, 1.9, 'OIII_1663'), (1.9, 2.6, 'OI_1304'), (2.6, 5.01, 'OVI')]
SLOPE_EDGES = [-2.0, -1.75, -1.5, -1.25, -1.0]
OXY_EDGES = [1.0, 4.0, 16.0]
FLOOR = 300
MIN_LINE_SNR = 3.0


LBIN_DEX = 0.4

def line_for(z, table):
    out = np.full(len(z), '', dtype=object)
    for lo, hi, l in table:
        out[(z >= lo) & (z < hi)] = l
    return out


def ew_bins_from_median(ew, err, med):
    """6 factor-2 bins around med (per object), fine only when ew/err >= MIN_LINE_SNR."""
    with np.errstate(all='ignore'):
        r = ew / med
        code = np.digitize(r, [0.25, 0.5, 1.0, 2.0, 4.0])          # 0..5
        ok = np.isfinite(ew) & np.isfinite(err) & (err > 0) & (ew / err >= MIN_LINE_SNR) & np.isfinite(med) & (med > 0)
    return np.where(ok, code, -1).astype(np.int8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='sdss', choices=['sdss', 'desi', 'all'])
    ap.add_argument('--floor', type=int, default=FLOOR)
    a = ap.parse_args()
    t0 = time.time()
    json.dump(dict(done=0, total=4, status='running', note='reading', started=t0), open(STATUS, 'w'))
    tabs = []
    for sv in (['sdss', 'desi'] if a.survey == 'all' else [a.survey]):
        p = f'{PROD}/qso_master_{sv}.parquet'
        if not os.path.exists(p):
            print('missing', p); continue
        cols = ['Z', 'ZWARN', 'SNR_ALL', 'QCUT', 'LOGFABS_NORM', 'LOGL_NORM', 'ALPHA_LAM', 'BAL_CLASS', 'Z_SYS', 'Z_SYS_SRC', 'SURVEY']
        lines = sorted({l for _, _, l in PRIMARY} | {l for _, _, l in OXYGEN})
        for l in lines:
            cols += [f'{l}_EW', f'{l}_EW_ERR']
        ids = ['PLATE', 'MJD', 'FIBERID'] if sv == 'sdss' else ['TARGETID', 'PROGRAM']
        d = pd.read_parquet(p, columns=[c for c in ids + cols] + (['SNR_SDSSPIX'] if sv == 'desi' else []))
        if sv == 'desi':      # S/N per SDSS-size pixel (dlog 1e-4), comparable with the SDSS SNR_ALL
            d['SNR_ALL'] = d.pop('SNR_SDSSPIX')
        tabs.append(d)
    d = pd.concat(tabs, ignore_index=True)
    train = d.QCUT.to_numpy() & (d.SNR_ALL.to_numpy() > 5) & np.isfinite(d.LOGFABS_NORM.to_numpy())
    d['TRAIN'] = train
    z = d.Z.to_numpy(dtype=float)
    # --- axis 1, 2
    d['ZBIN'] = (np.digitize(z, ZEDGES) - 1).astype(np.int8)
    d.loc[(z < 0) | (z >= 5), 'ZBIN'] = -1
    # luminosity (Nao, 2026-10-27): absolute flux density at the normalization window, 0.4-dex bins (= 1 mag)
    d['LBIN'] = np.where(np.isfinite(d.LOGFABS_NORM), np.floor(d.LOGFABS_NORM.fillna(0) / LBIN_DEX), -99).astype(np.int16)
    # --- axis 3: per-slice median of the primary line (training set, detected) -> factor-2 bins
    pl = line_for(z, PRIMARY); ol = line_for(z, OXYGEN)
    d['PRIMARY_LINE'], d['OXYGEN_LINE'] = pl, ol
    ew_p = np.full(len(d), np.nan); err_p = np.full(len(d), np.nan); ew_o = np.full(len(d), np.nan); err_o = np.full(len(d), np.nan)
    for l in set(pl):
        if l:
            m = pl == l; ew_p[m] = d.loc[m, f'{l}_EW']; err_p[m] = d.loc[m, f'{l}_EW_ERR']
    for l in set(ol):
        if l:
            m = ol == l; ew_o[m] = d.loc[m, f'{l}_EW']; err_o[m] = d.loc[m, f'{l}_EW_ERR']
    med = np.full(len(d), np.nan)
    med_table = []
    for zb in range(len(ZEDGES) - 1):
        m = (d.ZBIN.to_numpy() == zb)
        det = m & train & np.isfinite(ew_p) & (ew_p > 0) & (err_p > 0) & (ew_p / err_p >= MIN_LINE_SNR)
        if det.sum() >= 50:
            mv = float(np.median(ew_p[det])); med[m] = mv
        else:
            mv = np.nan
        med_table.append(dict(ZBIN=zb, z_lo=ZEDGES[zb], z_hi=ZEDGES[zb + 1], primary=pl[m][0] if m.any() else '', ew_median=mv, n_det=int(det.sum())))
    d['EBIN'] = ew_bins_from_median(ew_p, err_p, med)
    d['PRIMARY_EW'], d['PRIMARY_EW_MED'] = ew_p, med
    # --- axis 4, 5, 6
    al = d.ALPHA_LAM.to_numpy(dtype=float)
    d['SBIN'] = np.where(np.isfinite(al), np.digitize(al, SLOPE_EDGES), -1).astype(np.int8)
    with np.errstate(all='ignore'):
        oko = np.isfinite(ew_o) & np.isfinite(err_o) & (err_o > 0) & (ew_o / err_o >= MIN_LINE_SNR)
    d['OBIN'] = np.where(oko, np.digitize(ew_o, OXY_EDGES), -1).astype(np.int8)
    d['OXYGEN_EW'] = ew_o
    d['BBIN'] = d.BAL_CLASS.fillna(-1).astype(np.int8)
    json.dump(dict(done=2, total=4, status='running', note='merging cells', started=t0), open(STATUS, 'w'))
    # --- cells and merging (on the training set)
    key_cols = ['ZBIN', 'LBIN', 'EBIN', 'SBIN', 'OBIN', 'BBIN']
    lab = d[key_cols].copy()
    lab.loc[~train, :] = -99                        # non-training rows get labels but do not count
    base = d[key_cols].to_numpy().copy()
    cur = base.copy(); cur[~train] = -99
    def counts(arr):
        df = pd.DataFrame(arr, columns=key_cols); df['i'] = np.arange(len(df))
        c = df.groupby(key_cols)['i'].transform('size').to_numpy(); return c
    n6 = counts(cur)
    stage = np.zeros(len(d), np.int8)             # 0 full cell, 1 BAL collapsed, 2 + oxygen, 3 + slope, 4 pool (z,L,EW), 5 below floor even so
    for s_idx, col in enumerate([5, 4, 3], start=1):       # BBIN, OBIN, SBIN
        low = train & (n6 < a.floor)
        if not low.any():
            break
        cur[low, col] = -1; stage[low] = s_idx
        n6 = counts(cur)
    low = train & (n6 < a.floor)
    stage[low] = np.where(stage[low] < 4, 4, stage[low]); n6 = counts(cur)
    d['STACK_ZBIN'], d['STACK_LBIN'], d['STACK_EBIN'], d['STACK_SBIN'], d['STACK_OBIN'], d['STACK_BBIN'] = [cur[:, k] for k in range(6)]
    d['STACK_STAGE'] = np.where(train, stage, -1).astype(np.int8)
    d['STACK_N'] = np.where(train, n6, 0)
    d['STACK_CELL'] = np.where(train, ['z%d_L%d_E%d_S%d_O%d_B%d' % tuple(r) for r in cur], '')
    d['CELL'] = ['z%d_L%d_E%d_S%d_O%d_B%d' % tuple(r) for r in base]
    # --- stack-cell table
    t = d[train].groupby('STACK_CELL').agg(n=('Z', 'size'), z_med=('Z', 'median'), logL_med=('LOGL_NORM', 'median'), snr_med=('SNR_ALL', 'median'),
                                           ew_med=('PRIMARY_EW', 'median'), alpha_med=('ALPHA_LAM', 'median'), oxy_med=('OXYGEN_EW', 'median'),
                                           stage=('STACK_STAGE', 'max'), ZBIN=('STACK_ZBIN', 'first'), LBIN=('STACK_LBIN', 'first'),
                                           EBIN=('STACK_EBIN', 'first'), SBIN=('STACK_SBIN', 'first'), OBIN=('STACK_OBIN', 'first'), BBIN=('STACK_BBIN', 'first'),
                                           primary=('PRIMARY_LINE', 'first'), oxygen=('OXYGEN_LINE', 'first')).reset_index()
    t['ok'] = t.n >= a.floor
    tag = a.survey
    d.to_parquet(f'{PROD}/qso_cells_{tag}.parquet')
    t.to_csv(f'{PROD}/qso_stack_cells_{tag}.csv', index=False); t.to_csv(f'{REPO}/qso_stack_cells_{tag}.csv', index=False)
    pd.DataFrame(med_table).to_csv(f'{REPO}/qso_primary_ew_medians_{tag}.csv', index=False)
    # --- summary
    out = [f'survey {tag}: rows {len(d):,}; training (QCUT & S/N>5 & L_norm) {train.sum():,}; floor {a.floor}']
    full = pd.DataFrame(base[train], columns=key_cols)
    n_full = full.groupby(key_cols).size()
    out.append(f'six-axis cells occupied: {len(n_full):,}; with >= floor: {(n_full >= a.floor).sum():,} holding {int(n_full[n_full >= a.floor].sum()):,} quasars')
    ok = t[t.ok]
    out.append(f'stack cells after merging (>= floor): {len(ok):,} holding {int(ok.n.sum()):,} quasars ({100*ok.n.sum()/train.sum():.1f} %); below floor even as pool: {len(t)-len(ok):,} cells, {int(t[~t.ok].n.sum()):,} quasars')
    out.append('by stage: ' + ', '.join(f'{k}: {v}' for k, v in {0: 'full', 1: 'BAL collapsed', 2: '+oxygen', 3: '+slope', 4: 'pool'}.items() for v in [int((ok.stage == k).sum())]))
    out.append('axes with fine bins (training): EW %.0f %%, slope %.0f %%, oxygen %.0f %%, BAL known %.0f %%' % tuple(100 * (d.loc[train, c] >= 0).mean() for c in ['EBIN', 'SBIN', 'OBIN', 'BBIN']))
    rows = []
    for zb in range(len(ZEDGES) - 1):
        s = ok[ok.ZBIN == zb]; tt = t[t.ZBIN == zb]
        rows.append(dict(z=f'{ZEDGES[zb]:.1f}-{ZEDGES[zb+1]:.1f}', n_train=int(tt.n.sum()), stacks=len(s), full6=int((s.stage == 0).sum()), in_stacks=int(s.n.sum())))
    out.append(pd.DataFrame(rows).to_string(index=False))
    txt = '\n'.join(out); print(txt); open(f'{PROD}/qso_cells_{tag}.txt', 'w').write(txt)
    json.dump(dict(done=4, total=4, status='done', note=f'{len(ok)} stacks', started=t0), open(STATUS, 'w'))
    print(f'{time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
