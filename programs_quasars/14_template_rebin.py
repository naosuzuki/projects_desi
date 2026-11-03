"""Template amplitudes corrected for the window bias of the faint winnorm templates, and a mixed template set for the re-stack
(2026-11-01, after the template-normalized run 3).

Run 3 normalized every member by the amplitude of its stack's mean spectrum (05_stack.py --v3 --template, templates from the
window-normalized run 2, qso_stacks_winnorm.npz). The faint bins (M_AB(3025) > -20) still dive at 3000 A: they were filled by
quasars at z 0.19-0.3 whose 25-pixel 3000-3050 A window sits at the DESI blue edge and came out noisily low (selection on the
window noise; their template amplitudes are 2-7x brighter). The faint winnorm templates carry the same dive, so their unit mean
at 3000-3050 A is set on the artifact and the run-3 amplitudes of those cells are biased low by the same factor.

This script
  1. corrects the run-3 amplitudes per cell: NORM' = NORM x median over 3100-7000 A of T_winnorm / T_fitchain (f ~ NORM T_win
     ~ NORM' T_fit; the fit-chained run 1 stacks have no window selection, their faint shapes are smooth)
     -> v3/qso_norm_table_template_corr.parquet (cells without a fit-chained template keep the run-3 amplitude);
  2. writes the templates for run 4: the fit-chained stack where it exists, else the run-3 stack
     -> v3/qso_stacks_tplmix.npz (WMEAN, N, STACK_CELL, as 05_stack.load_templates reads them);
  3. prints the dry run of the new luminosity bins (04_cells.py --chain --template-bins --template-table <corr>).
Run 4 (run_v3_rebin.sh): 04_cells.py --chain --template-bins --template-table v3/qso_norm_table_template_corr.parquet,
then run_stack_v3.sh with QSO_TEMPLATE=v3/qso_stacks_tplmix.npz, then after_stack_v3.sh.
"""
import importlib.util, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qsolib as Q


def load(name):
    sp = importlib.util.spec_from_file_location(name, f'{HERE}/{name}.py'); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


def main():
    s5, c4 = load('05_stack'), load('04_cells')
    P = Q.PROD if hasattr(Q, 'PROD') else os.path.expanduser('~/data/sdssdesi_products/quasar')
    cells_path, win, fit = f'{P}/qso_cells_v3_all.parquet', f'{P}/v3/qso_stacks_winnorm.npz', f'{P}/v3/qso_stacks_fitchain.npz'
    run3 = f'{P}/v3/qso_stacks.npz'
    d = pd.read_parquet(cells_path, columns=['SURVEY', 'TARGETID', 'PROGRAM', 'PLATE', 'MJD', 'FIBERID', 'Z', 'LOGFABS_NORM', 'TRAIN', 'LBIN', 'CHAIN_SRC', 'LOGFABS_3025', 'STACK_CELL'])
    cells = sorted(set(d.STACK_CELL[d.TRAIN].dropna()))
    Tw, Tf = s5.load_templates(win, cells), s5.load_templates(fit, cells)
    body = (Q.QWAVE >= 3100) & (Q.QWAVE <= 7000)
    c = np.full(len(cells), np.nan)
    for k in range(len(cells)):
        ok = body & (Tw[k] != 0) & (Tf[k] != 0)
        if ok.sum() >= 500:
            c[k] = np.median(Tw[k][ok] / Tf[k][ok])
    cc = pd.DataFrame({'STACK_CELL': cells, 'C': c, 'L': [int(x.split('_')[0][1:]) for x in cells]})
    print('cells with a correction:', np.isfinite(c).sum(), 'of', len(cells))
    print(cc.groupby('L').C.agg(['size', 'median', 'min', 'max']).round(3).to_string())
    # --- 1. corrected amplitudes (members -> cell through the exact IDs of the measurement table)
    m = pd.read_parquet(c4.NORM_TABLE_V3, columns=['ID', 'SURVEY', 'TARGETID', 'PROGRAM', 'PLATE', 'MJD', 'FIBERID'])
    t = pd.read_parquet(c4.NORM_TABLE_TPL)
    m = m.merge(t[['ID', 'NORM', 'NORM_SRC', 'STORE_SURVEY']], on='ID')
    desi = (m.SURVEY == 'desi').to_numpy(); dt = d[d.TRAIN]
    a = m[desi].merge(dt[dt.SURVEY == 'desi'][['TARGETID', 'PROGRAM', 'STACK_CELL']].drop_duplicates(['TARGETID', 'PROGRAM']), on=['TARGETID', 'PROGRAM'], how='left')
    b = m[~desi].merge(dt[dt.SURVEY != 'desi'][['PLATE', 'MJD', 'FIBERID', 'STACK_CELL']].drop_duplicates(['PLATE', 'MJD', 'FIBERID']), on=['PLATE', 'MJD', 'FIBERID'], how='left')
    m = pd.concat([a, b]).merge(cc[['STACK_CELL', 'C']], on='STACK_CELL', how='left')
    print(f'members with a cell: {m.STACK_CELL.notna().sum():,}, with a correction: {m.C.notna().sum():,} of {len(m):,}')
    m['NORM'] = np.where(np.isfinite(m.C), m.NORM * m.C, m.NORM)
    corr = f'{P}/v3/qso_norm_table_template_corr.parquet'
    m[['ID', 'NORM', 'NORM_SRC', 'STORE_SURVEY']].to_parquet(corr, index=False); print('wrote', corr)
    # --- 2. mixed templates
    zf, z3 = np.load(fit, allow_pickle=True), np.load(run3, allow_pickle=True)
    lab = list(z3['STACK_CELL']); idx = {x: i for i, x in enumerate(zf['STACK_CELL'])}
    wm, n = z3['WMEAN'].copy(), z3['N'].copy(); nf = 0
    for i, x in enumerate(lab):
        if x in idx:
            wm[i], n[i] = zf['WMEAN'][idx[x]], zf['N'][idx[x]]; nf += 1
    mix = f'{P}/v3/qso_stacks_tplmix.npz'
    np.savez(mix, WAVE=z3['WAVE'], WMEAN=wm, N=n, STACK_CELL=np.array(lab)); print(f'wrote {mix}: {nf} fit-chained + {len(lab) - nf} run-3 templates')
    # --- 3. dry run of the bins
    z = d.Z.to_numpy(float)
    meas = c4.measured_ratio(d, tpl=corr)
    logr, lf, src = c4.chain_to_3025(z, d.LOGFABS_NORM.to_numpy(float), meas=meas, meas_all=True)
    lb = np.where(np.isfinite(lf), np.floor(np.nan_to_num(lf) / c4.LBIN_DEX), -99).astype(int)
    tr = d.TRAIN.to_numpy(); old = d.LBIN.to_numpy()[tr]
    dl = lb[tr] - old
    print('LBIN change (new - old), training:', {int(k): int(v) for k, v in zip(*np.unique(np.clip(dl, -3, 3), return_counts=True))})
    print(pd.DataFrame({'old': pd.Series(old).value_counts().sort_index(), 'new': pd.Series(lb[tr]).value_counts().sort_index()}).fillna(0).astype(int).to_string())
    f = pd.DataFrame({'z': z[tr], 'dlf': lf[tr] - d.LOGFABS_3025.to_numpy()[tr], 'old': old}); f = f[f.old <= -3]
    print('old faint bins (LBIN <= -3) by z: median shift in dex'); print(f.groupby(pd.cut(f.z, [0, 0.1, 0.19, 0.3, 0.5, 0.8, 5]), observed=True).dlf.agg(['size', 'median']).round(3).to_string())


if __name__ == '__main__':
    main()
