#!/usr/bin/env python3
"""Re-partition the training quasars into new stack cells from the v2 rest-frame store, without re-reading the coadds
(2026-11-04, for the BAL / non-BAL split: Nao 'two groups of eigenspectra').

The store (exdisk3, written by the v2 walk 2026-10-30) holds every training quasar's rest-frame spectrum on the 1e-4 grid,
normalized by its v2 window mean (NORM in the file), float16, without the telluric correction. Run 5 normalized each member
by its template amplitude NORM_5 (v3/qso_norm_table_template.parquet, exact IDs) and divided out the telluric factor
(05_stack.telluric_correct). Here each member is restored to the run-5 normalization, f = F_store x NORM_2 / NORM_5
(ivar x (NORM_5 / NORM_2)^2), the telluric factor is divided out with the member's own eigenspectrum fit as the model
(16_eigenfit.py: C0..C9 of the released basis, in the store's units), and the member is added to its cell of the new cells
table(s). Same accumulators and merge as 05_stack.py (WMEAN, WERR, N, ...), so 06_empca.py and 12_lum_stacks.py read the result.
Members without a run-5 amplitude (407 window fallbacks) or without a fit are skipped. The v2 walk wrote 1.6 % of the DESI
rows twice or more (the same TARGETID with a different spectrum, often a junk one 25x fainter); run 5 stacked one spectrum per
quasar, so per quasar the row whose window mean x chain ratio is closest to the run-5 amplitude is kept, and any row off by
more than 0.3 dex (on SDSS, without duplicates, 0.03 % are beyond 0.3) is a different spectrum and is dropped.

  python3 17_restack_store.py --cells qso_cells_v3bal_all --survey desi --shard 0/3 [--cells ... more tables]
  python3 17_restack_store.py --cells qso_cells_v3bal_all --merge
Accumulators in PROD/<cells>/acc/, result PROD/<cells>/qso_stacks.npz + qso_stacks_cells.csv (+ cells_index.csv).
With --cells qso_cells_v3_all (the run-5 cells) the result must reproduce v3/qso_stacks.npz (the store is float16: ~1e-3
per member, ~1e-4 per stack) - the validation of this path."""
import os, sys, glob, time, argparse, warnings, importlib.util
import numpy as np, pandas as pd, fitsio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qsolib as Q


def load(name):
    sp = importlib.util.spec_from_file_location(name, f'{HERE}/{name}.py'); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


S5, EF = load('05_stack'), load('16_eigenfit')
PROD, STORE = EF.PROD, EF.STORE


def cell_index(table):
    """Training quasars with a stack cell in <table> -> (sorted cell labels, DataFrame keyed by exact id with CELL_I)."""
    c = pd.read_parquet(f'{PROD}/{table}.parquet', columns=['SURVEY', 'TARGETID', 'PROGRAM', 'PLATE', 'MJD', 'FIBERID', 'TRAIN', 'STACK_CELL', 'CHAIN_LOGR'])
    c = c[c.TRAIN & (c.STACK_CELL != '')].copy()
    cells = np.array(sorted(c.STACK_CELL.unique())); c['CELL_I'] = np.searchsorted(cells, c.STACK_CELL.to_numpy()).astype(np.int32)
    desi = c.SURVEY == 'desi'
    c['ID'] = np.where(desi, c.TARGETID.fillna(0).astype(np.int64),
                       (c.PLATE.fillna(0).astype(np.int64) * 100000 + c.MJD.fillna(0).astype(np.int64)) * 10000 + c.FIBERID.fillna(0).astype(np.int64))
    c['KEY'] = np.where(desi, c.ID.astype(str) + '_' + c.PROGRAM.fillna('').astype(str), c.ID.astype(str))
    return cells, c.drop_duplicates('KEY').set_index('KEY')[['CELL_I', 'STACK_CELL', 'CHAIN_LOGR']]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', action='append', required=True, help='cells table(s) in PROD (name without .parquet); one accumulator set each')
    ap.add_argument('--survey', default='desi', choices=['desi', 'sdss']); ap.add_argument('--shard', default='0/1')
    ap.add_argument('--nfiles', type=int, default=0); ap.add_argument('--pattern', default='*.npz')
    ap.add_argument('--norm', default=f'{PROD}/v3/qso_norm_table_template.parquet', help='run-5 member amplitudes (exact IDs)')
    ap.add_argument('--fit', default=f'{PROD}/v3/qso_eigenfit_v6.parquet', help='16_eigenfit.py table (model for the telluric correction)')
    ap.add_argument('--basis', default=EF.BASIS); ap.add_argument('--no-telluric', action='store_true')
    ap.add_argument('--merge', action='store_true'); ap.add_argument('--save', type=int, default=200)
    ap.add_argument('--lrmax', type=float, default=0.3, help='max |log10(window norm x chain / run-5 amplitude)| of a kept row')
    a = ap.parse_args()
    if a.merge:
        for t in a.cells: merge(t)
        return
    si, sn = map(int, a.shard.split('/'))
    sets = []
    for t in a.cells:
        cells, idx = cell_index(t); accdir = f'{PROD}/{t}/acc'; os.makedirs(accdir, exist_ok=True)
        pd.DataFrame({'CELL_I': np.arange(len(cells)), 'STACK_CELL': cells}).to_csv(f'{accdir}/cells_index.csv', index=False)
        acc = S5.Acc(f'{accdir}/acc_{a.survey}_shard{si}of{sn}.npz', len(cells)); sets.append((t, idx, acc))
        print(f'{t}: {len(cells)} cells, {len(idx):,} members; {len(acc.done)} files already done', flush=True)
    norm = pd.read_parquet(a.norm).drop_duplicates('ID').set_index('ID')
    fit = pd.read_parquet(a.fit)
    nv = sum(c.startswith('C') and c[1:].isdigit() for c in fit.columns); ccols = [f'C{j}' for j in range(nv)]
    if a.survey == 'desi':
        fit = fit[fit.SURVEY == 'desi']; fit = fit[fit.TARGETID.notna()]
        fit['KEY'] = fit.TARGETID.astype(np.int64).astype(str) + '_' + fit.PROGRAM.astype(str)
        uniq, amb_v2, amb_v3 = EF.desi_lookup(); files = sorted(glob.glob(f'{STORE}/desi/{a.pattern}'))
    else:
        fit = fit[fit.SURVEY != 'desi']; fit['KEY'] = fit.ID.astype(np.int64).astype(str); look = EF.sdss_lookup()
        files = sorted(glob.glob(f'{STORE}/dr17/{a.pattern}')) + sorted(glob.glob(f'{STORE}/dr8/{a.pattern}'))
    fit = fit.drop_duplicates('KEY').set_index('KEY')[ccols + ['Z']]
    E, ok = EF.load_basis(a.basis, nv)
    if not a.no_telluric:
        sys.path.insert(0, os.path.join(HERE, '..', 'programs_telluric')); import tellapply; S5.TELL = tellapply
    files = [f for i, f in enumerate(files) if i % sn == si]
    if a.nfiles: files = files[:a.nfiles]
    done_all = set.intersection(*[s[2].done for s in sets]); todo = [f for f in files if os.path.basename(f) not in done_all]
    print(f'{a.survey} shard {si}/{sn}: {len(files)} files, {len(todo)} to do', flush=True)
    t0 = time.time(); nq = nskip = nunm = 0
    for k, fn in enumerate(todo):
        base = os.path.basename(fn)
        try:
            ids, ci, F, IV = EF.read_store(fn); z0 = np.load(fn); n2 = z0['NORM']; del z0
        except Exception as e:
            print('READ ERROR', base, repr(e), flush=True); continue
        if len(ids) == 0:
            for _, _, acc in sets: acc.done.add(base)
            continue
        if a.survey == 'desi':
            program = base.split('-')[2]; tid, z, _ = EF.match_desi(ids, ci, program, uniq, amb_v2, amb_v3)
            key = np.array(['' if pd.isna(t) else f'{int(t)}_{program}' for t in tid]); nid = np.array([0 if pd.isna(t) else int(t) for t in tid], np.int64)
            survey = np.full(len(ids), 'desi')
        else:
            nid = ids.astype(np.int64); key = nid.astype(str); hit = look.reindex(nid); z = hit.Z.to_numpy(float)
            survey = hit.SURVEY.fillna('').to_numpy()
        n5 = norm.NORM.reindex(nid).to_numpy(float); src5 = norm.NORM_SRC.reindex(nid).to_numpy()
        ft = fit.reindex(key); cf = ft[ccols].to_numpy(float)
        R = 10.0**sets[0][1].CHAIN_LOGR.reindex(key).to_numpy(float)
        with np.errstate(all='ignore'): lr = np.abs(np.log10(n2 * R / n5))
        use = (key != '') & np.isfinite(n5) & (src5 == 0) & np.isfinite(cf).all(axis=1) & np.isfinite(z) & (n2 > 0) & (lr <= a.lrmax)
        best = pd.DataFrame({'key': key, 'lr': np.where(use, lr, np.inf)}).groupby('key').lr.transform('idxmin').to_numpy() == np.arange(len(key))   # one row per quasar
        use &= best
        nunm += int((key == '').sum()); nskip += int((~use).sum())
        if not use.any():
            for _, _, acc in sets: acc.done.add(base)
            continue
        F, IV, cf, z, n2, n5, key, survey = F[use], IV[use], cf[use], z[use], n2[use], n5[use], key[use], survey[use]
        sc = (n2 / n5)[:, None]; f = F * sc; iv = IV / sc**2
        if not a.no_telluric:
            model = (cf @ E) * sc; model[:, ~ok] = 0.0
            for sv in np.unique(survey):
                r = survey == sv; f[r], iv[r], nt = S5.telluric_correct(f[r], iv[r], z[r], model[r], sv)
        for t, idx, acc in sets:
            m = idx.CELL_I.reindex(key).to_numpy()
            has = np.isfinite(m.astype(float))
            if has.any(): acc.add(m[has].astype(np.int32), f[has], iv[has])
            acc.nobj += int(has.sum()); acc.done.add(base)
        nq += int(use.sum())
        if (k + 1) % a.save == 0 or k == len(todo) - 1:
            for _, _, acc in sets: acc.save()
            dt = time.time() - t0
            print(f'{k + 1}/{len(todo)} files, {nq:,} members ({nskip} skipped, {nunm} unmatched), {dt / 60:.1f} min, {1e3 * dt / max(nq, 1):.2f} ms/member', flush=True)
    print('shard finished', flush=True)


def merge(table):
    accdir = f'{PROD}/{table}/acc'; out = f'{PROD}/{table}'
    cells = pd.read_csv(f'{accdir}/cells_index.csv'); tot = None; nobj = 0; parts = []
    for p in sorted(glob.glob(f'{accdir}/acc_*_shard*.npz')):
        z = np.load(p); parts.append(os.path.basename(p)); nobj += int(z['nobj'])
        if tot is None: tot = {k: z[k].astype(np.float64) for k in ('s1', 's2', 'n', 'wf', 'w')}
        else:
            for k in tot: tot[k] += z[k]
    n = tot['n']
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, tot['s1'] / n, np.nan); sd = np.sqrt(np.maximum(np.where(n > 1, tot['s2'] / n - mean**2, np.nan), 0))
        wmean = np.where(tot['w'] > 0, tot['wf'] / tot['w'], np.nan); werr = np.where(tot['w'] > 0, 1 / np.sqrt(tot['w']), np.nan)
    np.savez(f'{out}/qso_stacks.npz', WAVE=Q.QWAVE, MEAN=mean.astype(np.float32), SD=sd.astype(np.float32), N=n.astype(np.int32),
             WMEAN=wmean.astype(np.float32), WERR=werr.astype(np.float32), STACK_CELL=cells.STACK_CELL.to_numpy())
    ip = f'{PROD}/qso_stack_cells_{table.replace("qso_cells_", "")}.csv'
    info = pd.read_csv(ip) if os.path.exists(ip) else None
    cells['N_MAX'] = n.max(axis=1)
    cells = cells.merge(info, on='STACK_CELL', how='left') if info is not None else cells
    cells.to_csv(f'{out}/qso_stacks_cells.csv', index=False)
    print(f'{table}: merged {parts}: {nobj:,} members, {len(cells)} cells -> {out}/qso_stacks.npz')


if __name__ == '__main__':
    warnings.filterwarnings('ignore'); main()
