#!/usr/bin/env python3
"""Fit every training quasar with the released eigenspectra and measure its BAL troughs against the fitted model
(Nao 2026-11-04: 'use the current eigenspectra and identify all of the BALs').

Input: the v2 rest-frame store on exdisk3 (one npz per coadd / plate: normalized rest flux and ivar, float16, cropped
to the coverage; 1,284,634 training quasars), the basis (default data_eigenspectra/quasar_eigenspectra_v3.fits = v6).
Per quasar: weighted least squares of the first NVEC vectors over the good pixels redward of Lya (the forest is not in
the fit); then the pixels that sit below the model (absorption: the BAL velocity windows of C IV, Si IV, Mg II and
Al III where the 5-pixel smoothed flux / model < 0.9 by >= 2 sigma, and any pixel > 4 sigma below the model) are
masked and the fit repeated (NITER passes), so the model is the emission + continuum shape without the troughs.
qsolib.forest_and_bal then measures AI, BI, trough velocities and BAL_CLASS with the model as the continuum
(the pass-Q class used the per-spectrum power law; BAL_CLASS: 0 non-BAL, 1 AI > 0, 2 BI > 0, 3 LoBAL, -1 no coverage).

Store IDs: SDSS exact ((PLATE 1e5 + MJD) 1e4 + FIBERID); DESI TARGETIDs were stored as float64 and are rounded to a
multiple of 8 above 2^53, so a store row is matched to the v3 cells table by (rounded TARGETID, program from the file
name); the 0.7 % with two training quasars at the same rounded id are resolved through the v2 cells table (the store's
CELL_I is the index of the v2 stack-cell label) and the redshift.

Usage (shards as 05_stack.py; results in PROD/v3/eigenfit/, checkpointed per shard every --save files):
  python3 16_eigenfit.py --survey desi --shard 0/3 ; --survey sdss --shard 0/1 ; --merge
Merge -> PROD/v3/qso_eigenfit_<tag>.parquet: ID, SURVEY, PROGRAM, TARGETID/PLATE/MJD/FIBERID, Z, C0..C<n-1>, CHI2, NPIX,
NCLIP, BAL_* of the fitted continuum, BAL_CLASS (eigen) and BAL_CLASS_Q (pass Q)."""
import os, sys, glob, time, argparse, warnings
import numpy as np, pandas as pd, fitsio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qsolib as Q

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar')
STORE = '/Volumes/exdisk3/sdssdesi_qso_restframe_v2'
BASIS = os.path.join(HERE, '..', 'data_eigenspectra', 'quasar_eigenspectra_v3.fits')
LYA = 1215.67
BAL_LINES = ((1549.06, 'CIV'), (1399.8, 'SIIV'), (2798.75, 'MGII'), (1857.4, 'ALIII'))
BAL_COLS = [f'BAL_{k}_{t}' for t in ('CIV', 'SIIV', 'MGII', 'ALIII') for k in ('AI', 'NTROUGH', 'VMIN', 'VMAX')] + \
           [f'BAL_{k}_{t}' for t in ('CIV', 'SIIV') for k in ('BI', 'NBI')] + ['BAL_CLASS', 'FOREST_TRANS', 'FOREST_NPIX', 'REST_LAMMIN', 'REST_LAMMAX']


def load_basis(path, nvec):
    E = fitsio.read(path, ext='EIGVEC')[:nvec].astype(np.float64); ok = fitsio.read(path, ext='PIXMASK') > 0
    E[:, ~ok] = 0.0
    return E, ok


def read_store(fn):
    z = np.load(fn)
    ids, ci, starts, off = z['ID'], z['CELL_I'], z['START'], z['OFFSET']
    n = len(ids); F = np.zeros((n, Q.QNREST)); IV = np.zeros((n, Q.QNREST))
    flux, ivar = z['FLUX'].astype(np.float64), z['IVAR'].astype(np.float64)
    for i in range(n):
        a, b = off[i], off[i + 1]; s = starts[i]; F[i, s:s + b - a] = flux[a:b]; IV[i, s:s + b - a] = ivar[a:b]
    IV[~np.isfinite(F)] = 0; F[~np.isfinite(F)] = 0
    return ids, ci, F, IV


def wls(E, F, W):
    """Weighted least squares of each row of F on the rows of E (E [k, m], F, W [n, m]) -> coefficients [n, k]."""
    M = np.einsum('km,nm,lm->nkl', E, W, E, optimize=True); b = np.einsum('km,nm->nk', E, W * F, optimize=True)
    k = E.shape[0]; tr = np.trace(M, axis1=1, axis2=2); bad = ~(tr > 0)      # no usable pixel (e.g. nothing above Ly-alpha): no fit
    M += 1e-12 * tr[:, None, None] * np.eye(k)[None]; M[bad] = np.eye(k); b[bad] = 0.0
    try:
        c = np.linalg.solve(M, b[..., None])[..., 0]
    except np.linalg.LinAlgError:                                          # a near-degenerate row: solve one by one, least squares
        c = np.stack([np.linalg.lstsq(Mi, bi, rcond=None)[0] for Mi, bi in zip(M, b)])
    c[bad] = np.nan
    return c


def absorption_mask(F, IV, model):
    """True where a pixel sits below the model: BAL velocity windows with the smoothed flux/model < 0.9 by >= 2 sigma,
    and any pixel > 4 sigma below the model."""
    with np.errstate(all='ignore'):
        sig = 1.0 / np.sqrt(np.where(IV > 0, IV, np.nan))
        bad = (IV > 0) & (F < model - 4 * sig)
        k = np.ones(Q.BAL_SMOOTH) / Q.BAL_SMOOTH
        for l0, _ in BAL_LINES:
            v = (Q.QWAVE / l0 - 1) * Q.CKMS; s = (v >= Q.BAL_VRANGE[0]) & (v <= Q.BAL_VRANGE[1])
            fn = np.where(IV[:, s] > 0, F[:, s] / (Q.BAL_DEPTH * model[:, s]), 1.0); en = sig[:, s] / (Q.BAL_DEPTH * model[:, s])
            fs = np.apply_along_axis(lambda r: np.convolve(r, k, mode='same'), 1, fn)
            es = np.sqrt(np.apply_along_axis(lambda r: np.convolve(np.nan_to_num(r, nan=1e9)**2, k**2, mode='same'), 1, en))
            bad[:, s] |= (IV[:, s] > 0) & (fs + Q.BAL_NSIG * es < 1)
    return bad


def fit_file(E, ok, F, IV, niter):
    fitpix = ok & (Q.QWAVE > LYA)
    W = np.where(fitpix[None, :], IV, 0.0); clip = np.zeros(F.shape, bool)
    for it in range(niter + 1):
        if it:
            clip = absorption_mask(F, IV, model); W = np.where(clip, 0.0, np.where(fitpix[None, :], IV, 0.0))
        c = wls(E, F, W); model = c @ E
    with np.errstate(all='ignore'):
        chi2 = np.sum(W * (F - model)**2, axis=1); npix = (W > 0).sum(axis=1)
        cont = np.where((model > 0) & ok[None, :], model, np.nan)
        bal = Q.forest_and_bal(F, IV, cont)
    return c, chi2, npix, clip.sum(axis=1), bal


def desi_lookup():
    """(rounded TARGETID, program) -> exact TARGETID, Z; ambiguous ones -> list of (CELL_I of the v2 cell, Z, TARGETID)."""
    c = pd.read_parquet(f'{PROD}/qso_cells_v3_all.parquet', columns=['SURVEY', 'TARGETID', 'PROGRAM', 'Z', 'TRAIN', 'STACK_CELL', 'BAL_CLASS'])
    c = c[(c.SURVEY == 'desi') & c.TRAIN & (c.STACK_CELL != '')].copy()
    c['R'] = c.TARGETID.astype(np.int64).astype(np.float64).astype(np.int64)
    g = c.groupby(['R', 'PROGRAM']).size(); amb = g[g > 1].index
    v2 = pd.read_parquet(f'{PROD}/qso_cells_all.parquet', columns=['SURVEY', 'TARGETID', 'PROGRAM', 'Z', 'TRAIN', 'STACK_CELL'])
    v2 = v2[v2.TRAIN & (v2.STACK_CELL != '')].copy()
    cells = np.array(sorted(set(v2.STACK_CELL)))             # the store's CELL_I indexes the sorted cells of ALL surveys (2575)
    v2['CELL_I'] = np.searchsorted(cells, v2.STACK_CELL.to_numpy()).astype(np.int32); v2 = v2[v2.SURVEY == 'desi']
    v2['R'] = v2.TARGETID.astype(np.int64)
    amb_v2 = v2.set_index(['R', 'PROGRAM']).loc[v2.set_index(['R', 'PROGRAM']).index.isin(amb)].reset_index()
    amb_v3 = c.set_index(['R', 'PROGRAM']).loc[c.set_index(['R', 'PROGRAM']).index.isin(amb)].reset_index()
    uniq = c[~c.set_index(['R', 'PROGRAM']).index.isin(amb)].set_index(['R', 'PROGRAM'])[['TARGETID', 'Z', 'BAL_CLASS']].astype({'TARGETID': 'Int64'})
    return uniq, amb_v2, amb_v3


def match_desi(ids, ci, program, uniq, amb_v2, amb_v3):
    r = ids.astype(np.int64); key = pd.MultiIndex.from_arrays([r, np.full(len(r), program)])
    hit = uniq.reindex(key); tid = hit.TARGETID.array; z = hit.Z.to_numpy(dtype=float); bq = hit.BAL_CLASS.to_numpy()   # TARGETID kept Int64 (a float64 pass would round the ids)
    for i in np.where(pd.isna(tid))[0]:                       # ambiguous rounded id: v2 cell index -> Z -> v3 row
        a2 = amb_v2[(amb_v2.R == r[i]) & (amb_v2.PROGRAM == program) & (amb_v2.CELL_I == ci[i])]
        a3 = amb_v3[(amb_v3.R == r[i]) & (amb_v3.PROGRAM == program)]
        if len(a2) == 1 and len(a3):
            j = np.argmin(np.abs(a3.Z.to_numpy(float) - float(a2.Z.iloc[0])))
            if abs(a3.Z.iloc[j] - a2.Z.iloc[0]) < 1e-3:
                tid[i] = a3.TARGETID.iloc[j]; z[i] = a3.Z.iloc[j]; bq[i] = a3.BAL_CLASS.iloc[j]
    return tid, z, bq


def sdss_lookup():
    c = pd.read_parquet(f'{PROD}/qso_cells_v3_all.parquet', columns=['SURVEY', 'PLATE', 'MJD', 'FIBERID', 'Z', 'TRAIN', 'STACK_CELL', 'BAL_CLASS'])
    c = c[(c.SURVEY != 'desi') & c.TRAIN & (c.STACK_CELL != '')].copy()
    c['ID'] = (c.PLATE.astype(np.int64) * 100000 + c.MJD.astype(np.int64)) * 10000 + c.FIBERID.astype(np.int64)
    return c.drop_duplicates('ID').set_index('ID')[['PLATE', 'MJD', 'FIBERID', 'Z', 'SURVEY', 'BAL_CLASS']].rename(columns={'BAL_CLASS': 'BAL_CLASS_Q'})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='desi', choices=['desi', 'sdss'])
    ap.add_argument('--shard', default='0/1'); ap.add_argument('--nfiles', type=int, default=0)
    ap.add_argument('--basis', default=BASIS); ap.add_argument('--nvec', type=int, default=10); ap.add_argument('--niter', type=int, default=2)
    ap.add_argument('--tag', default='v6'); ap.add_argument('--save', type=int, default=300)
    ap.add_argument('--merge', action='store_true'); ap.add_argument('--pattern', default='*.npz', help='test: store files to use')
    a = ap.parse_args()
    out = f'{PROD}/v3/eigenfit'; os.makedirs(out, exist_ok=True)
    if a.merge:
        return merge(out, a.tag)
    si, sn = map(int, a.shard.split('/'))
    E, ok = load_basis(a.basis, a.nvec)
    if a.survey == 'desi':
        files = sorted(glob.glob(f'{STORE}/desi/{a.pattern}')); uniq, amb_v2, amb_v3 = desi_lookup()
    else:
        files = sorted(glob.glob(f'{STORE}/dr17/{a.pattern}')) + sorted(glob.glob(f'{STORE}/dr8/{a.pattern}')); look = sdss_lookup()
    files = [f for i, f in enumerate(files) if i % sn == si]
    if a.nfiles: files = files[:a.nfiles]
    part = f'{out}/eigenfit_{a.tag}_{a.survey}_shard{si}of{sn}.parquet'; done = set(); rows = []
    if os.path.exists(part):
        prev = pd.read_parquet(part); done = set(prev.FILE); rows = [prev]
    todo = [f for f in files if os.path.basename(f) not in done]
    print(f'{a.survey} shard {si}/{sn}: {len(files)} files, {len(todo)} to do; basis {os.path.basename(a.basis)} x {a.nvec}', flush=True)
    t0 = time.time(); nq = 0; nun = 0
    for k, fn in enumerate(todo):
        base = os.path.basename(fn)
        try:
            ids, ci, F, IV = read_store(fn)
        except Exception as e:
            print('READ ERROR', base, repr(e), flush=True); continue
        if len(ids) == 0: continue
        c, chi2, npix, nclip, bal = fit_file(E, ok, F, IV, a.niter)
        d = pd.DataFrame({f'C{j}': c[:, j] for j in range(a.nvec)})
        d['CHI2'] = chi2; d['NPIX'] = npix.astype(np.int32); d['NCLIP'] = nclip.astype(np.int32)
        for col in BAL_COLS: d[col] = bal[col]
        d = d.rename(columns={'BAL_CLASS': 'BAL_CLASS_E'})
        if a.survey == 'desi':
            program = base.split('-')[2]
            tid, z, bq = match_desi(ids, ci, program, uniq, amb_v2, amb_v3)
            d.insert(0, 'TARGETID', tid); d.insert(1, 'PROGRAM', program); d.insert(2, 'Z', z); d.insert(3, 'BAL_CLASS_Q', bq); d.insert(0, 'SURVEY', 'desi'); d['RID'] = ids.astype(np.int64); d['CELL_I'] = ci.astype(np.int32)
            nun += int(pd.isna(tid).sum())
        else:
            key = ids.astype(np.int64); hit = look.reindex(key)
            d.insert(0, 'ID', key)
            for col in ['PLATE', 'MJD', 'FIBERID', 'Z', 'SURVEY', 'BAL_CLASS_Q']: d[col] = hit[col].to_numpy()
            nun += int(hit.Z.isna().sum())
        d['FILE'] = base; rows.append(d); nq += len(d)
        if (k + 1) % a.save == 0 or k == len(todo) - 1:
            pd.concat(rows, ignore_index=True).to_parquet(part + '.tmp', index=False); os.replace(part + '.tmp', part); rows = [pd.read_parquet(part)]
            dt = time.time() - t0
            print(f'{k + 1}/{len(todo)} files, {nq:,} quasars ({nun} unmatched), {dt / 60:.1f} min, {1e3 * dt / max(nq, 1):.2f} ms/quasar', flush=True)
    print('shard finished', flush=True)


def merge(out, tag):
    parts = sorted(glob.glob(f'{out}/eigenfit_{tag}_*_shard*.parquet'))
    d = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    # re-match the DESI ids from the store's rounded id, program and v2 cell (the shard files written before 2026-11-04 10:00
    # carried rounded TARGETIDs for 17 % of the rows: a float64 pass through reindex); the fit itself does not depend on the id
    desi = (d.SURVEY == 'desi').to_numpy(); uniq, amb_v2, amb_v3 = desi_lookup(); tid = pd.array(np.zeros(len(d)), dtype='Int64'); zz = np.full(len(d), np.nan); bq = np.full(len(d), np.nan)
    for prog in d.PROGRAM[desi].dropna().unique():
        r = np.where(desi & (d.PROGRAM == prog).to_numpy())[0]
        t, z, b = match_desi(d.RID.to_numpy()[r], d.CELL_I.to_numpy()[r], prog, uniq, amb_v2, amb_v3); tid[r] = t; zz[r] = z; bq[r] = pd.array(b, dtype='float')
    changed = int((tid[desi] != d.TARGETID.array[desi]).sum()); d.loc[desi, 'TARGETID'] = tid[desi]; d.loc[desi, 'Z'] = zz[desi]; d.loc[desi, 'BAL_CLASS_Q'] = bq[desi]
    print(f'DESI ids re-matched: {changed:,} changed')
    print(f'{len(d):,} rows from {len(parts)} shards; unmatched {int(d.Z.isna().sum())}')
    d['BAL_CLASS_E'] = d.BAL_CLASS_E.astype(np.int8)
    o = f'{PROD}/v3/qso_eigenfit_{tag}.parquet'; d.to_parquet(o, index=False); print('wrote', o)
    print(pd.crosstab(d.BAL_CLASS_Q.fillna(-9), d.BAL_CLASS_E, rownames=['pass Q'], colnames=['eigen']))


if __name__ == '__main__':
    warnings.filterwarnings('ignore'); main()
