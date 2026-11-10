#!/usr/bin/env python3
"""Step 8b = Calibration 2 refinement (2026-10-29): the smooth DESI response beyond scale and tilt, per fiber and per petal
and period, from the residuals of the step-6 stars against their own stacks.

Step 8 (step8_response.py) solved A (mag at 6000 A) and B (mag/dex) per (EXPID, PETAL_LOC) from the photometric joint fit:
a straight line in log lambda. What a line cannot carry (curvature, the camera ends, fiber-to-fiber shape) shows up when a
star, calibrated exactly as in step 6 (A_FLUX, TILT, sub-pixel RV, dereddening, 10 pc), is divided by the stack of its cell:

  dm_s(lambda) = -2.5 log10(star / stack)  ->  25-pixel bins, strong and variable features masked (Balmer, Ca H&K, Na D,
  telluric bands) ->  Legendre coefficients c_0..c_5 in u = 2 (log lambda - log 3600) / (log 9800 - log 3600) - 1,
  errors from the weighted fit scaled by its chi2 (so they include line mismatch, not only photon noise).

Model, per component k, with the coadd's exposure weights w_e (EXPTIME x fraction of cameras in the coadd, normalized):

  c_sk = beta_k . (dCOL0, dMG0, dFEH, FEH_unknown)_s + sum_e w_e (F_k[FIBER_e] + P_k[PETAL_e, period_e]) + eps

  F = static per fiber (5,000; prior tau_F), P = per petal and half-year period (prior tau_P), beta = the intrinsic spread
  inside a cell (colour, luminosity, [M/H] against the cell median). Weighted sparse LSQ (lsqr), sigma^2 = err^2 + sig_int^2,
  4 sigma clipping, as in step 8. Hold-out: fit on even TARGETIDs, predict odd ones.

What this can NOT see: a response error common to every fiber is in the stacks themselves (each stack mixes hundreds of
fibers, petals and nights), so star / stack only measures the differences between fibers, petals and periods. The absolute
shape needs an outside reference (SDSS stacks of the same cells, Gaia XP, CALSPEC). Accordingly the published terms are
relative: F has mean zero inside each petal; P (k >= 2) has mean zero over petals and periods; c_0 and c_1 of P are not
published because A and B of step 8 already carry scale and tilt per exposure and petal.

Convention: CORR2(lambda) = 10^(-0.4 sum_k C_k P_k(u)), multiplied on top of step 8's CORR(lambda).

Outputs (~/data/sdssdesi_products/stars/):
  shape_star_residuals_desi.parquet   per star: TARGETID, PROG, STACK, C0..C5, E0..E5, RMS, NBIN, covariates
  response_shape_fiber.parquet        FIBER, PETAL_LOC, C0..C5, E0..E5, N_STARS (fiber terms, mean zero per petal)
  response_shape_petal_period.parquet PETAL_LOC, PERIOD, C2..C5, E2..E5, N_STARS (PERIOD = half-years from 2021-01-01)
  response_shape_coadd_desi/          per target of the exposure tables: TARGETID, PROGRAM, C0..C5 (exposure-weighted
                                      fiber + petal-period terms), NEXP, NEXP_SOLVED; one parquet per 2,000 healpix files
Usage: python3 step8b_shape.py [--residuals-only] [--from-residuals] [--workers 6]
"""
import os, sys, glob, time, json, argparse
import numpy as np, pandas as pd, fitsio
from numpy.polynomial import legendre
from scipy import sparse
from scipy.sparse.linalg import lsqr
import pyarrow as pa, pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
_spec = importlib.util.spec_from_file_location('step6', os.path.join(HERE, 'step6_stack.py')); S6 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(S6)
import star_hook

PROD = os.path.expanduser('~/data/sdssdesi_products/stars')
CACHE = S6.CACHE
EXPFIB = os.path.expanduser('~/data/sdssdesi_expfibermap')
TAG = 'desi_0.5x0.1_gspphot'
NK = 6; BIN = 25
LMIN, LMAX = np.log10(3600.), np.log10(9800.)
TAU = {'F': 0.05, 'P': 0.05}          # mag, priors toward 0 for fiber and petal-period terms
MJD0 = 59215.0                        # 2021-01-01; PERIOD = int((MJD - MJD0) / 182.6)
MASK_REST = [(3925, 3945), (3960, 3980), (4095, 4110), (4330, 4350), (4850, 4872), (5885, 5900), (6550, 6578)]
MASK_OBS = [(6860, 6960), (7160, 7340), (7590, 7700), (8130, 8350), (8940, 9800)]   # telluric (applied on the rest grid: z < 0.4 %)


def progress(done, total, status, note):
    try:
        json.dump({'done': done, 'total': total, 'status': status, 'note': note, 'started': T0},
                  open(os.path.expanduser('~/data/sdssdesi_status/step_star_response.json'), 'w'))
    except Exception:
        pass


# ---------- part 1: per-star residual curves against the stacks
def _setup():
    global STK, SIDX, CELLS, BINS, U, V
    F = fitsio.FITS(os.path.join(S6.OUTDIR, f'stacks_{TAG}.fits'))
    STK = F['FLUX'].read().astype(float); arr = F['CELLS'].read()
    CELLS = pd.DataFrame({c: (np.char.strip(arr[c].astype(str)) if arr[c].dtype.kind in 'SU' else arr[c].astype(arr[c].dtype.newbyteorder('='))) for c in arr.dtype.names})
    SIDX = {s: i for i, s in enumerate(CELLS.STACK)}
    lam = S6.SWAVE; ok = (np.log10(lam) >= LMIN) & (np.log10(lam) <= LMAX)
    for a, b in MASK_REST + MASK_OBS:
        ok &= ~((lam >= a) & (lam <= b))
    idx = np.flatnonzero(ok)
    BINS = [idx[(idx >= i) & (idx < i + BIN)] for i in range(0, S6.NPIX, BIN)]
    BINS = [b for b in BINS if len(b) >= BIN // 2]
    U = np.array([2 * (np.log10(lam[b]).mean() - LMIN) / (LMAX - LMIN) - 1 for b in BINS])
    V = legendre.legvander(U, NK - 1)


def _fit_star(s, w, stack):
    """dm in bins -> Legendre coefficients, errors scaled by chi2, rms, number of bins."""
    good = (w > 0) & (s > 0) & (stack > 0)
    y = np.full(len(BINS), np.nan); wt = np.zeros(len(BINS))
    with np.errstate(all='ignore'):
        dm = -2.5 * np.log10(s / stack)
    for j, b in enumerate(BINS):
        g = b[good[b]]
        if len(g) >= BIN // 3:
            ww = w[g]; y[j] = np.sum(ww * dm[g]) / ww.sum(); wt[j] = ww.sum()
    m = np.isfinite(y) & (wt > 0)
    if m.sum() < 0.6 * len(BINS):
        return None
    for it in range(2):
        A = V[m] * wt[m, None]; c, *_ = np.linalg.lstsq(A.T @ V[m], A.T @ y[m], rcond=None)
        r = y - V @ c; chi2 = np.sum(wt[m] * r[m]**2) / max(m.sum() - NK, 1)
        z = np.abs(r) * np.sqrt(wt / chi2)
        m2 = m & (z < 4)
        if m2.sum() == m.sum():
            break
        m = m2
    cov = np.linalg.inv((V[m] * wt[m, None]).T @ V[m]) * chi2
    return c, np.sqrt(np.diag(cov)), float(np.sqrt(np.mean(r[m]**2))), int(m.sum())


def _work(task):
    files, df = task
    out = []
    for base in files:
        p = os.path.join(star_hook.STORE, base.replace('.fits', '') + '.npz')
        if not os.path.exists(p):
            continue
        z = np.load(p); pos = {int(t): i for i, t in enumerate(z['TARGETID'])}; rid0 = int(z['RID0'])
        for row in df[df.FILE == base].itertuples(index=False):
            i = pos.get(int(row.TARGETID)); k = SIDX.get(row.STACK)
            if i is None or k is None:
                continue
            f = z['FLUX'][i].astype(float); iv = z['IVAR'][i].astype(float)
            s, w = S6.star_spectrum(row, f, iv, (iv <= 0).astype(np.int8), rid0)
            r = _fit_star(s, w, STK[k])
            if r is None:
                continue
            c, e, rms, nb = r
            out.append((row.TARGETID, row.PROG, row.STACK, *c, *e, rms, nb, row.SNR_ALL, row.DCOL, row.DMG, row.DFEH, row.FEHNA))
    return out


def residuals(workers):
    sel = pd.read_parquet(os.path.join(CACHE, f'step6_selection_{TAG}.parquet'))
    m3 = pd.read_parquet(star_hook.SEL300, columns=['FILE', 'TARGETID'])
    sel = sel.merge(m3, on=['FILE', 'TARGETID'])
    mem = pd.read_parquet(os.path.join(CACHE, f'members_{TAG}.parquet'), columns=['SID', 'MG0', 'COL0', 'FEH']).rename(columns={'FEH': 'FEH_S'})
    sel = sel.merge(mem, on='SID', how='left')
    _setup()
    c = CELLS.set_index('STACK')
    sel['DCOL'] = sel.COL0 - sel.STACK.map(c.COL0_MED); sel['DMG'] = sel.MG0 - sel.STACK.map(c.MG0_MED)
    sel['FEHNA'] = (~np.isfinite(sel.FEH_S)).astype(float)
    sel['DFEH'] = np.where(np.isfinite(sel.FEH_S), sel.FEH_S - sel.STACK.map(c.FEH_MED).fillna(sel.FEH_S), 0.0)
    sel['DFEH'] = sel.DFEH.fillna(0.0); sel[['DCOL', 'DMG']] = sel[['DCOL', 'DMG']].fillna(0.0)
    sel['PROG'] = (sel.PROGRAM.astype(str) == 'dark').astype(np.int8)
    print(f'residuals: {len(sel):,} stars in {sel.FILE.nunique():,} coadds, {sel.STACK.nunique()} stacks, {len(BINS)} bins; {time.time()-T0:.0f}s', flush=True)
    files = sorted(sel.FILE.unique()); chunks = [files[i:i + 200] for i in range(0, len(files), 200)]
    tasks = [(ch, sel[sel.FILE.isin(ch)]) for ch in chunks]
    import multiprocessing as mp
    rows = []
    with mp.get_context('fork').Pool(workers) as pool:
        for n, r in enumerate(pool.imap_unordered(_work, tasks)):
            rows += r
            if n % 10 == 0:
                progress(n, len(tasks) * 2, 'running', f'Calibration 2 refinement: star / stack residuals, {n*200:,} / {len(files):,} coadds')
                print(f'  {min((n+1)*200, len(files)):,}/{len(files):,} coadds, {len(rows):,} stars; {time.time()-T0:.0f}s', flush=True)
    cols = ['TARGETID', 'PROG', 'STACK'] + [f'C{k}' for k in range(NK)] + [f'E{k}' for k in range(NK)] + ['RMS', 'NBIN', 'SNR_ALL', 'DCOL', 'DMG', 'DFEH', 'FEHNA']
    res = pd.DataFrame(rows, columns=cols)
    res.to_parquet(os.path.join(PROD, 'shape_star_residuals_desi.parquet'))
    q = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
    print(f'residuals done: {len(res):,} stars; rms of the binned fit median {res.RMS.median():.4f} mag; coefficient spread (MAD) '
          + ' '.join(f'C{k} {q(res[f"C{k}"]):.4f}' for k in range(NK)) + f'; {time.time()-T0:.0f}s', flush=True)
    return res


# ---------- part 2: fiber and petal-period terms
def exposures(keys=None):
    """exposure rows (TARGETID, PROG, FIBER, PETAL, PERIOD, W) from the expfibermap tables; keys: restrict to these targets."""
    cols = ['TARGETID', 'PROGRAM', 'EXPID', 'MJD', 'PETAL_LOC', 'FIBER', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z']
    files = sorted(glob.glob(os.path.join(EXPFIB, '*', '*.parquet')))
    tid = None if keys is None else set(keys)
    parts = []
    for i, f in enumerate(files):
        d = pd.read_parquet(f, columns=cols)
        if tid is not None:
            d = d[d.TARGETID.isin(tid)]
            if not len(d):
                continue
        parts.append(_exp_rows(d))
        if i % 5000 == 0:
            print(f'  exposure tables {i:,}/{len(files):,}; {time.time()-T0:.0f}s', flush=True)
    e = pd.concat(parts, ignore_index=True)
    e['W'] = e.W / e.groupby(['TARGETID', 'PROG']).W.transform('sum')
    return e


def _exp_rows(d):
    frac = (d.IN_COADD_B.astype(np.float32) + d.IN_COADD_R.astype(np.float32) + d.IN_COADD_Z.astype(np.float32)) / 3
    e = pd.DataFrame({'TARGETID': d.TARGETID.to_numpy(np.int64), 'PROG': (d.PROGRAM.astype(str) == 'dark').to_numpy(np.int8),
                      'FIBER': d.FIBER.to_numpy(np.int32), 'PETAL': d.PETAL_LOC.to_numpy(np.int16),
                      'PERIOD': np.clip(((d.MJD.to_numpy(float) - MJD0) / 182.6).astype(np.int16), 0, 99),
                      'W': d.EXPTIME.to_numpy(np.float32) * frac.to_numpy()})
    return e[(e.W > 0) & (e.FIBER >= 0) & (e.FIBER < 5000)]


def fit_k(res, e, k, mask):
    """one coefficient: y = beta . X + sum_e w (F[fiber] + P[petal, period])."""
    s = res[mask].reset_index(drop=True); s['S'] = np.arange(len(s))
    j = e.merge(s[['TARGETID', 'PROG', 'S']], on=['TARGETID', 'PROG'], how='inner')
    pp = j.PETAL.astype(np.int64) * 100 + j.PERIOD.astype(np.int64)
    pcode, puniq = pd.factorize(pp)
    nF = 5000; nP = len(puniq); nB = 4
    X = s[['DCOL', 'DMG', 'DFEH', 'FEHNA']].to_numpy(float)
    y = s[f'C{k}'].to_numpy(float); err = s[f'E{k}'].to_numpy(float)
    rows = j.S.to_numpy(); w = j.W.to_numpy(float); fcol = j.FIBER.to_numpy(); pcol = nF + pcode
    keep = np.ones(len(y), bool); sig_int = 0.5 * float(np.std(y)); y0 = float(np.median(y))
    for it in range(3):
        sig = np.sqrt(err**2 + sig_int**2)
        kr = keep[rows]; ridx = np.searchsorted(np.flatnonzero(keep), rows[kr])
        nk_ = int(keep.sum()); xr = np.repeat(np.arange(nk_), nB); xc = np.tile(nF + nP + np.arange(nB), nk_)
        M = sparse.csr_matrix((np.r_[w[kr] / sig[rows[kr]], w[kr] / sig[rows[kr]], (X[keep] / sig[keep, None]).ravel()],
                               (np.r_[ridx, ridx, xr], np.r_[fcol[kr], pcol[kr], xc])), shape=(nk_, nF + nP + nB))
        pri = np.r_[np.full(nF, 1 / TAU['F']), np.full(nP, 1 / TAU['P']), np.full(nB, 1e-3)]
        Aug = sparse.vstack([M, sparse.diags(pri)]).tocsr(); rhs = np.r_[(y[keep] - y0) / sig[keep], np.zeros(nF + nP + nB)]
        th = lsqr(Aug, rhs, atol=1e-9, btol=1e-9, iter_lim=3000)[0]
        pred = np.bincount(rows, weights=w * (th[fcol] + th[pcol]), minlength=len(y)) + X @ th[nF + nP:] + y0
        r = y - pred; z = r / np.sqrt(err**2 + sig_int**2)
        mad = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        sig_int = float(np.sqrt(max(mad**2 - np.median(err[keep])**2, 1e-8))); keep = np.abs(z) < 4
    info = np.asarray(M.multiply(M).sum(axis=0)).ravel()
    therr = 1 / np.sqrt(info + pri**2)
    nst_f = np.bincount(fcol[keep[rows]], minlength=nF); nst_p = np.bincount(pcode[keep[rows]], minlength=nP)
    return dict(F=th[:nF], FE=therr[:nF], NF=nst_f, P=th[nF:nF + nP], PE=therr[nF:nF + nP], NP=nst_p, puniq=puniq,
                beta=th[nF + nP:], y0=y0, sig_int=sig_int, mad=mad, spread=1.4826 * np.median(np.abs(y - np.median(y))), n=len(y), kept=keep.mean())


def predict(sol, e, keys):
    j = e.merge(keys, on=['TARGETID', 'PROG'], how='inner')
    pmap = dict(zip(sol['puniq'], sol['P']))
    pp = (j.PETAL.astype(np.int64) * 100 + j.PERIOD.astype(np.int64)).map(pmap).fillna(0.0)
    j['P_'] = j.W * (sol['F'][j.FIBER.to_numpy()] + pp.to_numpy())
    return j.groupby(['TARGETID', 'PROG']).P_.sum().rename('PRED').reset_index()


def main():
    global T0
    ap = argparse.ArgumentParser(); ap.add_argument('--residuals-only', action='store_true'); ap.add_argument('--from-residuals', action='store_true')
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args(); T0 = time.time()
    progress(0, 4, 'running', 'Calibration 2 refinement: star / stack residuals')
    rp = os.path.join(PROD, 'shape_star_residuals_desi.parquet')
    res = pd.read_parquet(rp) if a.from_residuals else residuals(a.workers)
    if a.residuals_only:
        return
    progress(1, 4, 'running', 'Calibration 2 refinement: exposure tables of the stars')
    e = exposures(res.TARGETID.unique())
    print(f'exposure rows {len(e):,} for {e.groupby(["TARGETID","PROG"]).ngroups:,} star coadds, {e.FIBER.nunique():,} fibers, '
          f'{(e.PETAL.astype(int)*100+e.PERIOD).nunique()} petal-periods; {time.time()-T0:.0f}s', flush=True)
    s68 = lambda v: float(np.diff(np.percentile(v, [16, 84]))[0] / 2)
    half = (res.TARGETID % 2 == 0).to_numpy(); ok = np.all(np.isfinite(res[[f'C{k}' for k in range(NK)] + [f'E{k}' for k in range(NK)]]), axis=1).to_numpy()
    print('hold-out test (fit on even TARGETIDs, predict odd):')
    progress(2, 4, 'running', 'Calibration 2 refinement: hold-out test')
    test = res[ok & ~half].reset_index(drop=True); curve0 = np.zeros((len(test), len(U))); curve1 = np.zeros_like(curve0)
    for k in range(NK):
        sol = fit_k(res, e, k, ok & half)
        p = predict(sol, e, test[['TARGETID', 'PROG']]); t = test.merge(p, on=['TARGETID', 'PROG'], how='left')
        X = t[['DCOL', 'DMG', 'DFEH', 'FEHNA']].to_numpy(float); full = t.PRED.fillna(0).to_numpy() + X @ sol['beta'] + sol['y0']
        r0 = t[f'C{k}'] - sol['y0'] - X @ sol['beta']; r1 = t[f'C{k}'] - full
        curve0 += np.outer(r0, V[:, k]); curve1 += np.outer(r1, V[:, k])
        print(f'  C{k}: train {sol["n"]:,} stars, sig_int {sol["sig_int"]:.4f}; test sigma68 (after the intrinsic terms) {s68(r0):.4f} -> {s68(r1):.4f} '
              f'(median error {np.median(t[f"E{k}"]):.4f}); beta (dCOL, dMG, dFEH, FEH_na) ' + ' '.join(f'{b:+.3f}' for b in sol['beta']), flush=True)
    rms0 = np.sqrt(np.mean(curve0**2, axis=1)); rms1 = np.sqrt(np.mean(curve1**2, axis=1))
    print(f'  test curve rms over 3600-9800 A (sum of C0..C5): median {np.median(rms0):.4f} -> {np.median(rms1):.4f} mag; '
          f'C2..C5 only: see per-component lines', flush=True)

    progress(3, 4, 'running', 'Calibration 2 refinement: final solution and per-coadd table')
    print('final solution (all stars):')
    sols = [fit_k(res, e, k, ok) for k in range(NK)]
    fib = pd.DataFrame({'FIBER': np.arange(5000, dtype=np.int32), 'PETAL_LOC': (np.arange(5000) // 500).astype(np.int16)})
    for k, sol in enumerate(sols):
        F = sol['F'].copy()
        for p_ in range(10):
            m = (fib.PETAL_LOC == p_).to_numpy() & (sol['NF'] > 0); F[(fib.PETAL_LOC == p_).to_numpy()] -= np.average(F[m], weights=sol['NF'][m]) if m.any() else 0.0
        F[sol['NF'] == 0] = 0.0
        fib[f'C{k}'] = F.astype(np.float32); fib[f'E{k}'] = sol['FE'].astype(np.float32)
        print(f'  C{k}: {sol["n"]:,} stars, kept {sol["kept"]*100:.1f} %, spread {sol["spread"]:.4f} -> residual {sol["mad"]:.4f}; fiber terms '
              f'rms {np.sqrt(np.mean(F[sol["NF"] > 0]**2)):.4f}, petal-period terms rms {np.sqrt(np.mean(sol["P"]**2)):.4f}', flush=True)
    fib['N_STARS'] = sols[0]['NF'].astype(np.int32)
    fib.to_parquet(os.path.join(PROD, 'response_shape_fiber.parquet'))
    pp = pd.DataFrame({'K': sols[2]['puniq']})
    for k in range(2, NK):
        d = dict(zip(sols[k]['puniq'], sols[k]['P'])); de = dict(zip(sols[k]['puniq'], sols[k]['PE']))
        v = pp.K.map(d).to_numpy(); pp[f'C{k}'] = (v - np.average(v, weights=np.maximum(sols[k]['NP'], 1))).astype(np.float32); pp[f'E{k}'] = pp.K.map(de).astype(np.float32)
    pp['N_STARS'] = sols[2]['NP'].astype(np.int32)
    pp.insert(0, 'PETAL_LOC', (pp.K // 100).astype(np.int16)); pp.insert(1, 'PERIOD', (pp.K % 100).astype(np.int16)); pp = pp.drop(columns='K')
    pp.to_parquet(os.path.join(PROD, 'response_shape_petal_period.parquet'))
    print(f'written response_shape_fiber ({(fib.N_STARS > 0).sum():,} fibers with stars) and response_shape_petal_period ({len(pp)} petal-periods); {time.time()-T0:.0f}s', flush=True)
    coadd_table(fib, pp)
    progress(4, 4, 'done', f'refinement: shape terms for {(fib.N_STARS > 0).sum():,} fibers and {len(pp)} petal-periods (see step8b_shape.txt)')


def coadd_table(fib, pp):
    """per target of the exposure tables, streamed one healpix file at a time (no 53 M-row frame in memory)."""
    out = os.path.join(PROD, 'response_shape_coadd_desi'); os.makedirs(out, exist_ok=True)
    for f in glob.glob(os.path.join(out, '*.parquet')):
        os.remove(f)
    FC = np.stack([fib[f'C{k}'].to_numpy(float) for k in range(NK)], axis=1); solved = (fib.N_STARS > 0).to_numpy()
    PC = np.zeros((1000, NK))                                   # petal * 100 + period -> C0..C5 (C0, C1 of P not published: 0)
    for r in pp.itertuples():
        PC[int(r.PETAL_LOC) * 100 + int(r.PERIOD), 2:] = [getattr(r, f'C{k}') for k in range(2, NK)]
    files = sorted(glob.glob(os.path.join(EXPFIB, '*', '*.parquet'))); buf = []; part = 0; n = 0
    cols = ['TARGETID', 'PROGRAM', 'EXPID', 'MJD', 'PETAL_LOC', 'FIBER', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z']
    for i, f in enumerate(files):
        e = _exp_rows(pd.read_parquet(f, columns=cols))
        if not len(e):
            continue
        e['W'] = e.W / e.groupby(['TARGETID', 'PROG']).W.transform('sum')
        key = e.PETAL.astype(np.int64).to_numpy() * 100 + e.PERIOD.astype(np.int64).to_numpy()
        C = FC[e.FIBER.to_numpy()] + PC[np.clip(key, 0, 999)]
        d = pd.DataFrame(C * e.W.to_numpy()[:, None], columns=[f'C{k}' for k in range(NK)])
        d['TARGETID'] = e.TARGETID.to_numpy(); d['PROG'] = e.PROG.to_numpy(); d['NEXP'] = 1; d['NEXP_SOLVED'] = solved[e.FIBER.to_numpy()].astype(int)
        g = d.groupby(['TARGETID', 'PROG']).sum().reset_index()
        g.insert(1, 'PROGRAM', np.where(g.PROG == 1, 'dark', 'bright')); g = g.drop(columns='PROG')
        g[[f'C{k}' for k in range(NK)]] = g[[f'C{k}' for k in range(NK)]].astype(np.float32); g[['NEXP', 'NEXP_SOLVED']] = g[['NEXP', 'NEXP_SOLVED']].astype(np.int16)
        buf.append(g); n += len(g)
        if len(buf) >= 2000 or i == len(files) - 1:
            pd.concat(buf, ignore_index=True).to_parquet(os.path.join(out, f'part{part:03d}.parquet')); buf = []; part += 1
            print(f'  coadd table {i+1:,}/{len(files):,} files, {n:,} targets; {time.time()-T0:.0f}s', flush=True)
    print(f'response_shape_coadd_desi: {n:,} targets in {part} parts; {time.time()-T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
