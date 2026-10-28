#!/usr/bin/env python3
"""Step 8 = Calibration 2 (COORDINATION.md 3d, 2026-10-27): the DESI spectrophotometric response per exposure and petal,
from the stars, for every class.

Each DESI star coadd has from the joint fit (step2b_spcal.py) a scale a and a tilt t: the correction that brings its
spectrum onto the Gaia-anchored photometry, as Delta m(lambda) = a + t log10(lambda/5500). A coadd is the weighted mean of
its exposures, so for star s with exposures e (weights w_e = EXPTIME x fraction of cameras IN_COADD, normalized)

    y_s = sum_e w_e theta_k(e) + eps_s,      k = (EXPID, PETAL_LOC),  y = A (mag at 6000 A) or B (mag/dex)

solved by weighted least squares (scipy lsqr) with a weak prior toward the global median (tau_A = 0.2 mag, tau_B = 0.3
mag/dex) so exposures seen by few stars stay sane; sigma_s^2 = err^2 + sigma_int^2 with sigma_int (star-to-star scatter not
in the errors: photometry, type mismatch) from the residuals; 4 sigma clipping, two passes. B uses the TILT_OK stars only.
A hold-out test (solve on half the stars, predict the other half) reports how much of the tilt scatter the response removes.

Convention: multiply the pipeline flux by  CORR(lambda) = 10^(-0.4 (A_MAG6000 + B_MAG_DEX log10(lambda/6000)))
(the same sign as step 6: flux x A_FLUX, A_FLUX = 10^(-0.4 A_MAG)). A includes the point-source normalization; for extended
sources use A_REL = A - median(A) (the exposure-to-exposure throughput part) and B.

Outputs (~/data/sdssdesi_products/stars/):
  response_exp_petal.parquet   EXPID, PETAL_LOC, A_MAG6000, A_ERR, A_REL, B_MAG_DEX, B_ERR, N_STARS_A, N_STARS_B
  response_coadd_desi.parquet  TARGETID, PROGRAM, A_MAG6000, A_REL, B_MAG_DEX, A_ERR, B_ERR, NEXP, NEXP_SOLVED
                               (every target of the exposure tables: the weighted mean over that coadd's exposures)
Usage: python3 step8_response.py [--tag desi] [--min-stars 3]
"""
import os, sys, glob, time, argparse, json
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from scipy import sparse
from scipy.sparse.linalg import lsqr

CACHE = os.path.expanduser('~/data/sdssdesi_products/stars')
EXPFIB = os.path.expanduser('~/data/sdssdesi_expfibermap')
X6000 = np.log10(6000 / 5500.)
TAU = {'A': 0.2, 'B': 0.3}


def progress(done, total, status, note):
    try:
        json.dump({'done': done, 'total': total, 'status': status, 'note': note, 'started': T0},
                  open(os.path.expanduser('~/data/sdssdesi_status/step_star_response.json'), 'w'))
    except Exception:
        pass


def load_exposures():
    cols = ['TARGETID', 'PROGRAM', 'EXPID', 'PETAL_LOC', 'EXPTIME', 'IN_COADD_B', 'IN_COADD_R', 'IN_COADD_Z']
    files = sorted(glob.glob(os.path.join(EXPFIB, '*', '*.parquet')))
    parts = []
    for i, f in enumerate(files):
        d = pd.read_parquet(f, columns=cols)
        frac = (d.IN_COADD_B.astype(np.float32) + d.IN_COADD_R.astype(np.float32) + d.IN_COADD_Z.astype(np.float32)) / 3
        d = pd.DataFrame({'TARGETID': d.TARGETID.to_numpy(np.int64), 'PROG': (d.PROGRAM.astype(str) == 'dark').to_numpy(np.int8),
                          'K': d.EXPID.to_numpy(np.int64) * 10 + d.PETAL_LOC.to_numpy(np.int64),
                          'W': (d.EXPTIME.to_numpy(np.float32) * frac.to_numpy())})
        parts.append(d[d.W > 0])
        if i % 2000 == 0:
            print(f'  exposure tables {i:,}/{len(files):,}', flush=True)
    e = pd.concat(parts, ignore_index=True)
    e['W'] = e.W / e.groupby(['TARGETID', 'PROG']).W.transform('sum')
    return e


def solve(rows, kidx, w, y, err, nk, prior, tau, sig_int):
    """weighted LSQ for theta[nk]; rows/kidx/w: sparse design (star index, parameter index, weight)."""
    ns = len(y); sig = np.sqrt(err**2 + sig_int**2)
    M = sparse.csr_matrix((w / sig[rows], (rows, kidx)), shape=(ns, nk))
    R = sparse.identity(nk, format='csr') / tau
    Aug = sparse.vstack([M, R]).tocsr(); rhs = np.concatenate([(y - prior) / sig, np.zeros(nk)])
    theta = lsqr(Aug, rhs, atol=1e-8, btol=1e-8, iter_lim=2000)[0] + prior
    info = np.asarray(M.multiply(M).sum(axis=0)).ravel()            # diagonal of the normal matrix (data part)
    return theta, 1 / np.sqrt(info + 1 / tau**2)


def fit(st, e, comp, mask, label, holdout=None):
    """comp: 'A' or 'B'; st: stars with Y_A/Y_B/E_A/E_B; returns theta per K (Series) and diagnostics."""
    s = st[mask].reset_index(drop=True); s['S'] = np.arange(len(s))
    j = e.merge(s[['TARGETID', 'PROG', 'S']], on=['TARGETID', 'PROG'], how='inner')
    kcode, kuniq = pd.factorize(j.K)
    y = s[f'Y_{comp}'].to_numpy(float); err = s[f'E_{comp}'].to_numpy(float); prior = float(np.median(y))
    rows, w = j.S.to_numpy(), j.W.to_numpy(float)
    keep = np.ones(len(y), bool); sig_int = 0.5 * float(np.std(y))
    for it in range(3):
        good_rows = keep[rows]
        th, therr = solve(np.searchsorted(np.flatnonzero(keep), rows[good_rows]), kcode[good_rows], w[good_rows], y[keep], err[keep],
                          len(kuniq), prior, TAU[comp], sig_int)
        pred = np.bincount(rows, weights=w * th[kcode], minlength=len(y))
        res = y - pred; z = res / np.sqrt(err**2 + sig_int**2)
        mad = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        sig_int = float(np.sqrt(max(mad**2 - np.median(err[keep])**2, 1e-6)))
        keep = np.abs(z) < 4
    nst = np.bincount(kcode[keep[rows]], minlength=len(kuniq))
    out = pd.DataFrame({'K': kuniq, f'TH_{comp}': th, f'ERR_{comp}': therr, f'N_{comp}': nst})
    print(f'  {label} {comp}: {len(y):,} stars, {len(kuniq):,} exposure-petals, sigma_int {sig_int:.3f}, kept {keep.mean()*100:.1f} %; '
          f'scatter of y {1.4826*np.median(np.abs(y-np.median(y))):.3f} -> residual {mad:.3f}', flush=True)
    return out


def predict(e, th, comp, keys):
    j = e.merge(keys, on=['TARGETID', 'PROG'], how='inner').merge(th[['K', f'TH_{comp}']], on='K', how='left')
    j['WT'] = j.W * j[f'TH_{comp}'].notna(); j['P'] = j.W * j[f'TH_{comp}'].fillna(0)
    g = j.groupby(['TARGETID', 'PROG']).agg(P=('P', 'sum'), WT=('WT', 'sum')).reset_index()
    g[f'PRED_{comp}'] = np.where(g.WT > 0, g.P / g.WT, np.nan)
    return g[['TARGETID', 'PROG', f'PRED_{comp}']]


def main():
    global T0
    ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='desi'); ap.add_argument('--min-stars', type=int, default=3)
    a = ap.parse_args(); T0 = time.time()
    progress(0, 4, 'running', 'Calibration 2: reading the joint fit and the exposure tables')
    st = pd.read_parquet(os.path.join(CACHE, f'spcal_{a.tag}.parquet'), columns=['SID', 'A_MAG', 'A_ERR', 'TILT', 'TILT_ERR', 'TILT_OK', 'SNR_ALL'])
    sid = st.SID.str.split('-', n=1, expand=True)
    st['TARGETID'] = sid[0].astype(np.int64); st['PROG'] = (sid[1] == 'dark').astype(np.int8)
    st['Y_A'] = st.A_MAG + st.TILT * X6000; st['E_A'] = np.sqrt(st.A_ERR**2 + (st.TILT_ERR * X6000)**2)
    st['Y_B'] = st.TILT; st['E_B'] = st.TILT_ERR
    okA = np.isfinite(st.Y_A) & np.isfinite(st.E_A); okB = okA & st.TILT_OK.fillna(False).astype(bool) & np.isfinite(st.E_B)
    print(f'stars: {len(st):,} coadds, A usable {okA.sum():,}, B usable (TILT_OK) {okB.sum():,}; {time.time()-T0:.0f}s', flush=True)
    e = load_exposures()
    print(f'exposure rows {len(e):,} ({e.K.nunique():,} exposure-petals, {e.groupby(["TARGETID","PROG"]).ngroups:,} coadds); {time.time()-T0:.0f}s', flush=True)
    progress(1, 4, 'running', 'Calibration 2: hold-out test')

    # hold-out: fit on even TARGETIDs, predict the odd ones
    half = (st.TARGETID % 2 == 0).to_numpy()
    print('hold-out test (fit on even TARGETIDs, predict odd):')
    for comp, ok in (('A', okA), ('B', okB)):
        th = fit(st, e, comp, ok & half, 'train')
        test = st[ok & ~half][['TARGETID', 'PROG', f'Y_{comp}', f'E_{comp}']]
        p = predict(e, th[th[f'N_{comp}'] >= a.min_stars], comp, test[['TARGETID', 'PROG']]).merge(test, on=['TARGETID', 'PROG'])
        m = p[f'PRED_{comp}'].notna(); r0 = p[f'Y_{comp}'][m]; r1 = (p[f'Y_{comp}'] - p[f'PRED_{comp}'])[m]
        s68 = lambda v: float(np.diff(np.percentile(v, [16, 84]))[0] / 2)
        print(f'  test {comp}: {m.sum():,} of {len(p):,} predicted; sigma68 {s68(r0):.3f} -> {s68(r1):.3f} (median error {np.median(p[f"E_{comp}"][m]):.3f})', flush=True)

    progress(2, 4, 'running', 'Calibration 2: final solution on all stars')
    print('final solution (all stars):')
    thA = fit(st, e, 'A', okA, 'all'); thB = fit(st, e, 'B', okB, 'all')
    th = thA.merge(thB, on='K', how='outer')
    for c in ('A', 'B'):
        bad = th[f'N_{c}'].fillna(0) < a.min_stars
        th.loc[bad, [f'TH_{c}', f'ERR_{c}']] = np.nan
    resp = pd.DataFrame({'EXPID': (th.K // 10).astype(np.int64), 'PETAL_LOC': (th.K % 10).astype(np.int16),
                         'A_MAG6000': th.TH_A, 'A_ERR': th.ERR_A, 'A_REL': th.TH_A - np.nanmedian(th.TH_A),
                         'B_MAG_DEX': th.TH_B, 'B_ERR': th.ERR_B,
                         'N_STARS_A': th.N_A.fillna(0).astype(np.int32), 'N_STARS_B': th.N_B.fillna(0).astype(np.int32)})
    resp.to_parquet(os.path.join(CACHE, 'response_exp_petal.parquet'))
    q = lambda v: (np.nanmedian(v), 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v))))
    print(f'response_exp_petal: {len(resp):,} exposure-petals; solved A {resp.A_MAG6000.notna().sum():,}, B {resp.B_MAG_DEX.notna().sum():,}; '
          f'A median {q(resp.A_MAG6000)[0]:+.3f} spread {q(resp.A_MAG6000)[1]:.3f} mag; B median {q(resp.B_MAG_DEX)[0]:+.3f} spread {q(resp.B_MAG_DEX)[1]:.3f} mag/dex; '
          f'median N_STARS_B {np.median(resp.N_STARS_B[resp.N_STARS_B > 0]):.0f}', flush=True)
    for p_ in range(10):
        r = resp[resp.PETAL_LOC == p_]
        print(f'  petal {p_}: B median {np.nanmedian(r.B_MAG_DEX):+.3f}, A_REL median {np.nanmedian(r.A_REL):+.3f}')

    progress(3, 4, 'running', 'Calibration 2: per-coadd table for every target')
    tk = th.dropna(subset=['TH_A', 'TH_B'], how='all')[['K', 'TH_A', 'TH_B', 'ERR_A', 'ERR_B']]
    j = e.merge(tk, on='K', how='left')
    for c in ('A', 'B'):
        j[f'w{c}'] = j.W * j[f'TH_{c}'].notna(); j[f'p{c}'] = j.W * j[f'TH_{c}'].fillna(0); j[f'v{c}'] = (j.W * j[f'ERR_{c}'].fillna(0))**2
    j['solved'] = j.TH_B.notna() | j.TH_A.notna()
    g = j.groupby(['TARGETID', 'PROG']).agg(pA=('pA', 'sum'), wA=('wA', 'sum'), vA=('vA', 'sum'), pB=('pB', 'sum'), wB=('wB', 'sum'), vB=('vB', 'sum'),
                                             NEXP=('K', 'size'), NEXP_SOLVED=('solved', 'sum')).reset_index()
    with np.errstate(all='ignore'):
        co = pd.DataFrame({'TARGETID': g.TARGETID, 'PROGRAM': np.where(g.PROG == 1, 'dark', 'bright'),
                           'A_MAG6000': np.where(g.wA > 0, g.pA / g.wA, np.nan), 'B_MAG_DEX': np.where(g.wB > 0, g.pB / g.wB, np.nan),
                           'A_ERR': np.where(g.wA > 0, np.sqrt(g.vA) / g.wA, np.nan), 'B_ERR': np.where(g.wB > 0, np.sqrt(g.vB) / g.wB, np.nan),
                           'NEXP': g.NEXP.astype(np.int16), 'NEXP_SOLVED': g.NEXP_SOLVED.astype(np.int16)})
    co['A_REL'] = co.A_MAG6000 - np.nanmedian(resp.A_MAG6000)
    co.to_parquet(os.path.join(CACHE, 'response_coadd_desi.parquet'))
    print(f'response_coadd_desi: {len(co):,} coadds; with B {co.B_MAG_DEX.notna().mean()*100:.1f} %, with A {co.A_MAG6000.notna().mean()*100:.1f} %; {time.time()-T0:.0f}s', flush=True)
    progress(4, 4, 'done', f'{resp.B_MAG_DEX.notna().sum():,} exposure-petals solved; {len(co):,} coadds; B spread {q(resp.B_MAG_DEX)[1]:.3f} mag/dex')


if __name__ == '__main__':
    main()
