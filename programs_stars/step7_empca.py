#!/usr/bin/env python3
"""Step 7 (Nao, 2026-10-25): weighted EM-PCA eigenspectra from the HR-cell stacks (+ high-S/N individuals), per temperature class.

Normalization (decided 2026-10-25): every input is divided by its synthetic Gaia-G flux (G flux = 1; ivar x factor^2); the
factor is STORED (NORM table: NORM_GFLUX = photon-weighted mean F_lambda through the G band at 10 pc, M_G_SYN (Vega), M_G_SYN_AB,
M_G0_CELL) so absolute flux = normalized spectrum x NORM_GFLUX exactly. Linear flux; dereddened, RV-corrected inputs (the stacks
of step 6); weighted mean subtracted and kept as component 0; zero weight on uncovered pixels.
Classes by the cell's HR position: WD (M_G0 > 9 and (BP-RP)_0 < 0.8), HOT ((BP-RP)_0 < 0.3, not WD), AF (0.3-0.7), GK (0.7-2.0),
M (> 2.0). Weights: each stack counts once (its ivar renormalized to a median of 1 over its good pixels); --individuals K adds up to
K stars with S/N > --snmin per stack, read from the spPlates with step6_stack's machinery, with ivar scaled so that all
individuals together weigh as much as the stacks of that class.
EM-PCA (Bailey 2012) in numpy: alternate (E) per-object coefficients by weighted LSQ and (M) per-pixel eigenvector values by
weighted LSQ, orthonormalize, --niter iterations, --nvec components.
Output: ~/data/sdssdesi_stars/empca_<tag>_<class>.fits with HDUs MEAN [npix], EIGVEC [nvec, npix], COEFF [nobs, nvec], NORM table,
VARIANCE (R^2 of the weighted fit with the first k vectors), PIXMASK (pixels in the basis), header COEFF0/COEFF1; quick-look PNG per class.
Usage: python3 step7_empca.py --survey sdss --tag sdss_0.5x0.1_gspphot [--nvec 10] [--niter 25] [--individuals 5] [--snmin 20]
"""
import os, sys, argparse, time
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
passB = starlib.passB
OUTDIR = os.path.expanduser('~/data/sdssdesi_stars')
SWAVE = starlib.SWAVE; NPIX = len(SWAVE); DLOG = starlib.DLOG; SRID_MIN = starlib.SRID_MIN
CLASSES = ['WD', 'HOT', 'AF', 'GK', 'M']


def star_class(mg0, col0):
    if np.isfinite(mg0) and mg0 > 9 and col0 < 0.8: return 'WD'
    if col0 < 0.3: return 'HOT'
    if col0 < 0.7: return 'AF'
    if col0 < 2.0: return 'GK'
    return 'M'


def g_flux(flux, good):
    """synthetic Gaia-G flux of spectra on SWAVE (photon-weighted mean F_lambda through the G band; flat padding beyond the
    covered range as in starlib.gaia_mags). Returns NORM_GFLUX, M_G_SYN (Vega), M_G_SYN_AB."""
    n = flux.shape[0]
    # padding level = mean of the outermost 300 GOOD pixels of each spectrum (SDSS-I spectra start at 3800 A, not at the grid edge)
    lo = np.full(n, np.nan); hi = np.full(n, np.nan)
    for i in range(n):
        idx = np.where(good[i])[0]
        if len(idx) > 10:
            lo[i] = flux[i, idx[:300]].mean(); hi[i] = flux[i, idx[-300:]].mean()
    dw = float(np.median(np.diff(SWAVE)))
    pad_lo = np.arange(starlib.PAD_MIN, SWAVE[0] - 0.5 * dw, dw); pad_hi = np.arange(SWAVE[-1] + dw, starlib.PAD_MAX, dw)
    wz = np.concatenate([pad_lo, SWAVE, pad_hi])
    # inside the covered range, fill masked pixels by interpolation so the band integral is not biased by gaps
    filled = flux.copy()
    for i in range(n):
        g = good[i]
        if g.sum() > 10:
            filled[i, ~g] = np.interp(SWAVE[~g], SWAVE[g], flux[i, g])
    fz = np.concatenate([np.repeat(np.nan_to_num(lo)[:, None], len(pad_lo), axis=1), filled, np.repeat(np.nan_to_num(hi)[:, None], len(pad_hi), axis=1)], axis=1)
    gz = np.ones_like(fz, bool)
    flt = passB.filt('gaiadr3-G'); wt, den, inband = passB.band_weights(flt, wz)
    m_ab, _ = passB.ab_mags(flt, wz, fz, gz, min_range=0.0, min_good=0.0)
    gf = (fz * wt[None, :]).sum(axis=1) / wt.sum()                      # photon-weighted mean F_lambda in the band
    ok = np.isfinite(lo) & np.isfinite(hi) & (gf > 0)
    return np.where(ok, gf, np.nan), np.where(ok, m_ab - starlib.GAIA_AB_MINUS_VEGA['G'], np.nan), np.where(ok, m_ab, np.nan)


def empca(X, W, nvec, niter, min_cover=5, min_frac=0.5, min_wfrac=0.3, n_ref=None):
    """STABLE INTERFACE (imported by programs_templates/10_galaxy_empca.py): empca(X, W, nvec, niter, min_cover, min_frac,
    min_wfrac[, n_ref]) -> (P, C, pixok); explained(X, W, P) -> R^2 per truncation. Add keywords only with defaults.
    Weighted EM-PCA (Bailey 2012). X [n, p] data (mean already subtracted), W [n, p] weights (0 = missing).
    Pixels with weight in fewer than max(min_cover, min_frac*n) objects are excluded (W = 0, eigenvectors 0 there): the
    per-pixel M-step is ill-conditioned where only one or two objects contribute (e.g. the blue end covered by BOSS plates only).
    Init: SVD of the data with missing pixels set to 0. Returns P [nvec, p] orthonormal, C [n, nvec], pixmask [p]."""
    n, p = X.shape
    # a pixel enters the basis only where at least min_frac of the objects have it AND the median weight over all objects is
    # >= min_wfrac of the typical weight: pixels covered only by a few low-S/N inputs (the blue end from BOSS plates only,
    # the red end beyond 9200 A for SDSS-I) would otherwise dominate the unit-norm eigenvectors with noise
    # the pixel criterion is evaluated on the first n_ref rows (the stacks) only: the individuals carry small weights by design
    Wr = W[:n_ref] if n_ref else W; nr = Wr.shape[0]
    cover = (Wr > 0).sum(axis=0); wmed = np.median(Wr, axis=0); wtyp = np.median(Wr[Wr > 0]) if (Wr > 0).any() else 1.0
    pixok = (cover >= max(min_cover, int(min_frac * nr))) & (wmed >= min_wfrac * wtyp)
    W = np.where(pixok[None, :], W, 0.0); X = np.where(W > 0, X, 0.0)
    U, s, Vt = np.linalg.svd(X, full_matrices=False); P = Vt[:nvec].copy()
    C = np.zeros((n, nvec)); eye = np.eye(nvec)
    for it in range(niter):
        A = np.einsum('kp,np,lp->nkl', P, W, P); A += 1e-6 * np.einsum('nkk->n', A)[:, None, None] * eye[None] + 1e-12 * eye[None]
        C = np.linalg.solve(A, np.einsum('kp,np,np->nk', P, W, X)[..., None])[..., 0]
        A2 = np.einsum('nk,np,nl->pkl', C, W, C); A2 += 1e-6 * np.einsum('pkk->p', A2)[:, None, None] * eye[None] + 1e-12 * eye[None]
        P = np.linalg.solve(A2, np.einsum('nk,np,np->pk', C, W, X)[..., None])[..., 0].T
        P = np.where(pixok[None, :], P, 0.0)
        Q, R = np.linalg.qr(P.T); P = Q.T; C = C @ R.T
        order = np.argsort(-np.einsum('nk,nk->k', C, C)); P = P[order]; C = C[:, order]
    # final coefficients consistent with the final vectors
    A = np.einsum('kp,np,lp->nkl', P, W, P) + 1e-12 * eye[None]
    C = np.linalg.solve(A, np.einsum('kp,np,np->nk', P, W, X)[..., None])[..., 0]
    return P, C, pixok


def explained(X, W, P):
    """R^2 of the weighted fit using the first k vectors, k = 1..nvec, with coefficients refitted for each truncation."""
    X = np.where(W > 0, X, 0.0); tot = np.sum(W * X**2); out = []
    for k in range(1, P.shape[0] + 1):
        Pk = P[:k]; A = np.einsum('kp,np,lp->nkl', Pk, W, Pk) + 1e-12 * np.eye(k)[None]
        Ck = np.linalg.solve(A, np.einsum('kp,np,np->nk', Pk, W, X)[..., None])[..., 0]
        out.append(1 - np.sum(W * (X - Ck @ Pk)**2) / tot)
    return np.array(out)


def load_stacks(tag):
    F = fitsio.FITS(os.path.join(OUTDIR, f'stacks_{tag}.fits'))
    fl = F['FLUX'].read().astype(float); err = F['ERR'].read().astype(float); n = F['NSPEC'].read(); arr = F['CELLS'].read()
    cells = pd.DataFrame({c: (np.char.strip(arr[c].astype(str)) if arr[c].dtype.kind in 'SU' else arr[c].astype(arr[c].dtype.newbyteorder('='))) for c in arr.dtype.names})
    good = (n > 0) & np.isfinite(fl) & np.isfinite(err) & (err > 0)
    with np.errstate(all='ignore'):
        ivar = np.where(good, 1.0 / err**2, 0.0)
    return fl, ivar, good, cells


def load_individuals(tag, cells, k_per_stack, snmin, survey):
    """up to k stars per stack with S/N > snmin, read and calibrated exactly as in step 6 (absolute flux at 10 pc)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location('step6', os.path.join(HERE, 'step6_stack.py')); S6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(S6)
    mem = pd.read_csv(os.path.join(OUTDIR, f'members_{tag}.csv'))
    mem = mem[mem.SNR_ALL > snmin].sort_values('SNR_ALL', ascending=False).groupby('STACK').head(k_per_stack)
    selp = os.path.join(S6.CACHE, f'step6_selection_{tag}.parquet')
    if survey == 'desi' and os.path.exists(selp):          # the step-6 selection has every column; spectra from star_hook's store
        import star_hook
        df = pd.read_parquet(selp); df = df[df.SID.isin(mem.SID)]
        specs, ws, ids, stacks = [], [], [], []
        for row, f, iv, mk, rid0 in S6.iter_store(df, 0, star_hook.STORE):
            s, w = S6.star_spectrum(row, f, iv, mk, rid0); specs.append(s); ws.append(w); ids.append(row.SID); stacks.append(row.STACK)
        return (np.array(specs), np.array(ws), ids, stacks) if specs else (np.zeros((0, NPIX)), np.zeros((0, NPIX)), [], [])
    files = sorted(__import__('glob').glob(os.path.join(S6.PASSDIR[survey], '*', '*.parquet')))
    import pyarrow.parquet as pq
    need = ['FILE', 'Z', 'GAIA_DG', 'CAL_A', 'EBV_USED', 'EBV_SFD', 'ABS_SCALE', 'SNR_ALL'] + (['PLATE', 'MJD', 'FIBERID', 'SURVEY'] if survey == 'sdss' else ['TARGETID', 'PROGRAM'])
    files = [f for f in files if 'FILE' in pq.ParquetFile(f).schema_arrow.names]
    df = pd.concat([pd.read_parquet(f, columns=[c for c in need if c in pq.ParquetFile(f).schema_arrow.names]) for f in files], ignore_index=True)
    df['SID'] = (df.PLATE.astype(str) + '-' + df.MJD.astype(str) + '-' + df.FIBERID.astype(str)) if survey == 'sdss' else starlib.desi_sid(df)
    df = df[df.SID.isin(mem.SID)].merge(mem[['SID', 'STACK', 'CELL']], on='SID', how='inner')
    sp = os.path.join(S6.CACHE, f'spcal_{survey}.parquet')
    if os.path.exists(sp):
        df = df.merge(pd.read_parquet(sp, columns=['SID', 'A_FLUX', 'TILT', 'TILT_OK']), on='SID', how='left')
    else:
        df['A_FLUX'] = np.nan; df['TILT'] = np.nan; df['TILT_OK'] = False
    it = S6.iter_sdss if survey == 'sdss' else S6.iter_desi
    specs, ws, ids, stacks = [], [], [], []
    for row, f, iv, mk, rid0 in it(df, 0):
        s, w = S6.star_spectrum(row, f, iv, mk, rid0); specs.append(s); ws.append(w); ids.append(row.SID); stacks.append(row.STACK)
    if not specs:
        return np.zeros((0, NPIX)), np.zeros((0, NPIX)), [], []
    return np.array(specs), np.array(ws), ids, stacks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='sdss', choices=['sdss', 'desi']); ap.add_argument('--tag', required=True)
    ap.add_argument('--nvec', type=int, default=10); ap.add_argument('--niter', type=int, default=25)
    ap.add_argument('--individuals', type=int, default=0, help='stars per stack with S/N > --snmin added to the PCA'); ap.add_argument('--snmin', type=float, default=20.0)
    ap.add_argument('--classes', default=','.join(CLASSES)); ap.add_argument('--plot', action='store_true')
    ap.add_argument('--min-wfrac', type=float, default=0.3, help='weight cut of empca (meant for the BOSS-only blue end); 0 for DESI, where every pixel is covered: the cut drops the red end of M stars, whose photon-noise weights are low there at high S/N; SDSS M: 0.1 (keeps 3675-9202 A, masks the bright sky lines; 0 lets the noisy BOSS ends dominate)')
    a = ap.parse_args(); t0 = time.time()
    fl, ivar, good, cells = load_stacks(a.tag)
    cls = np.array([star_class(r.MG0_MED, r.COL0_MED) for _, r in cells.iterrows()])
    gf, mg, mg_ab = g_flux(fl, good)
    print(f'{a.tag}: {len(cells)} stacks; classes {dict(zip(*np.unique(cls, return_counts=True)))}; G-flux normalization ok for {np.isfinite(gf).sum()}; {time.time()-t0:.0f}s', flush=True)
    ind_f, ind_w, ind_ids, ind_stacks = (load_individuals(a.tag, cells, a.individuals, a.snmin, a.survey) if a.individuals else (np.zeros((0, NPIX)), np.zeros((0, NPIX)), [], []))
    if len(ind_ids):
        igood = ind_w > 0; igf, img, img_ab = g_flux(ind_f, igood)
        print(f'individuals: {len(ind_ids)} stars (<= {a.individuals} per stack, S/N > {a.snmin}); {time.time()-t0:.0f}s', flush=True)
    for c in a.classes.split(','):
        sel = np.where((cls == c) & np.isfinite(gf))[0]
        if len(sel) < 3:
            print(f'  {c}: {len(sel)} stacks - skipped'); continue
        X = fl[sel] / gf[sel, None]; W = ivar[sel] * gf[sel, None]**2 * good[sel]
        # each stack counts once: renormalize its weights to a median of 1 over its good pixels
        med = np.array([np.median(W[i][W[i] > 0]) if (W[i] > 0).any() else 1.0 for i in range(len(sel))]); W = W / med[:, None]
        norm_rows = [dict(ID=cells.STACK.iloc[i], KIND='stack', CLASS=c, NORM_GFLUX=gf[i], M_G_SYN=mg[i], M_G_SYN_AB=mg_ab[i], M_G0_CELL=cells.MG0_MED.iloc[i], COL0_CELL=cells.COL0_MED.iloc[i], FEH_CELL=cells.FEH_MED.iloc[i], N=int(cells.N.iloc[i])) for i in sel]
        if len(ind_ids):
            cmap = dict(zip(cells.STACK, cls)); isel = [j for j, s in enumerate(ind_stacks) if cmap.get(s) == c and np.isfinite(igf[j])]
            if isel:
                Xi = ind_f[isel] / igf[isel, None]; Wi = ind_w[isel] * igf[isel, None]**2
                Wi = Wi * (W.sum() / max(Wi.sum(), 1e-30))                     # all individuals together weigh as much as the stacks
                X = np.vstack([X, Xi]); W = np.vstack([W, Wi])
                norm_rows += [dict(ID=ind_ids[j], KIND='star', CLASS=c, NORM_GFLUX=igf[j], M_G_SYN=img[j], M_G_SYN_AB=img_ab[j], M_G0_CELL=np.nan, COL0_CELL=np.nan, FEH_CELL=np.nan, N=1) for j in isel]
        X = np.nan_to_num(X); W = np.where(np.isfinite(W), W, 0.0)
        mean = np.sum(W * X, axis=0) / np.where(W.sum(axis=0) > 0, W.sum(axis=0), np.nan)
        mean = np.where(np.isfinite(mean), mean, 0.0)
        Xc = X - mean[None, :]
        nvec = min(a.nvec, len(X) - 1)
        P, C, pixok = empca(Xc, W, nvec, a.niter, min_wfrac=a.min_wfrac, n_ref=len(sel))
        var = explained(np.where(pixok[None, :], Xc, 0.0), np.where(pixok[None, :], W, 0.0), P)
        mean = np.where(pixok, mean, 0.0)
        out = os.path.join(OUTDIR, f'empca_{a.tag}_{c}.fits')
        if os.path.exists(out): os.remove(out)
        hdr = {'COEFF0': SRID_MIN * DLOG, 'COEFF1': DLOG, 'CLASS': c, 'NORM': 'synthetic Gaia G flux = 1 (see NORM table)', 'NVEC': nvec, 'NITER': a.niter, 'NOBS': len(X), 'TAG': a.tag}
        with fitsio.FITS(out, 'rw') as F:
            F.write(mean.astype(np.float32), header=hdr, extname='MEAN'); F.write(P.astype(np.float32), extname='EIGVEC')
            F.write(C.astype(np.float32), extname='COEFF'); F.write(var.astype(np.float32), extname='VARIANCE'); F.write(pixok.astype(np.int8), extname='PIXMASK')
            F.write(pd.DataFrame(norm_rows).to_records(index=False), extname='NORM')
        print(f'  {c}: {len(sel)} stacks + {len(X)-len(sel)} stars, {pixok.sum()} of {NPIX} pixels -> {nvec} components; R2 with 1/2/3/5/{nvec}: ' +
              ' '.join(f'{v:.4f}' for v in var[[0, 1, 2, min(4, nvec-1), nvec-1]]) + f'; {out}  {time.time()-t0:.0f}s', flush=True)
        if a.plot:
            import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
            fig, axes = plt.subplots(5, 1, figsize=(11, 10), sharex=True)
            axes[0].plot(SWAVE, mean, 'k', lw=0.6); axes[0].set_title(f'{c}: mean (G flux = 1), {len(sel)} stacks + {len(X)-len(sel)} stars', fontsize=9, loc='left')
            for k in range(min(4, nvec)):
                axes[k + 1].plot(SWAVE, P[k], lw=0.6); axes[k + 1].set_title(f'eigenvector {k+1}: {var[k] - (var[k-1] if k else 0):.3f} of the weighted variance', fontsize=9, loc='left')
            axes[-1].set_xlabel('wavelength [A]'); fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, f'empca_{a.tag}_{c}.png'), dpi=100); plt.close(fig)


if __name__ == '__main__':
    main()
