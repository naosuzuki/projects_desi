#!/usr/bin/env python3
"""Telluric step 1 (Nao 2026-11-01): telluric-only eigenspectra, independently for SDSS DR8 (SDSS-I/II), DR17 (BOSS / eBOSS) and DESI, from
featureless white dwarfs (01_wd_residuals.py), quasars and galaxies (Phase 0 observed-frame residuals, docs/PLAN_telluric.md).
Per object: |R| > 0.3 clipped, a smooth baseline (weighted means of 200-A nodes on the pixels outside the bands, linear between them) removed so template
continuum errors do not enter. Sky pixels: where the median weight dips below 0.7 of its running median or the median chi2 exceeds 1.5 of its running median.
Sky-line residuals are also found on the observation means (spikes()). Only pixels inside the bands of the external atmosphere (telluric_guide.band_mask) enter: the eigenvectors are zero elsewhere. The guide says where, the data say what.
Per object and band: fit offset + slope + a tau (guide) and remove the offset and slope only. Without this the components are plate-level flux-calibration tilts and the red-end
ramp (first run 2026-11-01: T1-T6 broad, correlation with tau ~ 0, though shared between populations), not the atmosphere; the absorption depth itself stays in.
Objects seen together share the atmosphere and not their template errors, so the PCA runs on means per observation: (PLATE, MJD) on SDSS, the dominant TILEID on
DESI, with >= --minobj objects. Weighted EM-PCA (step7_empca.empca) about the weighted mean T0, with W = 1 / (1 / sum W + SIG_INT^2) so a few very bright plates do not dominate.
O2 and H2O separately (the guide assigns the pixels): they vary independently (airmass, water vapour), and in one PCA the red H2O band, whose observation means scatter 2-3 x their errors,
takes every component (first runs, 2026-11-01). Per pixel W is divided by the empirical chi2 of the observation means (smoothed) so noisy pixels do not lead.
Tests per component: variance against the same PCA with the objects randomly reassigned to observations (null); coefficient rank correlation between the populations (galaxy-only
against quasar-only means of the same observations; independent template errors, shared atmosphere); share of the vector in O2 and H2O pixels; correlation
with AIRMASS (plates-dr8) and with season.
Output: /Volumes/exdisk3/sdssdesi_telluric/basis/telluric_<survey>.fits (MEAN, EIGVEC [k, 4447], PIXMASK, COEFF table), plots_telluric_<survey>.png/.csv in the repo root.
Usage: python3 02_basis.py --survey dr8|dr17|desi [--no2 3] [--nh2o 4] [--minobj 5]
"""
import os, sys, glob, argparse
import numpy as np, pandas as pd, fitsio
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, 'programs_stars')]
import starlib
from step7_empca import empca
from telluric_guide import band_mask, species, smooth_continuum_removed
SWAVE = starlib.SWAVE; P0 = '/Volumes/exdisk3/sdssdesi_telluric/phase0'; WD = '/Volumes/exdisk3/sdssdesi_telluric/wd'; OUT = '/Volumes/exdisk3/sdssdesi_telluric/basis'
SIG_INT = 0.005; CLIP = 4.0; RED = {'dr8': 9200.0, 'dr17': 9880.0, 'desi': 9880.0}; NODE = 0.012
rng = np.random.default_rng(3)


def inputs(survey):
    fam = 'desi' if survey == 'desi' else 'sdss'; L = [(f'{P0}/galaxy_{fam}.npz', 'galaxy'), (f'{P0}/quasar_{fam}.npz', 'quasar'), (f'{WD}/wd_{survey}.npz', 'wd')]
    return [(p, k) for p, k in L if os.path.exists(p)]


def load(survey, px):
    """objects of the survey as arrays on the pixels px; R (clipped), W, PLATE, MJD, POP."""
    out = []
    for p, pop in inputs(survey):
        z = np.load(p); mjd = z['MJD'].astype(float); plate = z['PLATE'].astype(np.int64)
        s = np.ones(len(mjd), bool)
        if survey != 'desi': s = (mjd < 55025) if survey == 'dr8' else (mjd >= 55025)
        if pop == 'wd': s &= z['FEATURELESS']
        s &= plate >= 0
        R = z['R'][s][:, px].astype(np.float32); W = z['W'][s][:, px].astype(np.float32); W[~np.isfinite(R) | (np.abs(R) > 0.3)] = 0; R = np.where(W > 0, R, 0)
        print(f'  {pop}: {s.sum():,} objects', flush=True); out.append((R, W, plate[s], np.round(mjd[s]).astype(np.int64) if survey != 'desi' else np.zeros(s.sum(), np.int64), mjd[s], np.array([pop] * s.sum())))
    return [np.concatenate([o[i] for o in out]) for i in range(6)]


def skymask(R, W):
    """True for sky-affected pixels of the survey."""
    wn = W / np.maximum(np.nanmedian(np.where(W > 0, W, np.nan), axis=1, keepdims=True), 1e-30)
    mw = np.median(wn, axis=0); chi = np.nanmedian(np.where(W > 0, R**2 * W, np.nan), axis=0); n = 101
    return (mw < 0.7 * median_filter(mw, n)) | (chi > 1.5 * median_filter(chi, n))


def baseline(R, W, wave, clean):
    """per object the smooth baseline from node means on the clean pixels (linear interpolation, constant beyond the ends)."""
    lw = np.log10(wave); edges = np.arange(lw[0], lw[-1] + NODE, NODE); idx = np.digitize(lw, edges) - 1; nn = len(edges); B = np.zeros(R.shape, np.float32)
    A = np.zeros((nn, R.shape[1]), np.float32); A[idx, np.arange(R.shape[1])] = clean
    sw = W @ A.T; swr = (W * R) @ A.T; xn = np.array([lw[(idx == k) & clean].mean() if ((idx == k) & clean).any() else np.nan for k in range(nn)])
    for i in range(len(R)):
        ok = (sw[i] > 0) & np.isfinite(xn)
        if ok.sum() >= 2: B[i] = np.interp(lw, xn[ok], swr[i][ok] / sw[i][ok])
    return B


def detrend(R, W, w, tau, gap=80.0, minwidth=30.0, maxwidth=300.0):
    """per object and band (runs of pixels w with gaps < gap A): weighted fit tau, 1, x; subtract the 1 and x terms. Returns R and the kept-pixel mask (narrow bands dropped)."""
    edges = np.r_[0, np.where(np.diff(w) > gap)[0] + 1, len(w)]; R = R.copy(); keep = np.ones(len(w), bool); bands = []
    cut = [edges[0]]                                      # a band wider than maxwidth gets one offset + slope per part, else one slope has to follow 900 A of flux calibration
    for i0, i1 in zip(edges[:-1], edges[1:]):
        np_ = int(np.ceil((w[i1 - 1] - w[i0]) / maxwidth)); cut += [np.searchsorted(w, w[i0] + (w[i1 - 1] - w[i0]) * (j + 1) / np_) for j in range(np_ - 1)] + [i1]
    edges = np.array(cut)
    for i0, i1 in zip(edges[:-1], edges[1:]):
        if w[i1 - 1] - w[i0] < minwidth:
            keep[i0:i1] = False; continue
        bands.append((round(w[i0]), round(w[i1 - 1]))); x = (w[i0:i1] - w[i0:i1].mean()) / 100; D = np.vstack([tau[i0:i1], np.ones(i1 - i0), x]); r = R[:, i0:i1]; v = W[:, i0:i1]
        cov = (v > 0).mean(axis=1) >= 0.5; v[~cov] = 0; u = v / np.maximum(v.max(axis=1, keepdims=True), 1e-30)            # partly covered: the fit is singular, the band is dropped for that object
        A = np.einsum('kp,np,lp->nkl', D, u, D) + 1e-6 * np.eye(3)[None]; c = np.linalg.solve(A, np.einsum('kp,np,np->nk', D, u, r)[..., None])[..., 0]
        R[:, i0:i1] = r - c[:, 1:] @ D[1:]; R[~cov, i0:i1] = 0
    print(f'  bands {bands}', flush=True); return R, keep


def group(R, W, key, minobj):
    """weighted mean per observation. key int ids; returns ids, Xbar, Wsum."""
    u, inv, n = np.unique(key, return_inverse=True, return_counts=True); keep = n >= minobj; m = np.where(keep, np.cumsum(keep) - 1, -1)[inv]; s = m >= 0; ng = keep.sum()
    SW = np.zeros((ng, R.shape[1])); SR = np.zeros((ng, R.shape[1])); np.add.at(SW, m[s], W[s]); np.add.at(SR, m[s], (W * R)[s])
    return u[keep], np.where(SW > 0, SR / np.where(SW > 0, SW, 1), 0.0), SW


def spikes(R, W, key, minobj, w, factor=1.6, iters=3):
    """pixels whose observation means scatter more than their errors say, against the neighbouring pixels: sky-line residuals. They are additive, so they do not belong in a multiplicative basis.
    Grown by 1 pixel. A broad band excess is not touched (running median of 41 pixels)."""
    _, X, SW = group(R, W, key, minobj); m = np.nanmedian(np.where(SW > 0, X, np.nan), axis=0); chi = np.nanmean(np.where(SW > 0, np.minimum((X - m)**2 * SW, 100), np.nan), axis=0); chi = np.where(np.isfinite(chi), chi, np.inf); bad = np.zeros(len(w), bool)
    for _ in range(iters):
        ref = np.array([np.median(chi[max(0, i - 20):i + 21][~bad[max(0, i - 20):i + 21]]) if (~bad[max(0, i - 20):i + 21]).any() else np.inf for i in range(len(w))])
        bad = chi > factor * ref
    return bad | np.r_[bad[1:], False] | np.r_[False, bad[:-1]]


def pca(X, SW, sp, no2, nh2o):
    """about the weighted mean, O2 and H2O pixels apart; returns T0, P [no2 + nh2o, p], C, var per k, LAB, pixok."""
    W = np.where(SW > 0, 1 / (1 / np.where(SW > 0, SW, 1) + SIG_INT**2), 0.0); T0 = np.sum(W * X, axis=0) / np.maximum(W.sum(axis=0), 1e-30); Xc = np.where(W > 0, X - T0, 0.0)
    W = np.where(np.abs(Xc) * np.sqrt(W) > CLIP, 0.0, W); Xc = np.where(W > 0, Xc, 0.0)            # a few observations at a sky line
    chi = np.sum(W * Xc**2, axis=0) / np.maximum((W > 0).sum(axis=0), 1); W = W / np.maximum(uniform_filter1d(chi, 21), 1.0)[None]
    Ps, Cs, labs = [], [], []; ok = np.zeros(X.shape[1], bool)
    for code, k, lab in [(1, no2, 'O2'), (2, nh2o, 'H2O')]:
        m = sp == code; p, c, o = empca(Xc[:, m], W[:, m], k, 25, min_cover=5, min_frac=0.3, min_wfrac=0.0); P = np.zeros((k, X.shape[1])); P[:, m] = p; ok[m] = o
        Ps.append(P); Cs.append(c); labs += [f'{lab} {i + 1}' for i in range(k)]
    P = np.vstack(Ps); C = np.hstack(Cs); return T0, P, C, np.sum(C**2, axis=0) / len(C), labs, ok


def project(X, SW, T0, P):
    W = np.where(SW > 0, 1 / (1 / np.where(SW > 0, SW, 1) + SIG_INT**2), 0.0); Xc = np.where(W > 0, X - T0, 0.0); A = np.einsum('kp,np,lp->nkl', P, W, P) + 1e-12 * np.eye(len(P))[None]
    return np.linalg.solve(A, np.einsum('kp,np,np->nk', P, W, Xc)[..., None])[..., 0]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', required=True, choices=['dr8', 'dr17', 'desi']); ap.add_argument('--no2', type=int, default=3); ap.add_argument('--nh2o', type=int, default=4); ap.add_argument('--minobj', type=int, default=5)
    a = ap.parse_args(); S = a.survey; a.nvec = a.no2 + a.nh2o; print(S, [k for _, k in inputs(S)], flush=True)
    full = (SWAVE > 5600) & (SWAVE < RED[S]); px = np.where(full)[0]; wave = SWAVE[px]; bands = band_mask(wave)
    R, W, plate, mjdkey, mjd, pop = load(S, px); sky = skymask(R, W); clean = ~bands & ~sky
    print(f'  {len(R):,} objects, {full.sum()} pixels red of 5600 A, {bands.sum()} in bands, sky masks {int((sky & bands).sum())} of them', flush=True)
    R = R - baseline(R, W, wave, clean); use = bands & ~sky; ipx = np.where(use)[0]; Rb = R[:, ipx].astype(float); Wb = W[:, ipx].astype(float); w = wave[ipx]; del R, W
    key = plate * 100000 + (mjdkey - 50000 if S != 'desi' else 0)
    spk = spikes(Rb, Wb, key, a.minobj, w); Wb[:, spk] = 0; print(f'  sky spikes on the observation means: {int(spk.sum())} of {len(w)} band pixels', flush=True)
    Rb, keep = detrend(Rb, Wb, w, smooth_continuum_removed(w)); Rb, Wb, w, ipx = Rb[:, keep], Wb[:, keep], w[keep], ipx[keep]
    ids, X, SW = group(Rb, Wb, key, a.minobj); print(f'  {len(ids):,} observations with >= {a.minobj} objects ({np.isin(key, ids).sum():,} objects)', flush=True)
    sp = species(w); T0, P, C, var, labs, ok = pca(X, SW, sp, a.no2, a.nh2o)
    # null: objects randomly reassigned to the observations
    _, Xn, SWn = group(Rb, Wb, rng.permutation(key), a.minobj); _, _, _, varn, _, _ = pca(Xn, SWn, sp, a.no2, a.nh2o)
    # orientation: positive towards the external optical depth
    tau = smooth_continuum_removed(w); s = np.sign(P @ tau); s[s == 0] = 1; P = P * s[:, None]; C = C * s[None, :]
    corr_tau = [np.corrcoef(P[k][ok & (P[k] != 0)], tau[ok & (P[k] != 0)])[0, 1] for k in range(a.nvec)]
    # populations: coefficients from galaxy-only and quasar-only means of the same observations
    cross = {}
    for x, y in [('galaxy', 'quasar'), ('galaxy', 'wd'), ('quasar', 'wd')]:
        if (pop == x).sum() and (pop == y).sum():
            mn = 3 if 'wd' not in (x, y) else 1
            ia, Xa, Wa = group(Rb[pop == x], Wb[pop == x], key[pop == x], mn); ib, Xb, Wb_ = group(Rb[pop == y], Wb[pop == y], key[pop == y], mn); c, ai, bi = np.intersect1d(ia, ib, return_indices=True)
            if len(c) >= 10:
                ca = project(Xa[ai], Wa[ai], T0, P); cb = project(Xb[bi], Wb_[bi], T0, P); cross[f'{x}x{y}'] = (len(c), [spearmanr(ca[:, k], cb[:, k])[0] for k in range(a.nvec)])
    meta = pd.DataFrame({'KEY': key, 'PLATE': plate, 'MJD': mjd}).groupby('KEY').agg(PLATE=('PLATE', 'first'), MJD=('MJD', 'mean')).reindex(ids)
    season = {k: float(np.hypot(np.corrcoef(C[:, k], np.cos(2 * np.pi * meta.MJD / 365.25))[0, 1], np.corrcoef(C[:, k], np.sin(2 * np.pi * meta.MJD / 365.25))[0, 1])) for k in range(a.nvec)}
    air = {}
    if S == 'dr8' and os.path.exists(os.path.expanduser('~/sdssredux/plates-dr8.fits')):
        pl = fitsio.read(os.path.expanduser('~/sdssredux/plates-dr8.fits'), columns=['PLATE', 'MJD', 'AIRMASS']); am = pd.DataFrame({'PLATE': pl['PLATE'].astype(int), 'MJD': pl['MJD'].astype(int), 'AIRMASS': pl['AIRMASS'].astype(float)})
        m = meta.assign(MJD=np.round(meta.MJD).astype(int)).reset_index().merge(am.drop_duplicates(['PLATE', 'MJD']), on=['PLATE', 'MJD'], how='left'); g = np.isfinite(m.AIRMASS) & (m.AIRMASS > 0.9)
        air = {k: float(spearmanr(C[g, k], m.AIRMASS[g])[0]) for k in range(a.nvec)}; print(f'  airmass for {int(g.sum())} observations', flush=True)
    tab = pd.DataFrame({'T': labs, 'var': var, 'var_null': varn, 'var_over_null': (var / varn).round(2), 'rms_coeff': np.sqrt(var).round(4), 'corr_tau': np.round(corr_tau, 3), 'season_r': [round(season[k], 3) for k in range(a.nvec)]})
    if air: tab['airmass_r'] = [round(air[k], 3) for k in range(a.nvec)]
    for n, (nc, r) in cross.items(): tab[f'r_{n}_{nc}'] = np.round(r, 3)
    print(tab.to_string(index=False), flush=True)
    os.makedirs(OUT, exist_ok=True); full_idx = px[ipx]; M = np.zeros(len(SWAVE)); M[full_idx] = T0; E = np.zeros((a.nvec, len(SWAVE))); E[:, full_idx] = P; PM = np.zeros(len(SWAVE), np.int16); PM[full_idx] = ok
    out = f'{OUT}/telluric_{S}.fits'
    with fitsio.FITS(out, 'rw', clobber=True) as F:
        F.write(M.astype(np.float32), extname='MEAN', header={'COEFF0': 3.5502, 'COEFF1': 1e-4, 'SURVEY': S}); F.write(E.astype(np.float32), extname='EIGVEC'); F.write(PM, extname='PIXMASK')
        F.write(np.rec.fromarrays([ids] + [C[:, k] for k in range(a.nvec)], names=['KEY'] + [f'C{k + 1}' for k in range(a.nvec)]), extname='COEFF')
    png = os.path.join(ROOT, f'plots_telluric_{S}.png'); tab.to_csv(png.replace('.png', '.csv'), index=False)
    nr = a.nvec + 2; fig, ax = plt.subplots(nr, 1, figsize=(15, 1.7 * nr + 1), sharex=True); guide = smooth_continuum_removed(SWAVE); xl = (5650, min(RED[S] + 20, 9900))
    def show(y, c='k'):
        f = np.where(ok & (y != 0), y, np.nan); f[np.diff(w, prepend=w[0]) > 3 * np.median(np.diff(w))] = np.nan; return f
    ax[0].plot(SWAVE, np.exp(-guide), color='0.25', lw=0.7); ax[0].set_ylabel('T (guide)'); ax[0].set_ylim(0.3, 1.03)
    ax[1].plot(w, 100 * np.where(ok, T0, np.nan), lw=0.6, color='k'); ax[1].set_ylabel('T0: mean\nresidual [%]'); ax[1].axhline(0, color='0.7', lw=0.5)
    for k in range(a.nvec):
        amp = 100 * np.sqrt(var[k]); ax[k + 2].plot(w, amp * show(P[k]), lw=0.6, color=('C0' if labs[k].startswith('O2') else 'C3') if var[k] / varn[k] > 1.5 else '0.5'); ax[k + 2].axhline(0, color='0.7', lw=0.5); ax[k + 2].set_ylabel(f'{labs[k]}\nrms x vec [%]')
        cr = '  '.join(f'{n} r {rr[k]:+.2f}' for n, (_, rr) in cross.items())
        ax[k + 2].set_title(f'{labs[k]}: var / null {var[k] / varn[k]:.2f};  corr with guide tau {corr_tau[k]:+.2f};  season r {season[k]:.2f};  {cr}' + (f';  airmass r {air[k]:+.2f}' if air else ''), fontsize=9, loc='left', pad=2)
    for x in ax:
        x.set_xlim(*xl); x.grid(alpha=0.25)
    ax[-1].set_xlabel('observed wavelength [A]')
    fig.suptitle(f'{S.upper()}: telluric-only eigenspectra (observation means of featureless WDs, quasars, galaxies: {len(ids):,} observations, {np.isin(key, ids).sum():,} objects; ' + ', '.join(f'{p} {int(np.sum(pop == p)):,}' for p in np.unique(pop)) + ')\nblue O2, red H2O: variance above 1.5 x the null of randomly reassigned objects; grey: not above it', fontsize=11)
    fig.tight_layout(); fig.savefig(png, dpi=80); plt.close(fig); print(png, out, flush=True)


if __name__ == '__main__':
    main()
