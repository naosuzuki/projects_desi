#!/usr/bin/env python3
"""Telluric step 2 (Nao 2026-11-01), the one function all three threads call: per spectrum, observed frame, before stacking.
    b = tellapply.basis('dr17')                                  # 'dr8' | 'dr17' | 'desi'; tellapply.sdss_survey(mjd) for SDSS. Cached per process
    C, a, info = b.correct(wave_obs, flux, ivar, model)          # one spectrum; any grid with wave_obs in A, ivar = 0 on pixels you mask (lines)
    flux_corr = flux / C;  ivar_corr = ivar * C**2
Factor C = 1 + T0 + sum a_k T_k on DR8, C = 1 + T0 on DR17 and DESI (APPLY_A; a is still fitted and returned for the record), multiplicative because telluric absorption scales with the flux.
T0 holds across independent populations: T0 of the galaxies improves the quasars and the reverse, in all three surveys, and T0 of both improves the white dwarfs (scatter of the per-object fits against a shifted basis: 04_star_test.py, Basis.holdout). C = 1 outside the basis: between the bands and red of 9200 A (DR8) or 9280 A (DR17, DESI; the 9280-9880 A water band is dominated by sky residuals, see 02_basis.py).
Fit on the basis pixels of r = f / (s model) - 1 - T0, weights ivar x (s model)^2, s the weighted amplitude of the model on those pixels; |r| > 0.3 and, after a first fit, 4 sigma residuals out.
Fitted jointly: a_k with the prior N(0, amp_k^2), amp_k^2 = var - var_null of step 1 (floor 0.003), and a continuous piecewise-linear curve per band (one segment per 300 A) as nuisance terms. They never enter C: they carry
the template and flux-calibration error, as the detrending of 02_basis.py does. A segment with < 50 % of its pixels usable is dropped for that spectrum. Too few pixels in all: a = 0, C = 1 + T0.
Components: var / null >= 1.5 per survey (csv_files/plots_telluric_<survey>.csv) except DROP; the prior takes faint spectra smoothly to a = 0, so there is no S/N branch.
model_has_t0=True is for stars: telluric does not smear at z ~ 0, so a star's own stack already contains the mean T0. There r = f / (s model) - 1 is fitted without T0 and C still divides T0 out.
Held-out test of this recipe: Basis.holdout, 03_fit_test.py, 04_star_test.py.
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BASIS = '/Volumes/exdisk3/sdssdesi_telluric/basis'; COEFF0 = 3.5502; DLOG = 1e-4
APPLY_A = {'dr8': True, 'dr17': False, 'desi': False}   # the a_k enter C only where they beat the eigenvectors moved off the atmosphere on held-out pixels (2026-11-02); always fitted and returned
DROP = {'dr17': ('H2O 1',)}   # above the noise but sky-line spikes at 7300, 8300, 8950 A, additive
FLOOR = 0.003; RCLIP = 0.3; NSIG = 4.0; MINCOV = 0.5; MINPIX = 30; GAP = 80.0; MINWIDTH = 30.0; MAXWIDTH = 300.0; CMIN, CMAX = 0.5, 1.5; MAXGRIDS = 256


def sdss_survey(mjd):
    return 'dr8' if mjd < 55025 else 'dr17'


def bands_of(w):
    """node wavelengths in A of each band of the sorted basis pixels w, as 02_basis.detrend: split at gaps > GAP, one segment per MAXWIDTH, narrow bands dropped."""
    e = np.r_[0, np.where(np.diff(w) > GAP)[0] + 1, len(w)]
    return [np.linspace(w[a], w[b - 1], int(np.ceil((w[b - 1] - w[a]) / MAXWIDTH)) + 1) for a, b in zip(e[:-1], e[1:]) if w[b - 1] - w[a] >= MINWIDTH]


def fill(y, ok, w):
    """y across the masked pixels of a part: the fit skips sky pixels, the correction must not jump to 1 on them."""
    return np.interp(w, w[ok], y[ok]) if ok.sum() >= 2 else y * 0


class Basis:
    def __init__(self, survey, minratio=1.5, drop=(), path=None, csv=None, shift=0):
        p = path or f'{BASIS}/telluric_{survey}.fits'; self.survey = survey; self.lat = 10**(COEFF0 + DLOG * np.arange(4447)); self.grids = {}
        T0 = fitsio.read(p, ext='MEAN').astype(float); T = fitsio.read(p, ext='EIGVEC').astype(float); ok = fitsio.read(p, ext='PIXMASK') > 0
        t = pd.read_csv(csv or os.path.join(ROOT, 'csv_files', f'plots_telluric_{survey}.csv')); self.labels = t['T'].tolist(); ratio = (t['var'] / t['var_null']).to_numpy()
        self.amp2 = np.maximum((t['var'] - t['var_null']).to_numpy(), FLOOR**2); self.sel = np.array([k for k in range(len(T)) if ratio[k] >= minratio and self.labels[k] not in (drop or DROP.get(survey, ()))], int)
        if shift: T[:, ok] = np.roll(T[:, ok], shift, axis=1)     # control of the tests: off the atmosphere, same smoothness
        self.T0, self.T, self.ok = T0, T, ok; self.bands = bands_of(self.lat[ok]); self.parts = [(n[j], n[j + 1]) for n in self.bands for j in range(len(n) - 1)]
        self.Ffull = np.zeros((len(T) + 1, len(T0)))              # the same vectors filled across masked pixels inside each part, for C
        for n in self.bands:
            i = (self.lat >= n[0]) & (self.lat <= n[-1]); w = self.lat[i]; o = ok[i]; self.Ffull[0, i] = fill(T0[i], o, w)
            for k in range(len(T)): self.Ffull[k + 1, i] = fill(T[k, i], o, w)

    def grid(self, wave):
        """everything that depends only on the wavelength array, cached (a coadd file or an SDSS plate has one grid)."""
        key = (len(wave), round(float(wave[0]), 3), round(float(wave[-1]), 3)); g = self.grids.get(key)
        if g is not None: return g
        if len(self.grids) >= MAXGRIDS: self.grids.pop(next(iter(self.grids)))
        lw = np.log10(wave); ll = COEFF0 + DLOG * np.arange(4447); F = np.array([np.interp(lw, ll, v, left=0.0, right=0.0) for v in self.Ffull])
        okw = np.interp(lw, ll, self.ok.astype(float), left=0.0, right=0.0) > 0.999; pid = np.full(len(wave), -1, int)
        for j, (lo, hi) in enumerate(self.parts): pid[(wave >= lo) & (wave <= hi)] = j
        p = np.where(okw & (pid >= 0))[0]; x = pid[p]; ks = self.sel; wp = wave[p]; nh = sum(len(n) for n in self.bands); D = np.zeros((len(ks) + nh, len(p))); D[:len(ks)] = F[1 + ks][:, p]; j0 = len(ks)
        for n in self.bands:
            m = (wp >= n[0]) & (wp <= n[-1])
            for j in range(len(n)): D[j0 + j, m] = np.interp(wp[m], n, np.eye(len(n))[j])
            j0 += len(n)
        g = dict(F=F, p=p, x=x, D=D, C0=1 + F[0], t0=F[0][p]); self.grids[key] = g; return g

    def _solve(self, D, r, W, prior):
        A = (D * W) @ D.T; A[np.diag_indices_from(A)] += prior; return np.linalg.solve(A, (D * W) @ r)

    def correct(self, wave, flux, ivar, model, model_has_t0=False, sel=True, mincov=MINCOV):
        """C on wave, a [all 7 components, 0 where not used], info. sel=False: nuisance terms only, a = 0 (the T0-only reference of the tests)."""
        g = self.grid(wave); p, x, D = g['p'], g['x'], g['D']; ns = len(self.sel); a = np.zeros(len(self.T)); nan = dict(npix=0, s=np.nan, dchi2=0.0, r=None, W=None, c=None, p=p)
        if len(p) < MINPIX: return g['C0'].copy(), a, nan
        f = flux[p]; iv = ivar[p]; m = model[p]; good = (iv > 0) & (m > 0) & np.isfinite(f) & np.isfinite(m)
        if good.sum() < MINPIX: return np.clip(g['C0'], CMIN, CMAX), a, nan
        s = np.sum((iv * f * m)[good]) / np.sum((iv * m * m)[good])
        if not (s > 0): return np.clip(g['C0'], CMIN, CMAX), a, nan
        r = np.where(good, f / np.where(good, s * m, 1.0) - 1, 0.0)
        if not model_has_t0: r = r - g['t0']
        W = np.where(good & (np.abs(r) < RCLIP), iv * (s * m)**2, 0.0)
        n = np.bincount(x, minlength=len(self.parts)); u = np.bincount(x, weights=W > 0, minlength=len(self.parts)); W[(u / np.maximum(n, 1) < mincov)[x]] = 0
        if (W > 0).sum() < MINPIX: return np.clip(g['C0'], CMIN, CMAX), a, nan
        prior = np.r_[1 / self.amp2[self.sel] if sel else np.full(ns, 1e12), np.full(D.shape[0] - ns, 1e-6)]
        c = self._solve(D, r, W, prior); res = r - c @ D; bad = res**2 * W > NSIG**2
        if bad.any(): W = np.where(bad, 0.0, W); c = self._solve(D, r, W, prior); res = r - c @ D
        n0 = c[ns:] @ D[ns:]; a[self.sel] = c[:ns]
        return self.factor_of(g, a), a, dict(npix=int((W > 0).sum()), s=float(s), dchi2=float(np.sum(W * (r - n0)**2) - np.sum(W * res**2)), r=r, W=W, c=c, p=p)


    def factor_of(self, g, a):
        """C from fitted coefficients a [all components] on a grid; what correct() returns."""
        C = g['C0'] + a[self.sel] @ g['F'][1 + self.sel] if APPLY_A[self.survey] else g['C0']; return np.clip(C, CMIN, CMAX)

    def holdout(self, wave, flux, ivar, model, block=4, model_has_t0=False):
        """chi2 gain on held-out pixels, fitted on alternating blocks of the others, against the nuisance terms alone (a = 0). Returns gain, a, number of held-out pixels; nan when not fitted."""
        g = self.grid(wave); p, D = g['p'], g['D']; fit = (np.arange(len(p)) // block) % 2 == 0; iv = ivar.copy(); iv[p[~fit]] = 0; ns = len(self.sel)
        _, a, i1 = self.correct(wave, flux, iv, model, model_has_t0, mincov=0.25); _, _, i0 = self.correct(wave, flux, iv, model, model_has_t0, sel=False, mincov=0.25)
        if i1['c'] is None or i0['c'] is None: return np.nan, a, 0
        m = model[p]; h = ~fit & (ivar[p] > 0) & (m > 0) & np.isfinite(flux[p]); r = np.where(h, flux[p] / np.where(h, i1['s'] * m, 1.0) - 1, 0.0) - (0 if model_has_t0 else g['t0']); h &= np.abs(r) < RCLIP
        w = np.where(h, ivar[p] * (i1['s'] * m)**2, 0.0); c1 = np.minimum(w * (r - i1['c'] @ D)**2, 25.0); c0 = np.minimum(w * (r - i0['c'] @ D)**2, 25.0)
        return float(c0.sum() - c1.sum()), a, int(h.sum())


_cache = {}


def basis(survey, **kw):
    key = (survey, tuple(sorted(kw.items())))
    if key not in _cache: _cache[key] = Basis(survey, **kw)
    return _cache[key]


def factor(survey, wave, flux, ivar, model, model_has_t0=False):
    """rows = spectra. wave 1-D (shared) or 2-D (one per row, e.g. rest grid x (1 + z)); NaN model = no model at that pixel: not fitted, C still applies there. Returns C, a [n, 7]."""
    b = basis(survey); flux, ivar, model = np.atleast_2d(flux), np.atleast_2d(ivar), np.atleast_2d(model); wave = np.asarray(wave)
    C = np.empty(flux.shape); a = np.zeros((len(flux), len(b.T)))
    for i in range(len(flux)):
        C[i], a[i], _ = b.correct(wave if wave.ndim == 1 else wave[i], flux[i], ivar[i], model[i], model_has_t0=model_has_t0)
    return C, a
