#!/usr/bin/env python3
"""Telluric step 2, test before any re-read (Nao 2026-11-01 plan): does fitting the step 1 eigenspectra to individual spectra remove telluric, or fit noise and template error?
On the same residuals as 02_basis.py (same sky, spike and tilt treatment), R = T0 + sum a_k T_k per object. Coefficients come from one half of the band pixels (alternating blocks of
--block pixels) and are judged on the other half, whose noise they have not seen. Models: T0 only (reference); per spectrum free; per spectrum with a Gaussian prior a_k ~ N(0, amp_k^2), amp_k^2 = var - var_null
of the step 1 table, floor 0.003; per observation (all objects of the plate / tile, same prior). Components: var / null >= --minratio of plots_telluric_<survey>.csv, and each component alone with the prior.
Score: chi2 on the held-out pixels against T0 only, per object (positive = better), chi2 per pixel capped at 25 so a sky line does not decide; error from the object scatter.
The basis itself was built from these objects, on both pixel halves; with thousands of objects per component that leak is small but not zero.
--shift N is the control: the eigenvectors are moved N band pixels, so they keep their smoothness but no longer sit on the atmosphere. What they still gain is not telluric.
Output: plots_telluric_fittest_<survey>[_shift<N>].csv in the repo root.
Usage: python3 03_fit_test.py --survey dr8|dr17|desi [--block 4] [--minratio 1.5] [--minobj 5] [--shift 20]
"""
import os, sys, argparse, importlib
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
b2 = importlib.import_module('02_basis')
from telluric_guide import band_mask, smooth_continuum_removed
CAP = 25.0; FLOOR = 0.003


def coeffs(A, b, prior):
    return np.linalg.solve(A + np.diag(prior)[None], b[..., None])[..., 0]


def score(r, W, a, T, ho):
    """chi2 on the held-out pixels per object, capped per pixel."""
    d = r[:, ho] - (a @ T[:, ho]); return np.minimum(W[:, ho] * d**2, CAP).sum(axis=1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', required=True, choices=['dr8', 'dr17', 'desi']); ap.add_argument('--block', type=int, default=4); ap.add_argument('--minratio', type=float, default=1.5); ap.add_argument('--minobj', type=int, default=5); ap.add_argument('--shift', type=int, default=0)
    a = ap.parse_args(); S = a.survey; print(S, f'shift {a.shift}' if a.shift else '', flush=True)
    full = (b2.SWAVE > 5600) & (b2.SWAVE < b2.RED[S]); px = np.where(full)[0]; wave = b2.SWAVE[px]; bands = band_mask(wave)
    R, W, plate, mjdkey, mjd, pop = b2.load(S, px); sky = b2.skymask(R, W); clean = ~bands & ~sky
    R = R - b2.baseline(R, W, wave, clean); use = bands & ~sky; ipx = np.where(use)[0]; Rb = R[:, ipx].astype(float); Wb = W[:, ipx].astype(float); w = wave[ipx]; del R, W
    key = plate * 100000 + (mjdkey - 50000 if S != 'desi' else 0)
    spk = b2.spikes(Rb, Wb, key, a.minobj, w); Wb[:, spk] = 0
    Rb, keep = b2.detrend(Rb, Wb, w, smooth_continuum_removed(w)); Rb, Wb, w, ipx = Rb[:, keep], Wb[:, keep], w[keep], ipx[keep]; idx = px[ipx]
    p = f'{b2.OUT}/telluric_{S}.fits'; T0 = fitsio.read(p, ext='MEAN').astype(float)[idx]; T = fitsio.read(p, ext='EIGVEC').astype(float)[:, idx]; ok = fitsio.read(p, ext='PIXMASK')[idx] > 0
    if a.shift: T = np.roll(T, a.shift, axis=1)
    tab = pd.read_csv(os.path.join(ROOT, f'plots_telluric_{S}.csv')); labs = tab['T'].tolist(); amp2 = np.maximum(tab['var'] - tab['var_null'], FLOOR**2).to_numpy(); sel = np.where(tab['var'] / tab['var_null'] >= a.minratio)[0]
    Wb[:, ~ok] = 0; r = np.where(Wb > 0, Rb - T0, 0.0)
    blk = (np.arange(len(w)) // a.block) % 2 == 0; fit, ho = blk, ~blk
    A = np.einsum('kp,np,lp->nkl', T[:, fit], Wb[:, fit], T[:, fit]); B = np.einsum('kp,np,np->nk', T[:, fit], Wb[:, fit], r[:, fit])
    u, inv = np.unique(key, return_inverse=True); AO = np.zeros((len(u),) + A.shape[1:]); BO = np.zeros((len(u), B.shape[1])); np.add.at(AO, inv, A); np.add.at(BO, inv, B); nobs = np.bincount(inv)[inv]
    ref = score(r, Wb, np.zeros((len(r), len(T))), T, ho); npix = (Wb[:, ho] > 0).sum(axis=1)
    raw = np.minimum(Wb[:, ho] * Rb[:, ho]**2, CAP).sum(axis=1)
    def sub(ks, prior, obs=False):
        ks = np.asarray(ks); P = np.zeros(len(T)); c = np.zeros((len(r), len(T))); pr = 1 / amp2[ks] if prior else np.full(len(ks), 1e-8)
        A_, B_ = (AO, BO) if obs else (A, B); s = coeffs(A_[:, ks][:, :, ks], B_[:, ks], pr); c[:, ks] = s[inv] if obs else s; return score(r, Wb, c, T, ho)
    models = {'no correction': raw, 'per spectrum, free': sub(sel, False), 'per spectrum, prior': sub(sel, True), 'per observation, prior': sub(sel, True, True)}
    models.update({f'alone {labs[k]}, prior': sub([k], True) for k in range(len(T))})
    rows = []
    for name, c in models.items():
        for g in list(np.unique(pop)) + ['all']:
            m = (npix > 50) & ((pop == g) if g != 'all' else np.ones(len(pop), bool)); d = ref[m] - c[m]
            rows.append(dict(survey=S, model=name, pop=g, n=int(m.sum()), dchi2_per_obj=round(d.mean(), 2), err=round(d.std() / np.sqrt(max(m.sum(), 1)), 2), dchi2_per_pix=round(d.sum() / max(npix[m].sum(), 1), 5), chi2_per_pix_T0=round(ref[m].sum() / max(npix[m].sum(), 1), 4)))
    t = pd.DataFrame(rows); print(f'  components {[labs[k] for k in sel]}; {int(ho.sum())} held-out pixels', flush=True); print(t.pivot(index='model', columns='pop', values='dchi2_per_obj').reindex(list(models)).to_string()); print(t.pivot(index='model', columns='pop', values='err').reindex(list(models)).to_string())
    out = os.path.join(ROOT, f'plots_telluric_fittest_{S}' + (f'_shift{a.shift}' if a.shift else '') + '.csv'); t.to_csv(out, index=False); print(out, flush=True)


if __name__ == '__main__':
    main()
