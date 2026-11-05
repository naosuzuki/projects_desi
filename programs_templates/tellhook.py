#!/usr/bin/env python3
"""Galaxy side of the telluric correction (Nao 2026-11-02: "apply the telluric correction the Stellar thread is finding and redo
the stacks and PCA"). Step 2 of the telluric plan (docs/PLAN_telluric.md, memory telluric-three-step-plan): every spectrum is
divided by C(lambda_obs) = 1 + T0 + sum_k a_k T_k BEFORE stacking, with the per-survey telluric eigenspectra of
programs_telluric (02_basis.py) and the shared fit of programs_telluric/tellapply.py (Stars thread; a_k per spectrum with a
Gaussian prior, offset + slope per band as nuisance parameters, C = 1 where PIXMASK = 0).
What this module adds is the galaxy model the residual f / model - 1 needs: each spectrum's weighted least-squares fit with the
11-vector v3 galaxy basis (galaxy_empca_v3pre_G.fits, 1862-8119 A rest, lines masked +-LINE_KMS), on the REST grid. The
correction is fitted and applied on the rest grid with wave_obs = RWAVE (1 + z) per spectrum: f / model - 1 is the same in
both frames and the second sub-pixel rebin is saved. Where the basis does not cover (rest > 8119 A, i.e. the 8950-9800 A water
band at z < 0.1) the model is NaN and tellapply falls back to C = 1 + T0.
Record (Nao: keep every factor): per spectrum TELL_A<k> (a_k, all 7 components, 0 where not used), TELL_NFIT (pixels in the
basis fit), TELL_CHI2R (chi2 of the basis fit per pixel); the pass tables carry them. Hooked into 08_passC_stack.py (--telluric, DESI) and 11_passG_sdss.py
(--telluric, dr8 / dr17) in process_chunk / process_plate right after the rest-frame step, before lines and magnitudes.
Usage (test): python3 tellhook.py --survey desi --n 200   (fits stored v2 rest-frame spectra and prints the a_k statistics)
"""
import os, sys, importlib
import numpy as np
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'programs_telluric'))
BASIS = os.path.expanduser('~/data/sdssdesi_passC_v3/empca/galaxy_empca_v3pre_G.fits')
LINE_KMS = 1000.0
CKMS = 299792.458
SURVEY_KEY = {'desi': 'desi', 'dr8': 'dr8', 'dr17': 'dr17'}


class Telluric:
    def __init__(self, survey, rwave, basis=BASIS):
        self.survey = SURVEY_KEY[survey]
        with fitsio.FITS(basis) as f:
            P = f['EIGVEC'].read().astype(np.float64); w = f['WAVE'].read(); ok = f['PIXMASK'].read().astype(bool)
        assert len(w) == len(rwave) and np.allclose(w, rwave, rtol=1e-6), 'basis grid differs from the rest grid'
        self.P = P; self.ok = ok; self.rwave = rwave; self.nvec = len(P)
        pb = importlib.import_module('02_passB_measure')
        lam = np.array([x[1] for g in pb.LINES.values() for x in g])
        lw = np.log10(rwave)[:, None]; self.linemask = (np.abs(lw - np.log10(lam)[None, :]) < np.log10(1 + LINE_KMS / CKMS)).any(1)
        self.ta = importlib.import_module('tellapply')      # Stars thread: factor(survey, wave 2-D, flux, ivar, model) -> C, a
        self.tb = self.ta.basis(self.survey)

    def model(self, rflux, rivar):
        """weighted least squares of each spectrum on the 11 vectors (basis pixels, lines masked); returns model (NaN off the
        basis), chi2 per fitted pixel."""
        W = np.where(self.ok[None, :] & ~self.linemask[None, :] & np.isfinite(rflux), rivar, 0.0).astype(np.float64)
        F = np.where(W > 0, rflux, 0.0).astype(np.float64); P = self.P
        A = np.einsum('kp,np,lp->nkl', P, W, P); b = np.einsum('kp,np->nk', P, W * F)
        A[:, np.arange(self.nvec), np.arange(self.nvec)] += 1e-10 * np.trace(A, axis1=1, axis2=2)[:, None]
        c = np.linalg.solve(A, b[..., None])[..., 0]
        m = c @ P; m[:, ~self.ok] = np.nan
        npix = (W > 0).sum(1); chi2 = (W * (F - np.nan_to_num(m))**2).sum(1) / np.maximum(npix - self.nvec, 1)
        return m, chi2, npix

    def correct(self, rflux, rivar, z):
        """rflux, rivar (n, NREST) on the rest grid -> corrected copies and a dict of record columns."""
        m, chi2, npix = self.model(rflux, rivar)
        wobs = self.rwave[None, :] * (1.0 + np.asarray(z, float))[:, None]
        iv = np.where(self.linemask[None, :], 0.0, rivar)
        C, a = self.ta.factor(self.survey, wobs, rflux, iv, m)
        C = np.where(np.isfinite(C) & (C > 0), C, 1.0).astype(np.float32)
        rec = {f'TELL_A{k}': a[:, k].astype(np.float32) for k in range(a.shape[1])}
        rec['TELL_CHI2R'] = chi2.astype(np.float32); rec['TELL_NFIT'] = npix.astype(np.int32)
        return (rflux / C).astype(np.float32), (rivar * C**2).astype(np.float32), rec


if __name__ == '__main__':
    import argparse, pandas as pd
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', default='desi', choices=list(SURVEY_KEY)); ap.add_argument('--n', type=int, default=500)
    ap.add_argument('--file', default='bright/hpx7280', help='v2 rest-frame store file (DESI: prog/hpxN)')
    a = ap.parse_args()
    pb = importlib.import_module('02_passB_measure'); T = Telluric(a.survey, pb.RWAVE)
    zz = np.load(f'/Volumes/exdisk3/sdssdesi_restframe_v2/{a.file}.npz')
    tab = pd.read_parquet(os.path.expanduser(f'~/data/sdssdesi_passC_v2/tables/{a.file}.parquet'), columns=['TARGETID', 'Z']).set_index('TARGETID')
    ids, st, off, rid0 = zz['TARGETID'][:a.n], zz['START'], zz['OFFSET'], int(zz['RID_MIN'])
    full = np.zeros((len(ids), pb.NREST), np.float32); fiv = np.zeros_like(full)
    for i in range(len(ids)):
        j0 = rid0 + int(st[i]) - pb.RID_MIN; fl = zz['FLUX'][off[i]:off[i + 1]]; full[i, j0:j0 + len(fl)] = fl; fiv[i, j0:j0 + len(fl)] = zz['IVAR'][off[i]:off[i + 1]]
    z = tab.Z.reindex(ids).to_numpy(float)
    import time; t = time.time(); fc, ivc, rec = T.correct(full, fiv, z); dt = time.time() - t
    print(f'{len(ids)} spectra corrected in {dt:.2f}s ({dt / len(ids) * 1e3:.2f} ms each)')
    for k, v in rec.items():
        print(f'{k}: median {np.nanmedian(v):+.4f}, 16-84 % {np.nanpercentile(v, 16):+.4f}..{np.nanpercentile(v, 84):+.4f}')
    with np.errstate(all='ignore'):
        r = fc / full - 1
    print('median over spectra of max |1/C - 1|:', np.nanmedian(np.nanmax(np.abs(np.where(fiv > 0, r, np.nan)), 1)))
