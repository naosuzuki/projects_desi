#!/usr/bin/env python3
"""Pass B: per-spectrum measurements over the DESI DR2 loa healpix coadds.

Recipe: ms.tex Sec. 3.9 (SDSS-DESI paper) / programs_templates/README.md.
Per spectrum:
  Step 1  spectrophotometry: spectrum x MEAN_PSF_TO_FIBER_SPECFLUX vs FIBERFLUX_G/R
          (DECam or BASS/MzLS by PHOTSYS); scale correction per object (the two-band tilt is
          estimated and stored as CAL_B_EST but not applied unless --per-object-tilt)
  Step 2  Milky Way extinction: Fitzpatrick99, R_V=3.1, E(B-V) = 0.86 * EBV(SFD)
  Step 3  rest frame: flux-conserving rebin to log10 grid (dlog=1e-4, same pixel-id convention
          as sdsslib astrolib), integer-pixel shift by round(log10(1+z)/dlog), F_lam*(1+z)
  Step 4  synthetic rest-frame magnitudes (DECam r, g; SDSS u; NUV tophat 2300-3000) where the
          band is covered; K-corrections; absolute magnitudes with flat LCDM h=0.7, Om=0.3
          from the whole-galaxy FLUX_R (dereddened)
  Step 5  S/N (per camera, total, SDSS pixel scale); emission-line fluxes/EWs (linear LSQ,
          Gaussians at 3 widths, local linear continuum; [OII] doublet sum); Dn4000; fiber fraction
Output: one Parquet per coadd file in OUTDIR/<program>/hpx<N>.parquet (resumable: existing
files are skipped). All spectra of the file are written (stars are needed for the Step 1 check).
Memory: spectra are processed in chunks of CHUNK rows so a worker stays at a few GB whatever the
file size. Run workers with OPENBLAS_NUM_THREADS=1 (the work is elementwise numpy, not BLAS).

Usage:
  python3 02_passB_measure.py --nfiles 20 --seed 1            # random test subset
  python3 02_passB_measure.py --all --shard 0/2                # production, worker 0 of 2 (resumable)
  python3 02_passB_measure.py --files path1 path2 ...
"""
import os, sys, time, argparse, warnings, gc
import numpy as np
import pandas as pd
import fitsio
import scipy.sparse as sp
import extinction
from astropy.cosmology import FlatLambdaCDM
import speclite.filters as sfilters

warnings.filterwarnings('ignore')
MANIFEST = os.path.expanduser('~/github/projects_desi/programs_setups/manifests/healpix_main_darkbright.batch')
ZLOOKUP_COMPACT = os.path.expanduser('~/data/sdssdesi_products/zlookup_compact.npz')   # from 00_build_zlookup.py
OUTDIR = os.path.expanduser('~/data/sdssdesi_passB')
DLOG = 1e-4
REST_MIN, REST_MAX = 1400.0, 9900.0
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
LS_EXT = {'G': 3.214, 'R': 2.165, 'Z': 1.211}      # A / E(B-V)_SFD for LS photometry
FLUX_UNIT = 1e-17                                   # DESI flux unit, erg/s/cm2/A
AB_ZP = 3631e-23 * 2.99792458e18                    # erg/s/cm2/A * A^2 for F_nu = 3631 Jy
TILT_MAX_ERR = 0.03                                 # mag; both synthetic g and r must be this precise for a tilt estimate
PER_OBJECT_TILT = False                             # apply the per-object tilt (mode 2)? default: scale only
CHUNK = 500                                         # spectra per processing chunk (memory cap)
F32 = np.float32

# ---------------------------------------------------------------- filters (vectorized synthetic photometry)
_F = {}
def filt(name):
    if name not in _F:
        if name == 'nuv':
            w = np.arange(2200.0, 3100.1, 1.0); r = ((w >= 2300) & (w <= 3000)).astype(float)
            _F[name] = sfilters.FilterResponse(w, r, dict(group_name='tophat', band_name='nuv'))
        else:
            _F[name] = sfilters.load_filter(name)
    return _F[name]

def band_range(f):
    """wavelength range holding 99% of the response integral (0.5%-99.5%); speclite's DECam curves
    carry low-level tails out to ~11000 A that must not count as 'the band'."""
    r, w = f.response, f.wavelength
    c = np.cumsum(r * np.gradient(w)); c = c / c[-1]
    return w[np.searchsorted(c, 0.005)], w[min(np.searchsorted(c, 0.995), len(w) - 1)]

_BW = {}
def band_weights(f, wave):
    """Per-pixel weights so that  mag = -2.5 log10( sum(flux*wt) / den )  for flux in erg/s/cm2/A."""
    key = (f.name, len(wave), wave[0], wave[-1])
    if key not in _BW:
        r = np.interp(wave, f.wavelength, f.response, left=0, right=0)
        dw = np.gradient(wave)
        wt = r * wave * dw                                    # photon counting
        den = np.sum(AB_ZP / wave**2 * r * wave * dw)
        lo, hi = band_range(f)
        inband = (wave >= lo) & (wave <= hi)
        _BW[key] = (wt, den, inband)
    return _BW[key]

def ab_mags(f, wave, flux2d, good2d, inrange2d=None, ivar2d=None, min_range=0.98, min_good=0.7):
    """AB magnitudes (and errors) for all rows.
    inrange2d: pixels inside the spectrum's wavelength coverage (NaN if < min_range of the band);
    good2d: unmasked pixels; masked pixels inside the band are assigned the band-mean flux density
    (the integral is scaled by total weight / good weight), NaN if < min_good of the band is good."""
    wt, den, inband = band_weights(f, wave)
    n = flux2d.shape[0]
    if inrange2d is None:
        inrange2d = np.ones_like(good2d, dtype=bool)
    if not inband.any():
        return np.full(n, np.nan), np.full(n, np.nan)
    frac_range = inrange2d[:, inband].mean(axis=1)
    frac_good = good2d[:, inband].mean(axis=1)
    g = good2d & inrange2d
    wsum_good = g.astype(float) @ wt
    num_good = np.where(g, np.nan_to_num(flux2d), 0.0) @ wt
    with np.errstate(all='ignore'):
        num = num_good * (wt.sum() / np.where(wsum_good > 0, wsum_good, np.nan))
        m = -2.5 * np.log10(num / den)
        ok = (frac_range >= min_range) & (frac_good >= min_good) & (num > 0)
        if ivar2d is not None:
            var_good = np.where(g & (ivar2d > 0), 1.0 / np.where(ivar2d > 0, ivar2d, 1), 0.0) @ wt**2
            err = 1.0857 * np.sqrt(var_good) / np.where(num_good > 0, num_good, np.nan)
        else:
            err = np.full(n, np.nan)
    return np.where(ok, m, np.nan), np.where(ok, err, np.nan)

# ---------------------------------------------------------------- rebinning (sparse, flux-conserving)
_REBIN = {}
def rebin_matrix(wave):
    """Flux-conserving rebin of a fixed observed grid onto the log10 grid (centers 10^(rid*DLOG)).
    Returns (rid, M, M^2) with rebinned = M @ pixelvalues; rows of M sum to 1 (fractional overlaps)."""
    key = (len(wave), wave[0], wave[-1])
    if key in _REBIN:
        return _REBIN[key]
    edges = np.empty(len(wave) + 1)
    edges[1:-1] = np.sqrt(wave[:-1] * wave[1:])         # geometric midpoints, as in astrolib
    edges[0] = np.sqrt(wave[0]**3 / wave[1]); edges[-1] = np.sqrt(wave[-1]**3 / wave[-2])
    hmin = int(np.ceil((np.log10(edges[0]) / DLOG) + 0.5))
    hmax = int(np.floor((np.log10(edges[-1]) / DLOG) - 0.5))
    rid = np.arange(hmin, hmax + 1)
    bedges = 10.0**((rid - 0.5) * DLOG); redges = 10.0**((rid + 0.5) * DLOG)
    rows, cols, vals = [], [], []
    for k in range(len(rid)):
        h0 = np.searchsorted(edges, bedges[k], side='right') - 1
        h1 = np.searchsorted(edges, redges[k], side='left')
        for h in range(max(h0, 0), min(h1, len(wave))):
            ov = min(edges[h+1], redges[k]) - max(edges[h], bedges[k])
            if ov > 0:
                rows.append(k); cols.append(h); vals.append(ov / (redges[k] - bedges[k]))
    M = sp.csr_matrix((vals, (rows, cols)), shape=(len(rid), len(wave)))
    _REBIN[key] = (rid, M, M.multiply(M).tocsr())
    return _REBIN[key]

def rebin(wave, flux, ivar, mask):
    """flux [n,npix], ivar, mask -> (rid, rflux, rivar) on the log grid; ivar=0 where <50% good coverage."""
    rid, M, M2 = rebin_matrix(wave)
    good = (ivar > 0) & (mask == 0)
    f = np.where(good, flux, 0.0); var = np.where(good, 1.0 / np.where(ivar > 0, ivar, 1), 0.0)
    wsum = (M @ good.T.astype(float)).T                   # fraction of each log bin covered by good pixels
    rflux = (M @ f.T).T / np.where(wsum > 0, wsum, 1)
    rvar = (M2 @ var.T).T / np.where(wsum > 0, wsum, 1)**2
    rivar = np.where((wsum > 0.5) & (rvar > 0), 1.0 / np.where(rvar > 0, rvar, 1), 0.0)
    return rid, rflux, rivar

# ---------------------------------------------------------------- camera combination
def combine_cameras(f):
    """Read B/R/Z and combine onto one 0.8 A grid (ivar-weighted in the overlaps). float32 arrays."""
    waves = {c: f[f'{c}_WAVELENGTH'].read() for c in 'BRZ'}
    w = np.arange(3600.0, 9824.0 + 0.4, 0.8)
    n = f['B_FLUX'].get_dims()[0]
    num = np.zeros((n, len(w)), F32); den = np.zeros((n, len(w)), F32); anygood = np.zeros((n, len(w)), bool)
    percam = {}
    for c in 'BRZ':
        fl, iv, mk = f[f'{c}_FLUX'].read(), f[f'{c}_IVAR'].read(), f[f'{c}_MASK'].read()
        good = (iv > 0) & (mk == 0)
        i0 = int(round((waves[c][0] - 3600.0) / 0.8)); sl = slice(i0, i0 + len(waves[c]))
        num[:, sl] += np.where(good, fl * iv, 0).astype(F32); den[:, sl] += np.where(good, iv, 0).astype(F32); anygood[:, sl] |= good
        with np.errstate(all='ignore'):
            percam[c] = np.nanmedian(np.where(good, fl * np.sqrt(iv), np.nan), axis=1)
        del fl, iv, mk, good
    flux = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0).astype(F32)
    mask = (~anygood).astype(np.int8)
    return w, flux, den, mask, percam

# ---------------------------------------------------------------- emission lines
LINES = {   # group: list of (name, rest wavelength)
    'oii':  [('OII_3726', 3726.03), ('OII_3729', 3728.82)],
    'neiii':[('NEIII_3869', 3868.76)],
    'hd':   [('HDELTA', 4101.74)],
    'hg':   [('HGAMMA', 4340.47)],
    'hb':   [('HBETA', 4861.33)],
    'oiii': [('OIII_4959', 4958.91), ('OIII_5007', 5006.84)],
    'ha':   [('NII_6548', 6548.05), ('HALPHA', 6562.80), ('NII_6583', 6583.45)],
    'sii':  [('SII_6716', 6716.44), ('SII_6731', 6730.82)],
}
SIGMAS_KMS = [80.0, 150.0, 300.0]
SIDEBAND = 30.0   # A on each side of the group for the linear continuum
CKMS = 299792.458

def fit_lines(rwave, rflux, rivar):
    """Vectorized over spectra. rflux [n, npix] on common rest grid rwave. Returns dict of arrays."""
    n = rflux.shape[0]; out = {}
    for g, lines in LINES.items():
        lam = np.array([l[1] for l in lines]); lo, hi = lam.min(), lam.max()
        core = 3 * lam.max() * SIGMAS_KMS[-1] / CKMS
        sel = (rwave >= lo - core - SIDEBAND) & (rwave <= hi + core + SIDEBAND)
        w = rwave[sel]; y = rflux[:, sel].astype(float); wt = rivar[:, sel].astype(float)
        x = (w - w.mean()) / 100.0
        best_chi2 = np.full(n, np.inf); best_coef = None; best_cov = None; best_sig = np.full(n, np.nan)
        p = 2 + len(lam)
        for sig in SIGMAS_KMS:
            cols = [np.ones_like(w), x] + [np.exp(-0.5 * ((w - l) / (l * sig / CKMS))**2) for l in lam]
            A = np.stack(cols, axis=1)                                # [npix, p]
            AtW = A.T[None, :, :] * wt[:, None, :]                    # [n, p, npix]
            N = AtW @ A + 1e-9 * np.eye(p)[None]; bvec = np.einsum('npk,nk->np', AtW, y)
            coef = np.linalg.solve(N, bvec[..., None])[..., 0]
            cov = np.linalg.inv(N)
            model = coef @ A.T
            chi2 = np.sum(wt * (y - model)**2, axis=1)
            better = chi2 < best_chi2
            if best_coef is None:
                best_coef, best_cov = coef.copy(), cov.copy()
            else:
                best_coef[better] = coef[better]; best_cov[better] = cov[better]
            best_chi2[better] = chi2[better]; best_sig[better] = sig
        ok = np.sum(wt > 0, axis=1) > 10
        sig = np.where(np.isnan(best_sig), SIGMAS_KMS[0], best_sig)
        for j, (name, l) in enumerate(lines):
            amp = best_coef[:, 2 + j]; amp_err = np.sqrt(np.abs(best_cov[:, 2 + j, 2 + j]))
            width = l * sig / CKMS
            flux = amp * width * np.sqrt(2 * np.pi); flux_err = amp_err * width * np.sqrt(2 * np.pi)
            cont = best_coef[:, 0] + best_coef[:, 1] * (l - w.mean()) / 100.0
            with np.errstate(all='ignore'):
                ew = np.where(cont > 0, flux / cont, np.nan); ew_err = np.where(cont > 0, flux_err / np.abs(cont), np.nan)
            out[f'{name}_FLUX'] = np.where(ok, flux, np.nan); out[f'{name}_FLUX_ERR'] = np.where(ok, flux_err, np.nan)
            out[f'{name}_EW'] = np.where(ok, ew, np.nan); out[f'{name}_EW_ERR'] = np.where(ok, ew_err, np.nan)
        out[f'{g}_SIGMA'] = np.where(ok, best_sig, np.nan)
        if g == 'oii':
            # the doublet (2.8 A = 225 km/s) is degenerate at sigma >= 150 km/s: report the sum as the robust quantity
            width = lam * sig[:, None] / CKMS
            fl = best_coef[:, 2:4] * width * np.sqrt(2 * np.pi)
            var = (best_cov[:, 2, 2] * width[:, 0]**2 + best_cov[:, 3, 3] * width[:, 1]**2
                   + 2 * best_cov[:, 2, 3] * width[:, 0] * width[:, 1]) * 2 * np.pi
            lsum = 3727.42
            cont = best_coef[:, 0] + best_coef[:, 1] * (lsum - w.mean()) / 100.0
            with np.errstate(all='ignore'):
                out['OII_SUM_FLUX'] = np.where(ok, fl.sum(axis=1), np.nan)
                out['OII_SUM_FLUX_ERR'] = np.where(ok, np.sqrt(np.abs(var)), np.nan)
                out['OII_SUM_EW'] = np.where(ok & (cont > 0), fl.sum(axis=1) / cont, np.nan)
                out['OII_SUM_EW_ERR'] = np.where(ok & (cont > 0), np.sqrt(np.abs(var)) / np.abs(cont), np.nan)
    return out

def dn4000(rwave, rflux, rivar):
    fnu = rflux.astype(float) * rwave**2; good = rivar > 0
    def band(lo, hi):
        s = (rwave >= lo) & (rwave < hi)
        g = good[:, s]; n = g.sum(axis=1); nn = np.where(n > 0, n, 1)
        m = np.where(n > 5, np.sum(np.where(g, fnu[:, s], 0), axis=1) / nn, np.nan)
        var = np.where(n > 5, np.sum(np.where(g, rwave[s]**4 / np.where(rivar[:, s] > 0, rivar[:, s], 1), 0), axis=1) / nn**2, np.nan)
        return m, var
    b, vb = band(3850, 3950); r, vr = band(4000, 4100)
    with np.errstate(all='ignore'):
        d = r / b; err = d * np.sqrt(vr / r**2 + vb / b**2)
    return d, err

# ---------------------------------------------------------------- per-chunk measurements
RID_MIN = int(round(np.log10(REST_MIN) / DLOG)); RID_MAX = int(round(np.log10(REST_MAX) / DLOG))
NREST = RID_MAX - RID_MIN + 1
RWAVE = 10.0**(np.arange(RID_MIN, RID_MAX + 1) * DLOG)
A1_F99 = None   # F99 A(lambda) for A_V=1 on the observed grid, set on first use

def measure_chunk(df, w, flux, ivar, mask):
    """df: chunk of the fibermap+redshift table (n rows); flux/ivar/mask [n, npix] on grid w."""
    global A1_F99
    n = len(df)
    good = (ivar > 0) & (mask == 0)

    # --- Step 5a: S/N in the observed frame
    with np.errstate(all='ignore'):
        snr = np.where(good, flux * np.sqrt(ivar), np.nan)
        df['SNR_ALL'] = np.nanmedian(snr, axis=1)
        df['SNR_SDSSPIX'] = np.nanmedian(snr * np.sqrt(w * np.log(10) * DLOG / 0.8).astype(F32), axis=1)
        df['NGOOD'] = good.sum(axis=1)
    del snr

    # --- Step 1: spectrophotometry against fiber flux
    psf2fib = df['MEAN_PSF_TO_FIBER_SPECFLUX'].to_numpy(dtype=float)
    psf2fib = np.where(np.isfinite(psf2fib) & (psf2fib > 0), psf2fib, np.nan)
    scale0 = (np.nan_to_num(psf2fib, nan=1.0) * FLUX_UNIT)[:, None]
    ffib = flux * scale0                                    # fiber-scale spectrum, cgs (float64 via broadcast)
    ivar_fib = ivar / scale0**2
    south = df['PHOTSYS'].to_numpy() == 'S'
    names = {True: {'g': 'decam2014-g', 'r': 'decam2014-r', 'z': 'decam2014-z'},
             False: {'g': 'BASS-g', 'r': 'BASS-r', 'z': 'MzLS-z'}}
    zf_hi = max(band_range(filt('decam2014-z'))[1], band_range(filt('MzLS-z'))[1])
    wpad = np.arange(w[-1] + 0.8, zf_hi + 60, 0.8); wz = np.concatenate([w, wpad])
    def syn_mags(spec, ivar_s, covered):
        out, err = {}, {}
        padval = np.nanmean(np.where(covered[:, -250:], spec[:, -250:], np.nan), axis=1)
        specz = np.concatenate([spec, np.repeat(np.nan_to_num(padval)[:, None], len(wpad), axis=1)], axis=1)
        covz = np.concatenate([covered, np.ones((n, len(wpad)), bool)], axis=1)
        for b in 'grz':
            m = np.full(n, np.nan); e = np.full(n, np.nan)
            for s_flag in (True, False):
                rows = south == s_flag
                if rows.any():
                    if b == 'z':
                        m[rows], _ = ab_mags(filt(names[s_flag][b]), wz, specz[rows], covz[rows])
                    else:
                        m[rows], e[rows] = ab_mags(filt(names[s_flag][b]), w, spec[rows], covered[rows], ivar2d=ivar_s[rows])
            out[b] = m; err[b] = e
        return out, err
    syn, syn_err = syn_mags(ffib, ivar_fib, good)
    for b in 'grz':
        df[f'SYN_FIBER_{b.upper()}_RAW'] = np.where(np.isfinite(psf2fib), syn[b], np.nan)
    df['SYN_FIBER_G_ERR'], df['SYN_FIBER_R_ERR'] = syn_err['g'], syn_err['r']
    with np.errstate(all='ignore'):
        cat = {b: 22.5 - 2.5 * np.log10(df[f'FIBERFLUX_{b.upper()}'].to_numpy(dtype=float)) for b in 'grz'}
        dg = cat['g'] - df['SYN_FIBER_G_RAW'].to_numpy(); dr = cat['r'] - df['SYN_FIBER_R_RAW'].to_numpy()
    leff = {s_flag: (filt(names[s_flag]['g']).effective_wavelength.value, filt(names[s_flag]['r']).effective_wavelength.value) for s_flag in (True, False)}
    xg = np.where(south, np.log10(leff[True][0] / 6000), np.log10(leff[False][0] / 6000))
    xr = np.where(south, np.log10(leff[True][1] / 6000), np.log10(leff[False][1] / 6000))
    eg, er = syn_err['g'], syn_err['r']
    with np.errstate(all='ignore'):
        Rg, Rr = 10**(-0.4 * dg), 10**(-0.4 * dr)                       # catalog/synthetic flux ratios
        precise = np.isfinite(eg) & np.isfinite(er) & (eg < TILT_MAX_ERR) & (er < TILT_MAX_ERR)
        b2 = (Rr - Rg) / (xr - xg); a2 = Rr - b2 * xr
        wg = 1.0 / np.where(np.isfinite(eg) & (eg > 0), eg, 0.3)**2; wr = 1.0 / np.where(np.isfinite(er) & (er > 0), er, 0.3)**2
        wg = np.where(np.isfinite(Rg), wg, 0); wr = np.where(np.isfinite(Rr), wr, 0)
        a1 = (np.nan_to_num(Rg) * wg + np.nan_to_num(Rr) * wr) / np.where(wg + wr > 0, wg + wr, np.nan)
    if PER_OBJECT_TILT:
        a_t = np.where(precise, a2, a1); b_t = np.where(precise, b2, 0.0)
    else:
        a_t = a1; b_t = np.zeros(n)
    calok = np.isfinite(a_t) & (a_t > 0) & (np.abs(np.nan_to_num(dg, nan=0)) < 0.5) & (np.abs(np.nan_to_num(dr, nan=0)) < 0.5) & (np.isfinite(dg) | np.isfinite(dr))
    df['CAL_A'], df['CAL_B'] = np.where(calok, a_t, np.nan), np.where(calok, b_t, np.nan)
    df['CAL_B_EST'] = np.where(precise, b2, np.nan)
    df['CAL_MODE'] = np.where(calok, np.where(PER_OBJECT_TILT & precise, 2, 1), 0).astype(np.int8)
    df['CAL_DG'], df['CAL_DR'] = dg, dr
    df['CAL_OK'] = calok
    x = np.log10(w / 6000.0)
    corr = np.where(calok[:, None], np.nan_to_num(a_t, nan=1.0)[:, None] + np.nan_to_num(b_t)[:, None] * x[None, :], 1.0)
    corr = np.where(corr > 0, corr, 1.0)
    fcal = ffib * corr; icar = ivar_fib / corr**2
    icar = np.where(good & np.isfinite(icar), icar, 0.0)
    del ffib, ivar_fib
    synC, _ = syn_mags(fcal, icar, good)
    for b in 'grz':
        df[f'SYN_FIBER_{b.upper()}_CAL'] = np.where(calok, synC[b], np.nan)
    df['CAL_DZ_CHECK'] = cat['z'] - df['SYN_FIBER_Z_CAL']

    # --- Step 2: Milky Way extinction (Fitzpatrick99, R_V=3.1, E(B-V)=0.86*EBV)
    ebv = np.nan_to_num(df['EBV'].to_numpy(dtype=float), nan=0.0)
    if A1_F99 is None or len(A1_F99) != len(w):
        A1_F99 = extinction.fitzpatrick99(w, 1.0, 3.1)                 # A(lam) for A_V = 1 (linear in A_V)
    A = (3.1 * 0.86 * ebv)[:, None] * A1_F99[None, :]
    fder = fcal * 10**(0.4 * A); ider = icar * 10**(-0.8 * A)
    del fcal, icar, corr, A
    m_obs_r = np.full(n, np.nan)
    for s_flag in (True, False):
        rows = south == s_flag
        if rows.any():
            m_obs_r[rows], _ = ab_mags(filt(names[s_flag]['r']), w, fder[rows], good[rows])
    df['SYN_R_OBS_DERED'] = np.where(calok, m_obs_r, np.nan)

    # --- Step 3: rest frame
    rid_obs, rflux_o, rivar_o = rebin(w, fder, ider, mask)
    del fder, ider
    z = df['Z'].to_numpy(dtype=float); zok = np.isfinite(z) & (z > -0.005)
    rflux = np.zeros((n, NREST), F32); rivar = np.zeros((n, NREST), F32); rinrange = np.zeros((n, NREST), bool)
    shift = np.where(zok, np.round(np.log10(1 + np.clip(z, -0.004, None)) / DLOG), 0).astype(int)
    for i in np.where(zok)[0]:
        r0 = rid_obs - shift[i] - RID_MIN                              # rest-grid indices
        ok = (r0 >= 0) & (r0 < NREST)
        rflux[i, r0[ok]] = rflux_o[i, ok] * (1 + z[i]); rivar[i, r0[ok]] = rivar_o[i, ok] / (1 + z[i])**2
        rinrange[i, r0[ok]] = True
    del rflux_o, rivar_o
    rcov = rivar > 0

    # --- Step 4: rest-frame synthetic magnitudes, K, absolute magnitudes
    DM = np.where(zok & (z > 0), COSMO.distmod(np.clip(z, 1e-4, None)).value, np.nan)
    df['DM'] = DM
    bands = {'R': 'decam2014-r', 'G': 'decam2014-g', 'U': 'sdss2010-u', 'NUV': 'nuv'}
    with np.errstate(all='ignore'):
        m_r_tot = 22.5 - 2.5 * np.log10(df['FLUX_R'].to_numpy(dtype=float)) - LS_EXT['R'] * ebv   # whole galaxy, dereddened
        df['M_R_TOTAL_DERED'] = m_r_tot
        df['FIBER_FRACTION_R'] = df['FIBERFLUX_R'] / df['FLUX_R']
        for b, name in bands.items():
            mr, _ = ab_mags(filt(name), RWAVE, rflux, rcov, inrange2d=rinrange)
            mr = np.where(zok & calok, mr, np.nan)
            df[f'MREST_{b}_FIBER'] = mr
            df[f'K_{b}'] = df['SYN_R_OBS_DERED'] - mr                 # K = m_obs(r) - m_rest(band), fiber spectrum
            df[f'MABS_{b}'] = m_r_tot - DM - df[f'K_{b}']            # whole-galaxy absolute magnitude in band

    # --- Step 5b: lines and Dn4000 in the rest frame
    for k, v in fit_lines(RWAVE, rflux, rivar).items():
        df[k] = v
    d, de = dn4000(RWAVE, rflux, rivar)
    df['DN4000'], df['DN4000_ERR'] = d, de
    return df

def process_file(path, zl, program):
    t0 = time.time()
    f = fitsio.FITS(path)
    fm = f['FIBERMAP'].read(columns=['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'EBV', 'PHOTSYS',
                                      'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
                                      'MEAN_PSF_TO_FIBER_SPECFLUX', 'COADD_FIBERSTATUS'])
    w, flux, ivar, mask, percam = combine_cameras(f)
    f.close()
    n = len(fm)
    df = pd.DataFrame({c: (np.char.strip(fm[c].astype(str)) if fm[c].dtype.kind in 'SU' else fm[c].astype(fm[c].dtype.newbyteorder('='))) for c in fm.dtype.names})
    df['PROGRAM'] = program
    df = df.merge(zl, on=['TARGETID', 'PROGRAM'], how='left')
    assert len(df) == n
    df['SNR_B'], df['SNR_R'], df['SNR_Z'] = percam['B'], percam['R'], percam['Z']
    t_read = time.time() - t0
    parts = []
    for i0 in range(0, n, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, n))
        parts.append(measure_chunk(df.iloc[sl].copy(), w, flux[sl], ivar[sl], mask[sl]))
        gc.collect()
    out = pd.concat(parts, ignore_index=True)
    out['FILE'] = os.path.basename(path)
    out['T_READ'] = t_read
    out['PROC_SEC'] = time.time() - t0
    return out

# ---------------------------------------------------------------- driver
OVERFLOW_ROOTS = ['/Volumes/exdisk1/', '/Volumes/My Book/']   # the exFAT disk mounts under either name

def local_path(p):
    """manifest path (exdisk2) -> existing local path; 969 overflow files live on the exFAT disk."""
    if os.path.exists(p):
        return p
    for root in OVERFLOW_ROOTS:
        alt = p.replace('/Volumes/exdisk2/', root)
        if os.path.exists(alt):
            return alt
    return None

def main():
    global PER_OBJECT_TILT
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default=MANIFEST)
    ap.add_argument('--nfiles', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--files', nargs='*')
    ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--shard', default='0/1', help='i/n: process every n-th file starting at i (parallel workers)')
    ap.add_argument('--per-object-tilt', action='store_true', help='apply the two-band tilt per object (default: scale only)')
    a = ap.parse_args()
    PER_OBJECT_TILT = a.per_object_tilt
    if a.files:
        paths = a.files
    else:
        paths = [l.split()[1] for l in open(a.manifest)]
        if not a.all:
            rng = np.random.default_rng(a.seed); paths = [paths[i] for i in rng.choice(len(paths), a.nfiles, replace=False)]
    si, sn = (int(v) for v in a.shard.split('/'))
    paths = paths[si::sn]
    t0 = time.time()
    # compact lookup (00_build_zlookup.py): plain numpy arrays sorted by (program, healpix), no conversion transient
    zl = dict(np.load(ZLOOKUP_COMPACT))
    offsets = {int(k): (int(s), int(e)) for k, s, e in zip(zl['OFF_KEY'], zl['OFF_START'], zl['OFF_STOP'])}
    print(f'zlookup loaded: {len(zl["TARGETID"]):,} rows, {sum(v.nbytes for v in zl.values())/1e9:.1f} GB, {time.time()-t0:.0f}s', flush=True)
    STYPE = np.array(['GALAXY', 'QSO', 'STAR', 'OTHER'])
    done = 0; missing = []
    for k, p in enumerate(paths):
        prog = 'dark' if '/dark/' in p else 'bright'
        hpx = os.path.basename(p).split('-')[-1].replace('.fits', '')
        out = os.path.join(a.outdir, prog, f'hpx{hpx}.parquet')
        if os.path.exists(out) and not a.overwrite:
            continue
        # if the source disk has dropped (USB disconnect), wait for it instead of walking the manifest as MISSING
        waited = 0
        while p.startswith('/Volumes/exdisk2/') and not os.path.ismount('/Volumes/exdisk2'):
            if waited == 0:
                print(f'DISK UNMOUNTED /Volumes/exdisk2 at {time.strftime("%H:%M:%S")} - waiting for remount', flush=True)
            time.sleep(60); waited += 60
        if waited:
            print(f'DISK BACK after {waited//60} min at {time.strftime("%H:%M:%S")}', flush=True)
        lp = local_path(p)
        if lp is None:
            print('MISSING', p, flush=True); missing.append(p); continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        try:
            s, e = offsets.get((0 if prog == 'dark' else 1) * 10**6 + int(hpx), (0, 0))
            zsub = pd.DataFrame({'TARGETID': zl['TARGETID'][s:e], 'Z': zl['Z'][s:e].astype(float), 'ZERR': zl['ZERR'][s:e].astype(float),
                                 'ZWARN': zl['ZWARN'][s:e], 'SPECTYPE': STYPE[zl['SPECTYPE_CODE'][s:e]], 'DELTACHI2': zl['DELTACHI2'][s:e].astype(float),
                                 'ZCAT_PRIMARY': zl['ZCAT_PRIMARY'][s:e]})
            zsub['PROGRAM'] = prog
            df = process_file(lp, zsub, prog)
            df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
            done += 1
            print(f'{k+1}/{len(paths)} {os.path.basename(lp)} n={len(df)} read {df.T_READ.iloc[0]:.1f}s total {df.PROC_SEC.iloc[0]:.1f}s  elapsed {time.time()-t0:.0f}s', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print('ERROR', p, repr(e), flush=True)
    print(f'done: {done} files, missing {len(missing)}, {time.time()-t0:.0f}s')
    if missing:
        with open(os.path.join(a.outdir, f'missing_shard{si}of{sn}.txt'), 'w') as fh:
            fh.write('\n'.join(missing) + '\n')

if __name__ == '__main__':
    main()
