#!/usr/bin/env python3
"""qsolib: quasar per-spectrum measurements shared by the SDSS and DESI drivers (programs_quasars/01_, 02_).

Imports pass B (programs_templates/02_passB_measure.py) for the common pieces (filters, synthetic magnitudes, rebin,
extinction constants, galaxy line list, Dn4000) and adds the quasar measurements on a 600-9900 A rest grid:
  power-law continuum slope (ALPHA_LAM, line-free windows, chi2 in linear flux), broad and oxygen lines (flux, EW,
  sigma, centroid VOFF and peak VPEAK velocities), narrow [O III]/[O I], forest coverage/transmission, AI/BI troughs
  for C IV / Si IV / Mg II / Al III and BAL_CLASS, line redshifts and the systemic redshift Z_SYS (peak of [O III] ->
  Mg II -> O I 1304 -> pipeline z).
Entry points:
  quasar_rest(df, w, flux, ivar, mask)           DESI: linear observed grid -> rest grid (uses CAL_A/CAL_B, EBV from df)
  rest_from_loggrid(rid0, flux, ivar, mask, z, scale, ebv)   SDSS: spPlate log10 grid (dlog=1e-4) -> rest grid
  measure_rest(df, rflux, rivar)                 all quasar measurements + galaxy lines/Dn4000 on the rest spectrum
See docs/HANDOVER_quasars.md sec. 8 for the design decisions (Nao, 2026-10-24).
"""
import os, sys, importlib.util
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('passB', os.path.join(HERE, '..', 'programs_templates', '02_passB_measure.py'))
passB = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(passB)

# SDSS and DESI wavelengths are in vacuum; pass B's optical line list is in AIR wavelengths (found 2026-10-24 from the
# +76 km/s offset of our [O III] peak against DR16Q). Override here with vacuum values (Morton 1991 / NIST) without
# touching 02_passB_measure.py, which is mid-run for the galaxies (that fix is for the galaxy thread to schedule).
VAC = {'OII_3726': 3727.09, 'OII_3729': 3729.88, 'NEIII_3869': 3869.86, 'HDELTA': 4102.89, 'HGAMMA': 4341.68, 'HBETA': 4862.68,
       'OIII_4959': 4960.30, 'OIII_5007': 5008.24, 'NII_6548': 6549.86, 'HALPHA': 6564.61, 'NII_6583': 6585.27,
       'SII_6716': 6718.29, 'SII_6731': 6732.67}
passB.LINES = {g: [(n, VAC.get(n, l)) for n, l in lines] for g, lines in passB.LINES.items()}

DLOG = passB.DLOG
CKMS = passB.CKMS
F32 = np.float32
QREST_MIN, QREST_MAX = 600.0, 9900.0
QRID_MIN = int(round(np.log10(QREST_MIN) / DLOG)); QRID_MAX = int(round(np.log10(QREST_MAX) / DLOG))
QNREST = QRID_MAX - QRID_MIN + 1
QWAVE = 10.0**(np.arange(QRID_MIN, QRID_MAX + 1) * DLOG)

# line-free continuum windows (rest A), after Vanden Berk+01 Table 2 / Shen+11 / Kamble+20
CONT_WINDOWS = [(1275, 1290), (1350, 1365), (1445, 1465), (1690, 1705), (1770, 1800), (2155, 2190), (2225, 2250),
                (2650, 2700), (2950, 3000), (3300, 3350), (3775, 3825), (4200, 4230), (5100, 5130), (5600, 5700),
                (6000, 6100), (6900, 7000), (7400, 7500), (8000, 8100)]
CONT_MIN_PIX = 10          # good pixels per window to count it
CONT_MIN_RATIO = 1.25      # lam_max/lam_min of the windows used must exceed this

# broad-line groups: name -> (list of (line, rest A), sigma grid km/s)
BROAD_SIG = [600.0, 1200.0, 2500.0, 5000.0]
QLINES = {
    'ovi':   ([('OVI', 1033.82)], BROAD_SIG),                 # O VI 1031.9+1037.6 blend (+Lybeta 1025.7 in the blue wing)
    'lya':   ([('LYA', 1215.67), ('NV', 1240.81)], BROAD_SIG),
    'oi':    ([('OI_1304', 1304.35)], BROAD_SIG[:2], 0.012),   # O I 1302+1305 + Si II 1306 blend; narrow sidebands (N V wing, C II 1335)
    'oiii_sf': ([('OIII_1663', 1663.5)], BROAD_SIG[:2], 0.012), # O III] 1660.8+1666.1 semi-forbidden (He II 1640, Al II 1671 nearby)
    'siiv':  ([('SIIV', 1399.8)], BROAD_SIG),
    'civ':   ([('CIV', 1549.06)], BROAD_SIG),
    'ciii':  ([('ALIII', 1857.4), ('CIII', 1908.73)], BROAD_SIG),
    'mgii':  ([('MGII', 2798.75)], BROAD_SIG),
    'hb_b':  ([('HBETA_B', 4862.68)], BROAD_SIG),
    'ha_b':  ([('HALPHA_B', 6564.61)], BROAD_SIG),
    'oiii_n': ([('OIII_N', 5008.24)], [100.0, 200.0, 400.0, 600.0]),
    'oi_n':   ([('OI_6300', 6302.05)], [100.0, 200.0, 400.0, 600.0]),
}
OXYGEN_NOTE = '[O II] 3727 and [O III] 4959/5007 fluxes/EWs come from the pass-B galaxy line list (OII_SUM_*, OIII_5007_*).'
SIDEBAND_FRAC = 0.025      # width of each continuum sideband as a fraction of lambda (~ 7500 km/s); a group may override (3rd item)
LINE_Z = {'Z_MGII': 'MGII', 'Z_OIII': 'OIII_N', 'Z_CIV': 'CIV', 'Z_HBETA': 'HBETA_B'}

FOREST = (1040.0, 1190.0)
BAL_VRANGE = (-25000.0, 0.0)      # km/s, relative to the line (AI); BI uses -25000 < v < -3000 and >= 2000 km/s
BAL_MIN_WIDTH = 450.0             # km/s (AI definition of Hall+02 as used by Trump+06, Lyke+20)
BAL_DEPTH = 0.9
BAL_SMOOTH = 5                    # pixels (~350 km/s) boxcar before the depth test
BAL_NSIG = 2.0                    # smoothed flux must be this many sigma below 0.9 C


# ---------------------------------------------------------------- rest-frame spectrum on the quasar grid
def quasar_rest(df, w, flux, ivar, mask):
    """Re-apply the pass-B calibration (CAL_A, CAL_B from df) and extinction, rebin and shift to the quasar
    rest grid QWAVE. Returns rflux, rivar [n, QNREST] (F_lam * (1+z), erg/s/cm2/A) and the coverage mask."""
    n = len(df)
    good = (ivar > 0) & (mask == 0)
    psf2fib = df['MEAN_PSF_TO_FIBER_SPECFLUX'].to_numpy(dtype=float)
    psf2fib = np.where(np.isfinite(psf2fib) & (psf2fib > 0), psf2fib, 1.0)
    a = np.nan_to_num(df['CAL_A'].to_numpy(dtype=float), nan=1.0); b = np.nan_to_num(df['CAL_B'].to_numpy(dtype=float))
    x = np.log10(w / 6000.0)
    corr = a[:, None] + b[:, None] * x[None, :]
    corr = np.where(corr > 0, corr, 1.0) * (psf2fib * passB.FLUX_UNIT)[:, None]
    ebv = np.nan_to_num(df['EBV'].to_numpy(dtype=float), nan=0.0)
    if passB.A1_F99 is None or len(passB.A1_F99) != len(w):
        passB.A1_F99 = passB.extinction.fitzpatrick99(w, 1.0, 3.1)
    A = (3.1 * 0.86 * ebv)[:, None] * passB.A1_F99[None, :]
    f = flux * corr * 10**(0.4 * A)
    iv = np.where(good, ivar / corr**2 * 10**(-0.8 * A), 0.0)
    del corr, A
    rid_obs, rflux_o, rivar_o = passB.rebin(w, f, iv, mask)
    del f, iv
    z = df['Z'].to_numpy(dtype=float); zok = np.isfinite(z) & (z > -0.005)
    rflux = np.zeros((n, QNREST), F32); rivar = np.zeros((n, QNREST), F32)
    shift = np.where(zok, np.round(np.log10(1 + np.clip(z, -0.004, None)) / DLOG), 0).astype(int)
    for i in np.where(zok)[0]:
        r0 = rid_obs - shift[i] - QRID_MIN
        ok = (r0 >= 0) & (r0 < QNREST)
        rflux[i, r0[ok]] = rflux_o[i, ok] * (1 + z[i]); rivar[i, r0[ok]] = rivar_o[i, ok] / (1 + z[i])**2
    return rflux, rivar


# ---------------------------------------------------------------- power-law continuum
ALPHA_GRID = np.arange(-3.5, 1.5001, 0.02)

def fit_powerlaw(rflux, rivar):
    """F_lam = F3000 (lam/3000)^alpha fitted in LINEAR flux space (chi2 with the pixel ivar) over the line-free
    windows: alpha on a grid, amplitude solved analytically; error from the chi2 curvature. Vectorized."""
    n = rflux.shape[0]
    inwin = np.zeros(len(QWAVE), bool); widx = np.full(len(QWAVE), -1)
    for k, (lo, hi) in enumerate(CONT_WINDOWS):
        s = (QWAVE >= lo) & (QWAVE < hi); inwin |= s; widx[s] = k
    sel = np.where(inwin)[0]
    y = rflux[:, sel].astype(float); iv = rivar[:, sel].astype(float); ws = widx[sel]
    lam = QWAVE[sel] / 3000.0
    use = iv > 0
    nwin = len(CONT_WINDOWS)
    cnt = np.zeros((n, nwin), int)
    for k in range(nwin):
        cnt[:, k] = use[:, ws == k].sum(axis=1)
    winok = cnt >= CONT_MIN_PIX
    use &= winok[:, ws]
    wt = np.where(use, iv, 0.0)
    chi2 = np.full((n, len(ALPHA_GRID)), np.inf); amp = np.zeros((n, len(ALPHA_GRID)))
    wy2 = (wt * y**2).sum(axis=1)
    for j, al in enumerate(ALPHA_GRID):
        A = lam**al
        sAA = (wt * A**2).sum(axis=1); sAy = (wt * A * y).sum(axis=1)
        with np.errstate(all='ignore'):
            a = np.where(sAA > 0, sAy / sAA, 0.0)
        amp[:, j] = a; chi2[:, j] = wy2 - a * sAy          # = sum wt (y - aA)^2 at the optimum
    jb = np.argmin(chi2, axis=1); ii = np.arange(n)
    alpha = ALPHA_GRID[jb]; a0 = amp[ii, jb]; cmin = chi2[ii, jb]
    # curvature from the neighbours (parabola through 3 points): sigma_alpha = sqrt(2 / d2chi2)
    jm = np.clip(jb - 1, 0, len(ALPHA_GRID) - 1); jp = np.clip(jb + 1, 0, len(ALPHA_GRID) - 1)
    d2 = (chi2[ii, jp] - 2 * cmin + chi2[ii, jm]) / (ALPHA_GRID[1] - ALPHA_GRID[0])**2
    with np.errstate(all='ignore'):
        alpha_err = np.where(d2 > 0, np.sqrt(2.0 / d2), np.nan)
    ndof = use.sum(axis=1) - 2
    lam_used = np.where(use, QWAVE[sel][None, :], np.nan)
    with np.errstate(all='ignore'):
        lmin = np.nanmin(lam_used, axis=1); lmax = np.nanmax(lam_used, axis=1)
    nw = winok.sum(axis=1)
    ok = (nw >= 2) & np.isfinite(lmin) & (lmax / lmin >= CONT_MIN_RATIO) & (a0 > 0) & (jb > 0) & (jb < len(ALPHA_GRID) - 1)
    with np.errstate(all='ignore'):
        out = {'ALPHA_LAM': np.where(ok, alpha, np.nan), 'ALPHA_LAM_ERR': np.where(ok, alpha_err, np.nan),
               'CONT_LOGF3000': np.where(ok, np.log10(np.where(a0 > 0, a0, 1)), np.nan), 'CONT_NWIN': nw.astype(np.int16),
               'CONT_LAMMIN': np.where(ok, lmin, np.nan), 'CONT_LAMMAX': np.where(ok, lmax, np.nan),
               'CONT_CHI2NU': np.where(ok & (ndof > 0), cmin / np.where(ndof > 0, ndof, 1), np.nan)}
        cont = np.where(ok[:, None], a0[:, None] * (QWAVE[None, :] / 3000.0)**alpha[:, None], np.nan)   # [n, QNREST]
    return out, cont


# ---------------------------------------------------------------- broad lines
def fit_qlines(rflux, rivar):
    n = rflux.shape[0]; out = {}
    for g, spec in QLINES.items():
        lines, sigmas = spec[0], spec[1]; sideband = spec[2] if len(spec) > 2 else SIDEBAND_FRAC
        lam = np.array([l[1] for l in lines]); lo, hi = lam.min(), lam.max()
        core = 3 * sigmas[-1] / CKMS
        sel = (QWAVE >= lo * (1 - core - sideband)) & (QWAVE <= hi * (1 + core + sideband))
        w = QWAVE[sel]; y = rflux[:, sel].astype(float); wt = rivar[:, sel].astype(float)
        wmid = np.sqrt(w[0] * w[-1]); x = (w - wmid) / (0.1 * wmid)
        p = 2 + len(lam)
        best_chi2 = np.full(n, np.inf); best_coef = np.zeros((n, p)); best_cov = np.zeros((n, p, p)); best_sig = np.full(n, np.nan)
        for sig in sigmas:
            cols = [np.ones_like(w), x] + [np.exp(-0.5 * ((w - l) / (l * sig / CKMS))**2) for l in lam]
            A = np.stack(cols, axis=1)
            AtW = A.T[None, :, :] * wt[:, None, :]
            N = AtW @ A + 1e-9 * np.eye(p)[None]; bvec = np.einsum('npk,nk->np', AtW, y)
            coef = np.linalg.solve(N, bvec[..., None])[..., 0]
            cov = np.linalg.inv(N)
            chi2 = np.sum(wt * (y - coef @ A.T)**2, axis=1)
            better = chi2 < best_chi2
            best_coef[better] = coef[better]; best_cov[better] = cov[better]; best_chi2[better] = chi2[better]; best_sig[better] = sig
        # enough good pixels in the core of the line and in both sidebands
        core_px = (w >= lo * (1 - core)) & (w <= hi * (1 + core))
        ok = ((wt[:, core_px] > 0).sum(axis=1) > 10) & ((wt[:, w < lo * (1 - core)] > 0).sum(axis=1) > 5) \
             & ((wt[:, w > hi * (1 + core)] > 0).sum(axis=1) > 5)
        sig = np.where(np.isnan(best_sig), sigmas[0], best_sig)
        cont_all = best_coef[:, 0:1] + best_coef[:, 1:2] * x[None, :]
        resid = y - cont_all
        for j, (name, l) in enumerate(lines):
            amp = best_coef[:, 2 + j]; amp_err = np.sqrt(np.abs(best_cov[:, 2 + j, 2 + j]))
            width = l * sig / CKMS
            flux = amp * width * np.sqrt(2 * np.pi); flux_err = amp_err * width * np.sqrt(2 * np.pi)
            cont = best_coef[:, 0] + best_coef[:, 1] * (l - wmid) / (0.1 * wmid)
            with np.errstate(all='ignore'):
                ew = np.where(cont > 0, flux / cont, np.nan); ew_err = np.where(cont > 0, flux_err / np.abs(cont), np.nan)
                # centroid of the continuum-subtracted data within +-2.5 sigma (only for the dominant line of a group)
                win = np.abs(w[None, :] - l) <= 2.5 * width[:, None]
                if len(lines) > 1:   # keep the window away from the other line(s) of the group
                    for l2 in lam:
                        if l2 != l:
                            win &= np.abs(w[None, :] - l) < np.abs(w[None, :] - l2)
                m = win & (wt > 0) & (resid > 0)
                num = np.sum(np.where(m, resid * w[None, :], 0), axis=1); den = np.sum(np.where(m, resid, 0), axis=1)
                lc = np.where(den > 0, num / den, np.nan)
                voff = (lc / l - 1) * CKMS
                snr = flux / flux_err
                # peak of the line: continuum-subtracted data smoothed with a Gaussian of sigma_best/2, argmax within
                # +-2 sigma, parabolic refinement over 3 pixels (rest-frame log grid => constant velocity step)
                lpk = np.full(n, np.nan)
                for sg in sigmas:
                    rows = np.where((sig == sg) & ok)[0]
                    if len(rows) == 0:
                        continue
                    spx = max(1.0, (sg / 2) / (DLOG * np.log(10) * CKMS))          # kernel sigma in pixels
                    kx = np.arange(-int(3 * spx), int(3 * spx) + 1); ker = np.exp(-0.5 * (kx / spx)**2); ker /= ker.sum()
                    r = np.where(wt[rows] > 0, resid[rows], 0.0)
                    sm = np.apply_along_axis(lambda v: np.convolve(v, ker, mode='same'), 1, r)
                    win2 = np.abs(w[None, :] - l) <= 2 * l * sg / CKMS
                    if len(lines) > 1:
                        for l2 in lam:
                            if l2 != l:
                                win2 &= np.abs(w[None, :] - l) < np.abs(w[None, :] - l2)
                    sm = np.where(win2, sm, -np.inf)
                    j = np.argmax(sm, axis=1); jm = np.clip(j - 1, 0, len(w) - 1); jp = np.clip(j + 1, 0, len(w) - 1)
                    ii = np.arange(len(rows))
                    y0, ym, yp = sm[ii, j], sm[ii, jm], sm[ii, jp]
                    den2 = ym - 2 * y0 + yp
                    dx = np.where((den2 < 0) & np.isfinite(ym) & np.isfinite(yp), 0.5 * (ym - yp) / den2, 0.0)
                    lpk[rows] = np.where(np.isfinite(y0) & (y0 > 0), w[j] * 10**(dx * DLOG), np.nan)
                vpeak = (lpk / l - 1) * CKMS
            good_line = ok & np.isfinite(flux)
            out[f'{name}_VPEAK'] = np.where(good_line & (snr > 3), vpeak, np.nan)
            out[f'{name}_FLUX'] = np.where(good_line, flux, np.nan); out[f'{name}_FLUX_ERR'] = np.where(good_line, flux_err, np.nan)
            out[f'{name}_EW'] = np.where(good_line, ew, np.nan); out[f'{name}_EW_ERR'] = np.where(good_line, ew_err, np.nan)
            out[f'{name}_CONT'] = np.where(good_line, cont, np.nan); out[f'{name}_SNR'] = np.where(good_line, snr, np.nan)
            out[f'{name}_VOFF'] = np.where(good_line & (snr > 3), voff, np.nan)
            out[f'{name}_SIGMA'] = np.where(good_line, best_sig, np.nan)
        out[f'{g}_CHI2'] = np.where(ok, best_chi2, np.nan)
    return out


# ---------------------------------------------------------------- forest and BAL proxies
def forest_and_bal(rflux, rivar, cont):
    n = rflux.shape[0]; out = {}
    good = rivar > 0
    s = (QWAVE >= FOREST[0]) & (QWAVE < FOREST[1])
    g = good[:, s]; nf = g.sum(axis=1)
    with np.errstate(all='ignore'):
        fsum = np.sum(np.where(g, rflux[:, s], 0), axis=1); csum = np.sum(np.where(g, cont[:, s], 0), axis=1)
        out['FOREST_NPIX'] = nf.astype(np.int32)
        out['FOREST_TRANS'] = np.where((nf >= 20) & (csum > 0), fsum / csum, np.nan)
    lam_good = np.where(good, QWAVE[None, :], np.nan)
    with np.errstate(all='ignore'):
        out['REST_LAMMIN'] = np.nanmin(lam_good, axis=1); out['REST_LAMMAX'] = np.nanmax(lam_good, axis=1)
    # troughs: AI (Hall+02: >=450 km/s, -25000<v<0) and BI (Weymann+91: >=2000 km/s, -25000<v<-3000) on the
    # smoothed normalized flux, requiring the smoothed flux to be >1 sigma below 0.9 C. Lines: C IV, Si IV (HiBAL),
    # Mg II, Al III (LoBAL). BAL_CLASS: 0 non-BAL, 1 weak HiBAL (AI>0, BI=0), 2 strong HiBAL (BI_CIV>0 or BI_SIIV>0),
    # 3 LoBAL (AI_MGII>0 or AI_ALIII>0, any C IV); -1 = C IV region not covered (z<1.5) and no LoBAL test possible.
    k = np.ones(BAL_SMOOTH) / BAL_SMOOTH
    pix_kms = DLOG * np.log(10) * CKMS
    def troughs(fn, es, gd, vs, vlo, vhi, minwidth):
        """fn, es smoothed normalized flux and error [n, m]; returns index (km/s), ntrough, vmin, vmax of deepest."""
        below = gd & np.isfinite(fn) & (fn + BAL_NSIG * es < 1) & (vs >= vlo)[None, :] & (vs <= vhi)[None, :]
        minpx = int(np.ceil(minwidth / pix_kms)); dv = np.gradient(vs)
        idx = np.zeros(n); ntr = np.zeros(n, np.int16); vmin = np.full(n, np.nan); vmax = np.full(n, np.nan)
        for i in np.where(below.sum(axis=1) >= minpx)[0]:
            d = np.diff(np.concatenate([[0], below[i].astype(np.int8), [0]]))
            starts = np.where(d == 1)[0]; stops = np.where(d == -1)[0]; deep = 0.0
            for st, sp_ in zip(starts, stops):
                if sp_ - st < minpx:
                    continue
                seg = slice(st, sp_); a_seg = np.sum((1 - fn[i, seg]) * dv[seg])
                idx[i] += a_seg; ntr[i] += 1
                if a_seg > deep:
                    deep = a_seg; vmin[i] = vs[st]; vmax[i] = vs[sp_ - 1]
        return idx, ntr, vmin, vmax
    cov = {}
    for tag, l0 in (('CIV', 1549.06), ('SIIV', 1399.8), ('MGII', 2798.75), ('ALIII', 1857.4)):
        v = (QWAVE / l0 - 1) * CKMS
        s = (v >= BAL_VRANGE[0]) & (v <= BAL_VRANGE[1]); vs = v[s]
        with np.errstate(all='ignore'):
            fn = rflux[:, s] / (BAL_DEPTH * cont[:, s])
            fn_err = 1.0 / np.sqrt(np.where(rivar[:, s] > 0, rivar[:, s], np.nan)) / (BAL_DEPTH * cont[:, s])
            fs = np.apply_along_axis(lambda r: np.convolve(np.nan_to_num(r, nan=1.0), k, mode='same'), 1, fn)
            es = np.sqrt(np.apply_along_axis(lambda r: np.convolve(np.nan_to_num(r, nan=1e9)**2, k**2, mode='same'), 1, fn_err))
            fn = np.where(np.isfinite(fn), fs, np.nan)
        gd = good[:, s]
        covered = (gd.sum(axis=1) >= 0.8 * s.sum()) & np.isfinite(cont[:, s]).all(axis=1)
        cov[tag] = covered
        ai, nai, vmin, vmax = troughs(fn, es, gd, vs, BAL_VRANGE[0], 0.0, BAL_MIN_WIDTH)
        out[f'BAL_AI_{tag}'] = np.where(covered, ai, np.nan); out[f'BAL_NTROUGH_{tag}'] = np.where(covered, nai, -1).astype(np.int16)
        out[f'BAL_VMIN_{tag}'] = np.where(covered, vmin, np.nan); out[f'BAL_VMAX_{tag}'] = np.where(covered, vmax, np.nan)
        if tag in ('CIV', 'SIIV'):
            bi, nbi, _, _ = troughs(fn, es, gd, vs, BAL_VRANGE[0], -3000.0, 2000.0)
            out[f'BAL_BI_{tag}'] = np.where(covered, bi, np.nan); out[f'BAL_NBI_{tag}'] = np.where(covered, nbi, -1).astype(np.int16)
    hi_cov = cov['CIV'] | cov['SIIV']; lo_cov = cov['MGII'] | cov['ALIII']
    ai_hi = np.nan_to_num(out['BAL_AI_CIV']) + np.nan_to_num(out['BAL_AI_SIIV']); bi_hi = np.nan_to_num(out['BAL_BI_CIV']) + np.nan_to_num(out['BAL_BI_SIIV'])
    ai_lo = np.nan_to_num(out['BAL_AI_MGII']) + np.nan_to_num(out['BAL_AI_ALIII'])
    cls = np.where(hi_cov | lo_cov, 0, -1).astype(np.int8)
    cls = np.where(hi_cov & (ai_hi > 0), 1, cls); cls = np.where(hi_cov & (bi_hi > 0), 2, cls)
    # LoBAL: a Mg II / Al III trough together with a C IV/Si IV trough when the C IV region is covered (LoBALs are a
    # subset of HiBALs); without C IV coverage (z<1.5) only a strong Mg II trough (AI > 1000 km/s) counts; Al III alone is not
    # used there (at z~1 it sits at the blue edge of the spectrum where the continuum is an extrapolation)
    ai_mg = np.nan_to_num(out['BAL_AI_MGII'])
    cls = np.where(lo_cov & (ai_lo > 0) & ((hi_cov & (ai_hi > 0)) | (~hi_cov & (ai_mg > 1000))), 3, cls)
    out['BAL_CLASS'] = cls
    return out




# ---------------------------------------------------------------- SDSS: spectra already on the log10 grid
def rest_from_loggrid(rid0, flux, ivar, mask, z, scale, ebv):
    """spPlate arrays [n, npix] on the log10 grid starting at pixel id rid0 (= round(COEFF0/DLOG); COEFF1 must be
    1e-4). scale [n]: multiplicative calibration (CAL_A); ebv [n]: SFD E(B-V). Flux unit 1e-17 -> cgs.
    Returns rflux, rivar [n, QNREST] = F_lam * (1+z), dereddened (F99, R_V=3.1, 0.86 E(B-V)), shifted by Z."""
    n, npix = flux.shape
    w = 10.0**((rid0 + np.arange(npix)) * DLOG)
    good = (ivar > 0) & (mask == 0)
    corr = (np.nan_to_num(scale, nan=1.0) * passB.FLUX_UNIT)[:, None]
    A1 = passB.extinction.fitzpatrick99(w, 1.0, 3.1)
    A = (3.1 * 0.86 * np.nan_to_num(ebv))[:, None] * A1[None, :]
    f = flux * corr * 10**(0.4 * A); iv = np.where(good, ivar / corr**2 * 10**(-0.8 * A), 0.0)
    rflux = np.zeros((n, QNREST), F32); rivar = np.zeros((n, QNREST), F32)
    zok = np.isfinite(z) & (z > -0.005)
    shift = np.where(zok, np.round(np.log10(1 + np.clip(z, -0.004, None)) / DLOG), 0).astype(int)
    rid_obs = rid0 + np.arange(npix)
    for i in np.where(zok)[0]:
        r0 = rid_obs - shift[i] - QRID_MIN
        ok = (r0 >= 0) & (r0 < QNREST)
        rflux[i, r0[ok]] = f[i, ok] * (1 + z[i]); rivar[i, r0[ok]] = iv[i, ok] / (1 + z[i])**2
    return rflux, rivar


# ---------------------------------------------------------------- all measurements on the rest spectrum
_PB_SL = slice(passB.RID_MIN - QRID_MIN, passB.RID_MAX - QRID_MIN + 1)   # pass-B rest grid (1400-9900) inside ours

def measure_rest(df, rflux, rivar):
    """Adds the quasar columns (and the pass-B galaxy lines + Dn4000 when not already present) to df."""
    if 'HBETA_FLUX' not in df.columns:
        for k, v in passB.fit_lines(passB.RWAVE, rflux[:, _PB_SL], rivar[:, _PB_SL]).items():
            df[k] = v
        d, de = passB.dn4000(passB.RWAVE, rflux[:, _PB_SL], rivar[:, _PB_SL])
        df['DN4000'], df['DN4000_ERR'] = d, de
    pl, cont = fit_powerlaw(rflux, rivar)
    for k, v in pl.items():
        df[k] = v
    for k, v in fit_qlines(rflux, rivar).items():
        df[k] = v
    for k, v in forest_and_bal(rflux, rivar, cont).items():
        df[k] = v
    z = df['Z'].to_numpy(dtype=float)
    for zc, line in LINE_Z.items():
        df[zc] = z + (1 + z) * df[f'{line}_VOFF'].to_numpy(dtype=float) / CKMS          # centroid-based
        df[zc + '_PK'] = z + (1 + z) * df[f'{line}_VPEAK'].to_numpy(dtype=float) / CKMS  # peak-based
    df['Z_OI_PK'] = z + (1 + z) * df['OI_1304_VPEAK'].to_numpy(dtype=float) / CKMS
    # systemic redshift (Nao, 2026-10-24): peak of [O III] 5007 (z<0.9), else peak of Mg II (0.9-2.3), else peak of
    # O I 1304 when detected (z>2.3), else the pipeline z. Z_SYS_SRC: 1 OIII, 2 MgII, 3 OI, 0 pipeline.
    zsys = z.copy(); src = np.zeros(len(df), np.int8)
    for code, col in ((3, 'Z_OI_PK'), (2, 'Z_MGII_PK'), (1, 'Z_OIII_PK')):
        v = df[col].to_numpy(dtype=float); use = np.isfinite(v) & (np.abs(v - z) < 0.01 * (1 + z))
        zsys = np.where(use, v, zsys); src = np.where(use, code, src)
    df['Z_SYS'], df['Z_SYS_SRC'] = zsys, src
    return df


def measure_qso_chunk_desi(df, w, flux, ivar, mask):
    """DESI: pass-B chain (steps 1-5) then the quasar measurements on the quasar rest grid."""
    df = passB.measure_chunk(df, w, flux, ivar, mask)
    rflux, rivar = quasar_rest(df, w, flux, ivar, mask)
    return measure_rest(df, rflux, rivar)
