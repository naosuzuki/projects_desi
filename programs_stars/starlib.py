"""starlib: per-star measurements for pass S (stellar templates), shared by the DESI driver (02_passS_measure.py) and,
later, the SDSS spPlate driver. Imports the common pieces of pass B (filters, rebin, extinction constants).

Inputs are the pipeline spectra of a coadd (flux in 1e-17 erg/s/cm2/A, calibrated to the total flux of a point source,
which for a star IS the star's flux) on the common 0.8 A grid of passB.combine_cameras, plus the fibermap/catalog rows.

Measurements
  1. Gaia check: synthetic Gaia DR3 G, BP, RP (Vega) from the pipeline spectrum vs GAIA_PHOT_*_MEAN_MAG.
     G and RP reach beyond 9824 A and BP/G below 3600 A, so the spectrum is padded flat (mean of the outer 250 good
     pixels); the fraction of the band weight that lies in the padding is stored (GAIA_EXTRAP_G/BP/RP). Catalog minus
     synthetic: GAIA_DG, GAIA_DBP, GAIA_DRP; color residual GAIA_DBPRP = (BP-RP)_cat - (BP-RP)_syn (tilt handle).
  2. RV frame: integer-pixel shift by the Redrock z on the log10 grid (dlog=1e-4 = 69 km/s), RV_KMS = c z,
     RV_RESID_KMS = residual after the integer shift (stored for a sub-pixel stacking pass if Nao chooses it).
  3. Spectral indices (A, pseudo-continuum EW: positive = absorption; band edges converted from air to vacuum) on the RV-corrected spectrum: Ca II K, HdeltaA,
     G band 4300, HgammaA, Hbeta, Mg b, Fe5270, Fe5335, He II 4686, He I 5876, Na D, Halpha (negative = emission),
     Ca II triplet (Cenarro+01, 3 lines + sum); band ratios TiO5, CaH2, CaH3 (Reid+95; M-dwarf type and gravity).
     Errors from the pixel ivar inside the feature band. NaN when < 80 % of the pixels of any band are good.
  4. Continuum colors on the RV-corrected spectrum (mag): C_BLUE = -2.5 log F(4100)/F(5600), C_RED = -2.5 log F(5600)/F(8100),
     each raw (_RAW) and corrected for Milky Way extinction with F99, 0.86 x E(B-V)_SFD (_DERED).
"""
import os, importlib.util
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('passB', os.path.join(HERE, '..', 'programs_templates', '02_passB_measure.py'))
passB = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(passB)

DLOG = passB.DLOG; CKMS = passB.CKMS; F32 = np.float32; FLUX_UNIT = passB.FLUX_UNIT
# RV-corrected grid: the observed 3600-9824 A range shifted by up to +-1000 km/s
SREST_MIN, SREST_MAX = 3550.0, 9880.0
SRID_MIN = int(round(np.log10(SREST_MIN) / DLOG)); SRID_MAX = int(round(np.log10(SREST_MAX) / DLOG))
SNREST = SRID_MAX - SRID_MIN + 1
SWAVE = 10.0**(np.arange(SRID_MIN, SRID_MAX + 1) * DLOG)

# Gaia DR3 zero points (Riello+21): m_AB - m_Vega = ZP_AB - ZP_Vega
GAIA_AB_MINUS_VEGA = {'G': 25.8010 - 25.6874, 'BP': 25.3540 - 25.3385, 'RP': 25.1040 - 24.7479}
GAIA_FILTERS = {'G': 'gaiadr3-G', 'BP': 'gaiadr3-BP', 'RP': 'gaiadr3-RP'}
PAD_MIN, PAD_MAX = 3000.0, 11000.0        # the Gaia bands span 3300-10500 A; pad flat from the spectrum's ends to here
PAD_NPIX = 250

# index definitions: name -> (feature band, list of continuum bands)  [A]
INDICES = {
    'CAK':      ((3925.0, 3945.0),     [(3903.0, 3923.0), (4000.0, 4020.0)]),
    'HDELTA_A': ((4083.50, 4122.25),   [(4041.60, 4079.75), (4128.50, 4161.00)]),      # Worthey & Ottaviani 97
    'G4300':    ((4281.375, 4316.375), [(4266.375, 4282.625), (4318.875, 4335.125)]),  # Lick
    'HGAMMA_A': ((4319.75, 4363.50),   [(4283.50, 4319.75), (4367.25, 4419.75)]),
    'HEII_4686':((4676.0, 4696.0),     [(4650.0, 4670.0), (4705.0, 4725.0)]),
    'HBETA_LK': ((4847.875, 4876.625), [(4827.875, 4847.875), (4876.625, 4891.625)]),  # Lick Hbeta
    'MGB':      ((5160.125, 5192.625), [(5142.625, 5161.375), (5191.375, 5206.375)]),
    'FE5270':   ((5245.650, 5285.650), [(5233.150, 5248.150), (5285.650, 5318.150)]),
    'FE5335':   ((5312.125, 5352.125), [(5304.625, 5315.875), (5353.375, 5363.375)]),
    'HEI_5876': ((5866.0, 5886.0),     [(5840.0, 5860.0), (5920.0, 5940.0)]),
    'NAD':      ((5876.875, 5909.375), [(5860.625, 5875.625), (5922.125, 5948.125)]),
    'HALPHA_ST':((6554.0, 6574.0),     [(6500.0, 6540.0), (6590.0, 6630.0)]),
    'CAT1':     ((8484.0, 8513.0),     [(8474.0, 8484.0), (8563.0, 8577.0), (8619.0, 8642.0), (8700.0, 8725.0), (8776.0, 8792.0)]),  # Cenarro+01
    'CAT2':     ((8522.0, 8562.0),     [(8474.0, 8484.0), (8563.0, 8577.0), (8619.0, 8642.0), (8700.0, 8725.0), (8776.0, 8792.0)]),
    'CAT3':     ((8642.0, 8682.0),     [(8474.0, 8484.0), (8563.0, 8577.0), (8619.0, 8642.0), (8700.0, 8725.0), (8776.0, 8792.0)]),
}
# band ratios (Reid+95): name -> (numerator band, denominator band); F_num / F_den
RATIOS = {'TIO5': ((7126.0, 7135.0), (7042.0, 7046.0)), 'CAH2': ((6814.0, 6846.0), (7042.0, 7046.0)), 'CAH3': ((6960.0, 6990.0), (7042.0, 7046.0))}
COLOR_BANDS = {'4100': (4000.0, 4200.0), '5600': (5500.0, 5700.0), '8100': (8000.0, 8200.0)}
MIN_GOOD_FRAC = 0.8


def air_to_vac(w):
    """air -> vacuum wavelength (A), IAU standard (Morton 1991 / Ciddor 1996 form used by SDSS)."""
    w = np.asarray(w, float); s = 1e4 / w
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.00015997409212 / (38.92568793293 - s**2)
    return w * n


# the index definitions above are in AIR wavelengths (Lick, Worthey & Ottaviani 97, Cenarro+01, Reid+95); the spectra are in
# VACUUM (SDSS and DESI), so the band edges are converted once here (+1.1 A at 4000 A, +2.4 A at 8600 A) - 2026-10-25
INDICES = {k: (tuple(air_to_vac(v[0])), [tuple(air_to_vac(b)) for b in v[1]]) for k, v in INDICES.items()}
RATIOS = {k: (tuple(air_to_vac(v[0])), tuple(air_to_vac(v[1]))) for k, v in RATIOS.items()}


# ---------------------------------------------------------------- 1. Gaia synthetic photometry
def gaia_mags(w, flux, ivar, mask):
    """Synthetic Gaia DR3 G/BP/RP (Vega) of the pipeline spectrum, padded flat beyond the DESI range.
    Returns dict of arrays: SYN_GAIA_{G,BP,RP}, GAIA_EXTRAP_{G,BP,RP} (band-weight fraction in the padding)."""
    n = flux.shape[0]
    good = (ivar > 0) & (mask == 0)
    f = flux * FLUX_UNIT
    with np.errstate(all='ignore'):
        lo = np.nanmean(np.where(good[:, :PAD_NPIX], f[:, :PAD_NPIX], np.nan), axis=1)
        hi = np.nanmean(np.where(good[:, -PAD_NPIX:], f[:, -PAD_NPIX:], np.nan), axis=1)
    dw = float(np.median(np.diff(w)))
    pad_lo = np.arange(PAD_MIN, w[0] - 0.5 * dw, dw); pad_hi = np.arange(w[-1] + dw, PAD_MAX, dw)
    wz = np.concatenate([pad_lo, w, pad_hi])
    fz = np.concatenate([np.repeat(np.nan_to_num(lo)[:, None], len(pad_lo), axis=1), np.where(good, f, 0.0),
                         np.repeat(np.nan_to_num(hi)[:, None], len(pad_hi), axis=1)], axis=1)
    gz = np.concatenate([np.ones((n, len(pad_lo)), bool), good, np.ones((n, len(pad_hi)), bool)], axis=1)
    inspec = np.zeros(len(wz), bool); inspec[len(pad_lo):len(pad_lo) + len(w)] = True
    out = {}
    for b, name in GAIA_FILTERS.items():
        flt = passB.filt(name)
        wt, den, inband = passB.band_weights(flt, wz)
        m, _ = passB.ab_mags(flt, wz, fz, gz, min_range=0.0, min_good=0.5)
        out[f'SYN_GAIA_{b}'] = m - GAIA_AB_MINUS_VEGA[b]
        out[f'GAIA_EXTRAP_{b}'] = np.full(n, float(wt[~inspec].sum() / wt.sum()))
        ok = np.isfinite(lo) & np.isfinite(hi)
        out[f'SYN_GAIA_{b}'] = np.where(ok, out[f'SYN_GAIA_{b}'], np.nan)
    return out


# ---------------------------------------------------------------- 1b. generic synthetic AB magnitudes (any speclite filter)
def synth_mags(w, flux, ivar, mask, names, prefix):
    """AB magnitudes of the raw pipeline spectrum through the speclite filters `names` (dict band -> filter name), with the
    spectrum padded flat beyond its ends (as gaia_mags) so bands that overhang the spectrum still get a value; the band-weight
    fraction in the padding is stored as <prefix>_EXTRAP_<band> (use it to down-weight or drop the band in the joint fit)."""
    n = flux.shape[0]
    good = (ivar > 0) & (mask == 0)
    f = flux * FLUX_UNIT; iv = ivar / FLUX_UNIT**2
    with np.errstate(all='ignore'):
        lo = np.nanmean(np.where(good[:, :PAD_NPIX], f[:, :PAD_NPIX], np.nan), axis=1)
        hi = np.nanmean(np.where(good[:, -PAD_NPIX:], f[:, -PAD_NPIX:], np.nan), axis=1)
    dw = float(np.median(np.diff(w)))
    pad_lo = np.arange(PAD_MIN, w[0] - 0.5 * dw, dw); pad_hi = np.arange(w[-1] + dw, PAD_MAX, dw)
    wz = np.concatenate([pad_lo, w, pad_hi])
    fz = np.concatenate([np.repeat(np.nan_to_num(lo)[:, None], len(pad_lo), axis=1), np.where(good, f, 0.0),
                         np.repeat(np.nan_to_num(hi)[:, None], len(pad_hi), axis=1)], axis=1)
    gz = np.concatenate([np.ones((n, len(pad_lo)), bool), good, np.ones((n, len(pad_hi)), bool)], axis=1)
    ivz = np.concatenate([np.zeros((n, len(pad_lo))), np.where(good, iv, 0.0), np.zeros((n, len(pad_hi)))], axis=1)
    inspec = np.zeros(len(wz), bool); inspec[len(pad_lo):len(pad_lo) + len(w)] = True
    out = {}
    for b, name in names.items():
        flt = passB.filt(name); wt, den, inband = passB.band_weights(flt, wz)
        m, e = passB.ab_mags(flt, wz, fz, gz, min_range=0.0, min_good=0.5, ivar2d=ivz)
        ok = np.isfinite(lo) & np.isfinite(hi)
        out[f'{prefix}_{b}'] = np.where(ok, m, np.nan); out[f'{prefix}_{b}_ERR'] = np.where(ok, e, np.nan)
        out[f'{prefix}_EXTRAP_{b}'] = np.full(n, float(wt[~inspec].sum() / wt.sum()))
    return out


SDSS_FILTERS = {'U': 'sdss2010-u', 'G': 'sdss2010-g', 'R': 'sdss2010-r', 'I': 'sdss2010-i', 'Z': 'sdss2010-z'}
LS_FILTERS = {True: {'G': 'decam2014-g', 'R': 'decam2014-r', 'Z': 'decam2014-z'}, False: {'G': 'BASS-g', 'R': 'BASS-r', 'Z': 'MzLS-z'}}


# ---------------------------------------------------------------- 2. RV frame
def rv_frame(df, w, flux, ivar, mask, rid0=None):
    """Rebin the pipeline spectrum (cgs) to the log grid (or, for SDSS spPlates already on the log10 grid with
    dlog = 1e-4, take the pixels as they are: rid0 = round(COEFF0/DLOG)) and shift by the integer pixel closest to the
    pipeline z. Returns rflux, rivar [n, SNREST] on SWAVE (F_lam, no (1+z) factor: stars), plus RV_KMS, RV_RESID_KMS."""
    n = len(df)
    f = flux * FLUX_UNIT; iv = ivar / FLUX_UNIT**2
    if rid0 is None:
        rid_obs, rflux_o, rivar_o = passB.rebin(w, f, iv, mask)
    else:
        good = (ivar > 0) & (mask == 0)
        rid_obs = rid0 + np.arange(flux.shape[1]); rflux_o = np.where(good, f, 0.0); rivar_o = np.where(good, iv, 0.0)
    z = np.nan_to_num(df['Z'].to_numpy(dtype=float), nan=0.0)
    z = np.clip(z, -0.004, 0.004)                                   # |RV| < 1200 km/s: anything else is not a star's RV
    shift = np.round(np.log10(1 + z) / DLOG).astype(int)
    rflux = np.zeros((n, SNREST), F32); rivar = np.zeros((n, SNREST), F32)
    for s in np.unique(shift):
        rows = np.where(shift == s)[0]
        r0 = rid_obs - s - SRID_MIN; ok = (r0 >= 0) & (r0 < SNREST)
        rflux[np.ix_(rows, r0[ok])] = rflux_o[rows][:, ok]; rivar[np.ix_(rows, r0[ok])] = rivar_o[rows][:, ok]
    rv = df['Z'].to_numpy(dtype=float) * CKMS
    resid = (z - (10.0**(shift * DLOG) - 1)) * CKMS
    return rflux, rivar, rv, resid


# ---------------------------------------------------------------- 3. indices and ratios
def _band_mean(rflux, rivar, lo, hi):
    """mean flux and its variance in [lo, hi) over good pixels; NaN if < MIN_GOOD_FRAC of the band is good."""
    s = (SWAVE >= lo) & (SWAVE < hi); g = rivar[:, s] > 0
    n = g.sum(axis=1); nn = np.where(n > 0, n, 1); ok = n >= MIN_GOOD_FRAC * s.sum()
    with np.errstate(all='ignore'):
        m = np.sum(np.where(g, rflux[:, s], 0), axis=1) / nn
        v = np.sum(np.where(g, 1.0 / np.where(g, rivar[:, s], 1), 0), axis=1) / nn**2
    return np.where(ok, m, np.nan), np.where(ok, v, np.nan), ok


def indices(rflux, rivar):
    n = rflux.shape[0]; out = {}
    rflux = rflux.astype(float); rivar = rivar.astype(float)
    for name, (feat, conts) in INDICES.items():
        # pseudo-continuum: straight line through the band means (2 bands: exact; >2: least squares)
        # pseudo-continuum through the good band means (bands with < MIN_GOOD_FRAC good pixels are dropped; at least
        # 2 bands, or 60 % of them for the 5-band Ca II triplet definition, must remain: sky residuals mask many pixels there)
        xs, ys, oks = [], [], []
        for lo, hi in conts:
            m, v, ok = _band_mean(rflux, rivar, lo, hi); xs.append(0.5 * (lo + hi)); ys.append(m); oks.append(ok)
        xs = np.array(xs); ys = np.stack(ys, axis=1); okb = np.stack(oks, axis=1)
        need = max(2, int(np.ceil(0.6 * len(conts))))
        wb = okb.astype(float); nb = wb.sum(axis=1); oks = nb >= need
        with np.errstate(all='ignore'):
            xm = (wb * xs[None, :]).sum(axis=1) / nb; ym = np.nansum(wb * ys, axis=1) / nb
            sxx = (wb * (xs[None, :] - xm[:, None])**2).sum(axis=1)
            slope = np.where(sxx > 0, np.nansum(wb * (xs[None, :] - xm[:, None]) * (ys - ym[:, None]), axis=1) / np.where(sxx > 0, sxx, 1), 0.0)
            icpt = ym - slope * xm
        s = (SWAVE >= feat[0]) & (SWAVE < feat[1]); wf = SWAVE[s]
        g = rivar[:, s] > 0; okf = g.sum(axis=1) >= MIN_GOOD_FRAC * s.sum()
        cont = icpt[:, None] + slope[:, None] * wf[None, :]
        dlam = np.gradient(wf)
        with np.errstate(all='ignore'):
            frac = np.where(g, 1.0 - rflux[:, s] / cont, 0.0)
            # masked pixels inside the feature: assume the mean depth of the good ones (scale by band width / good width)
            ew = np.sum(frac * dlam[None, :], axis=1) * (dlam.sum() / np.where(g.sum(axis=1) > 0, np.sum(np.where(g, dlam[None, :], 0), axis=1), np.nan))
            var = np.sum(np.where(g, dlam[None, :]**2 / (np.where(g, rivar[:, s], 1) * cont**2), 0), axis=1)
        ok = oks & okf & np.isfinite(ew) & (icpt + slope * wf.mean() > 0)
        out[f'{name}_EW'] = np.where(ok, ew, np.nan); out[f'{name}_EW_ERR'] = np.where(ok, np.sqrt(var), np.nan)
    with np.errstate(all='ignore'):
        out['CAT_EW'] = out['CAT1_EW'] + out['CAT2_EW'] + out['CAT3_EW']
        out['CAT_EW_ERR'] = np.sqrt(out['CAT1_EW_ERR']**2 + out['CAT2_EW_ERR']**2 + out['CAT3_EW_ERR']**2)
    for name, (num, den) in RATIOS.items():
        mn, vn, okn = _band_mean(rflux, rivar, *num); md, vd, okd = _band_mean(rflux, rivar, *den)
        with np.errstate(all='ignore'):
            r = mn / md; err = np.abs(r) * np.sqrt(vn / mn**2 + vd / md**2)
        ok = okn & okd & (md > 0)
        out[name] = np.where(ok, r, np.nan); out[f'{name}_ERR'] = np.where(ok, err, np.nan)
    return out


# ---------------------------------------------------------------- 4. colors
def colors(rflux, rivar, ebv):
    out = {}
    means = {k: _band_mean(rflux.astype(float), rivar.astype(float), lo, hi)[0] for k, (lo, hi) in COLOR_BANDS.items()}
    with np.errstate(all='ignore'):
        cb = -2.5 * np.log10(means['4100'] / means['5600']); cr = -2.5 * np.log10(means['5600'] / means['8100'])
    out['C_BLUE_RAW'], out['C_RED_RAW'] = cb, cr
    # F99 (R_V=3.1) A(lam)/A_V at the band centers; A_V = 3.1 * 0.86 * EBV
    a = {k: float(passB.extinction.fitzpatrick99(np.array([0.5 * (lo + hi)]), 1.0, 3.1)[0]) for k, (lo, hi) in COLOR_BANDS.items()}
    av = 3.1 * 0.86 * np.nan_to_num(ebv, nan=0.0)
    out['C_BLUE_DERED'] = cb - av * (a['4100'] - a['5600']); out['C_RED_DERED'] = cr - av * (a['5600'] - a['8100'])
    return out


# ---------------------------------------------------------------- all together
# Gaia extinction per unit E(B-V) at R_V = 3.1 (Casagrande & VandenBerg 2018, solar-type star): A_G, A_BP, A_RP
AG_PER_EBV, ABP_PER_EBV, ARP_PER_EBV = 2.740, 3.374, 2.035
EBPRP_PER_EBV = ABP_PER_EBV - ARP_PER_EBV


def distance_columns(df):
    """Nao's steps 3b and 4: distance from the Gaia parallax, Milky Way extinction, scale to absolute flux (10 pc).
    Parallax: PLX_CORR (SDSS catalog) or PARALLAX + 0.017 mas (DESI fibermap). EBV_3D is NaN until a 3D dust map
    (Bayestar19 / Edenhofer+24 through dustmaps) is installed; EBV_USED then falls back to 0.86 x SFD (EBV_SRC = 1; 2 = 3D map).
    The dereddening itself is a per-star factor applied at stack time, so the tables can be re-dereddened without re-reading."""
    n = len(df); out = {}
    if 'PLX_CORR' in df:
        plx = df['PLX_CORR'].to_numpy(dtype=float)
    elif 'PARALLAX' in df:
        plx = df['PARALLAX'].to_numpy(dtype=float) + 0.017
    else:
        plx = np.full(n, np.nan)
    with np.errstate(all='ignore'):
        plx = np.where(np.isfinite(plx) & (plx > 0), plx, np.nan)
        d = 1000.0 / plx                                           # pc
        out['DIST_PC'] = d; out['DM'] = 5 * np.log10(d / 10.0)
        out['ABS_SCALE'] = (d / 10.0)**2                          # multiply F_lam by this to put the star at 10 pc
        ebv_sfd = np.nan_to_num(df['EBV'].to_numpy(dtype=float), nan=0.0) if 'EBV' in df else np.zeros(n)
        out['EBV_SFD'] = ebv_sfd
        out['EBV_3D'] = df['EBV_3D'].to_numpy(dtype=float) if 'EBV_3D' in df else np.full(n, np.nan)   # step3b_ebv3d.py
        out['EBV_USED'] = np.where(np.isfinite(out['EBV_3D']), out['EBV_3D'], 0.86 * ebv_sfd)
        out['EBV_SRC'] = np.where(np.isfinite(out['EBV_3D']), 2, 1).astype(np.int8)
        # A_G from the same E(B-V) as the spectrum (policy 2026-10-24: the 3D map is the correction; GSP-Phot A_G,
        # carried in GAIA_AG_GSPPHOT where present, is a per-star cross-check only)
        out['A_G_USED'] = AG_PER_EBV * out['EBV_USED']; out['EBPRP_USED'] = EBPRP_PER_EBV * out['EBV_USED']
    return out


def measure_stars(df, w, flux, ivar, mask, rid0=None):
    """df: STAR rows (with Z, EBV, GAIA_PHOT_*_MEAN_MAG, a parallax column); flux/ivar/mask [n, npix] pipeline spectra
    on grid w (rid0: log-grid pixel id of w[0] when the spectra are already on the log grid, e.g. spPlates).
    Adds the pass-S columns to df and returns it."""
    for k, v in gaia_mags(w, flux, ivar, mask).items():
        df[k] = v
    with np.errstate(all='ignore'):
        for b in ('G', 'BP', 'RP'):
            cat = df[f'GAIA_PHOT_{b}_MEAN_MAG'].to_numpy(dtype=float) if f'GAIA_PHOT_{b}_MEAN_MAG' in df else np.full(len(df), np.nan)
            cat = np.where(np.isfinite(cat) & (cat > 0) & (cat < 30), cat, np.nan)
            df[f'GAIA_D{b}'] = cat - df[f'SYN_GAIA_{b}'].to_numpy(dtype=float)
        df['GAIA_DBPRP'] = df['GAIA_DBP'] - df['GAIA_DRP']
    rflux, rivar, rv, resid = rv_frame(df, w, flux, ivar, mask, rid0=rid0)
    df['RV_KMS'], df['RV_RESID_KMS'] = rv, resid
    for k, v in distance_columns(df).items():
        df[k] = v
    with np.errstate(all='ignore'):
        # absolute synthetic G from the spectrum (steps 3b + 4 in one number): compare with Gaia M_G0
        df['M_G_SYN0'] = df['SYN_GAIA_G'] - df['DM'] - df['A_G_USED']
    for k, v in indices(rflux, rivar).items():
        df[k] = v
    for k, v in colors(rflux, rivar, df['EBV'].to_numpy(dtype=float)).items():
        df[k] = v
    good = rivar > 0
    lam = np.where(good, SWAVE[None, :], np.nan)
    with np.errstate(all='ignore'):
        df['RV_LAMMIN'] = np.nanmin(lam, axis=1); df['RV_LAMMAX'] = np.nanmax(lam, axis=1)
    return df
