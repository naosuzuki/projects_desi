"""Shared rest-frame resampling for the galaxy, quasar and star threads (v2 of the read, Nao 2026-10-29): sub-pixel
de-redshift onto the common rest-frame log10 lattice (dlog = 1e-4, 69.0 km/s per pixel, pixel id rid <-> log10 lambda = rid * DLOG,
the passB / sdsslib astrolib convention), in ONE resampling step from the pipeline spectrum.

Before (v1): DESI was rebinned onto the OBSERVED log lattice, then shifted by round(log10(1+z)/dlog) whole pixels; SDSS was shifted
by whole pixels on its native lattice. Rounding misregisters each spectrum by up to +-0.5 px (+-35 km/s; ~20 km/s rms blur in a stack).

DESI (desi_rest): native 0.8 A linear pixels (combine_cameras grid). Rest bin rid spans 10^((rid -+ 0.5) DLOG) in the rest frame,
i.e. 10^((rid -+ 0.5) DLOG) (1+z) observed; its value is the flux-conserving average of the native pixels by fractional overlap
(edges at geometric midpoints, as passB.rebin), ivar from sum (overlap^2 var). Same 50 % good-coverage rule as passB.rebin.
At z = 0 the result is identical to passB.rebin. No interpolation kernel: one binning step, exact sub-pixel registration.

SDSS (sdss_rest): native log lattice (COEFF1 = dlog), so a rest pixel falls a fixed fraction a of a pixel from a native pixel:
Lanczos-3 interpolation (6 pixels, normalized; variance sum L^2 var), linear where one of the 6 is masked but the two neighbours are
good, the nearer pixel where one of the two neighbours is masked (as programs_stars/step6_stack.rest_shift, which showed linear interpolation blurs a stack more than integer rounding while
Lanczos-3 keeps the line width). a = 0 gives the native pixel unchanged. COEFF0 is used exactly (no rounding to the lattice).

Both return OBSERVED-frame F_lambda units on the rest lattice; the caller applies its own (1+z) convention, calibration and
extinction (extinction belongs to the observed wavelength, so apply it to the native pixels before calling).
"""
import numpy as np

DLOG = 1e-4
F32 = np.float32
# SDSS observed-frame bad pixels (Nao 2026-10-30), per survey, on the observed lattice rid = log10(lambda_vac) / 1e-4.
# DR17 (BOSS, v5_13_2): a sky-subtraction wiggle at 6322.7-6335.8 A, peak 6331.40 A (rid 38015). In 28,394 sky fibres of 300
#   plates the mean chi is +0.24 at 6331.40 with -0.22 / +0.14 lobes, and only 2 % of those pixels are flagged. In galaxies and
#   quasars the median relative residual is +2.4 % at 6331.40, present at every MJD 55000-59000. Absent in DR8.
# DR8: none. A Halpha sky residual at 6563-6566 A is only chi +0.1 in the sky fibres, about 0.14 % of the object
#   continuum, so it is not masked; the DR17 dip there was stellar Halpha absorption of the stars on the plates.
# Masked before the rest-frame shift: good &= sdss_badpix_ok(coeff0, npix, survey).
SDSS_BAD_RID = {'dr17': tuple(range(38009, 38019)), 'dr8': ()}


def sdss_badpix_ok(coeff0, npix, survey, dlog=DLOG):
    """[npix] bool, False on SDSS_BAD_RID[survey] of the native log lattice log10 lambda = coeff0 + i dlog."""
    ok = np.ones(npix, bool); r0 = int(round(coeff0 / dlog))
    for r in SDSS_BAD_RID[survey.lower()]:
        if 0 <= r - r0 < npix:
            ok[r - r0] = False
    return ok


def native_edges(w):
    """pixel edges of a native grid at geometric midpoints (passB.rebin_matrix convention)."""
    e = np.empty(len(w) + 1)
    e[1:-1] = np.sqrt(w[:-1] * w[1:]); e[0] = np.sqrt(w[0]**3 / w[1]); e[-1] = np.sqrt(w[-1]**3 / w[-2])
    return e


try:
    from numba import njit
except ImportError:                                   # numpy fallback below (25x slower)
    njit = None

if njit is not None:
    @njit(cache=True)
    def _sweep(E, f, var, g, lo, hi, rf, ri):
        """one spectrum: native edges E [npix+1]; rest bins with observed edges lo/hi (increasing) -> rf, ri. Two-pointer sweep."""
        npix = len(f); h = 0
        for k in range(len(lo)):
            a = lo[k]; b = hi[k]
            if a < E[0] or b > E[npix]:
                continue
            while h < npix - 1 and E[h + 1] <= a:
                h += 1
            num = 0.0; ws = 0.0; vs = 0.0; W = b - a; j = h
            while j < npix and E[j] < b:
                ov = min(E[j + 1], b) - max(E[j], a)
                if ov > 0:
                    t = ov / W
                    num += t * f[j]; ws += t * g[j]; vs += t * t * var[j]
                j += 1
            if ws > 0.5:
                rf[k] = num / ws
                v = vs / (ws * ws)
                if v > 0:
                    ri[k] = 1.0 / v


def desi_rest(w, flux, ivar, mask, z, rid_min, nrest, chunk=64, dlog=DLOG):
    """w [npix] native wavelengths (vacuum, A); flux/ivar/mask [n, npix]; z [n] -> (rflux, rivar) [n, nrest] float32 on rest pixels
    rid_min .. rid_min + nrest - 1 of the lattice log10 lambda = rid * dlog (rid_min in units of dlog). Rows with non-finite z come
    back empty."""
    flux = np.atleast_2d(flux); ivar = np.atleast_2d(ivar); mask = np.atleast_2d(mask)
    n, npix = flux.shape
    E = native_edges(np.asarray(w, float)); dpix = np.diff(E)
    good = (ivar > 0) & (mask == 0)
    f = np.where(good, flux, 0.0); var = np.where(good, 1.0 / np.where(ivar > 0, ivar, 1.0), 0.0); g = good.astype(float)
    rid = np.arange(rid_min, rid_min + nrest)
    lo0 = 10.0**((rid - 0.5) * dlog); hi0 = 10.0**((rid + 0.5) * dlog)
    z = np.atleast_1d(np.asarray(z, float))
    if njit is not None:
        rflux = np.zeros((n, nrest), F32); rivar = np.zeros((n, nrest), F32)
        rid = np.arange(rid_min, rid_min + nrest); lo0 = 10.0**((rid - 0.5) * dlog); hi0 = 10.0**((rid + 0.5) * dlog)
        for i in range(n):
            if np.isfinite(z[i]):
                _sweep(E, f[i].astype(float), var[i], g[i], lo0 * (1 + z[i]), hi0 * (1 + z[i]), rflux[i], rivar[i])
        return rflux, rivar
    zfin = z[np.isfinite(z)]
    zmax = max(float(zfin.max()), 0.0) if len(zfin) else 0.0
    K = int(np.ceil((hi0[-1] - lo0[-1]) * (1 + zmax) / dpix.min())) + 2     # native pixels a rest bin can touch
    rflux = np.zeros((n, nrest), F32); rivar = np.zeros((n, nrest), F32)
    for s in range(0, n, chunk):
        sl = slice(s, min(s + chunk, n)); zz = z[sl]; m = len(zz)
        zf = np.where(np.isfinite(zz), zz, np.nan)
        lo = lo0[None, :] * (1 + zf[:, None]); hi = hi0[None, :] * (1 + zf[:, None])          # observed-frame bin edges
        inside = (lo >= E[0]) & (hi <= E[-1])
        h0 = np.clip(np.searchsorted(E, np.nan_to_num(lo, nan=E[0]), side='right') - 1, 0, npix - 1)
        num = np.zeros((m, nrest)); wsum = np.zeros((m, nrest)); vsum = np.zeros((m, nrest))
        width = hi - lo
        rows = np.arange(m)[:, None]
        for k in range(K):
            h = h0 + k; hv = np.minimum(h, npix - 1)
            ov = np.minimum(E[hv + 1], hi) - np.maximum(E[hv], lo)
            wt = np.where((h < npix) & (ov > 0) & inside, ov, 0.0) / np.where(inside, width, 1.0)
            gi = g[sl][rows, hv]
            num += wt * f[sl][rows, hv]; wsum += wt * gi; vsum += wt**2 * var[sl][rows, hv]
        ok = wsum > 0.5
        with np.errstate(all='ignore'):
            rflux[sl] = np.where(ok, num / np.where(ok, wsum, 1), 0.0)
            rv = vsum / np.where(ok, wsum, 1)**2
            rivar[sl] = np.where(ok & (rv > 0), 1.0 / np.where(rv > 0, rv, 1), 0.0)
    return rflux, rivar


def lanczos_shift(f, var, good, x0, nrest, step=1.0):
    """one spectrum on a unit-spaced lattice; rest pixel k sits at native position x0 + k * step (x0 real; step = rest pixel /
    native pixel, 1 for the same dlog) -> (out, ivar) [nrest]. var = inf at bad pixels."""
    x = x0 + np.arange(nrest) * step
    j0 = np.floor(x + 1e-9).astype(int); a = np.clip(x - j0, 0.0, None); npix = len(f)
    out = np.zeros(nrest); iv = np.zeros(nrest)
    on = a < 1e-6                                                        # rest pixel on a native pixel: copied
    ok = on & (j0 >= 0) & (j0 < npix); ok[ok] = good[j0[ok]]
    out[ok] = f[j0[ok]]; iv[ok] = 1.0 / var[j0[ok]]
    k6 = np.flatnonzero(~on & (j0 - 2 >= 0) & (j0 + 3 < npix))
    J = j0[k6][:, None] + np.arange(-2, 4)[None, :]
    g6 = good[J].all(axis=1); k = k6[g6]; J = J[g6]
    t = np.arange(-2, 4)[None, :] - a[k][:, None]; L = np.sinc(t) * np.sinc(t / 3); L /= L.sum(axis=1, keepdims=True)
    out[k] = (f[J] * L).sum(axis=1); iv[k] = 1.0 / (var[J] * L**2).sum(axis=1)
    k2 = np.flatnonzero(~on & (j0 >= 0) & (j0 + 1 < npix) & (iv == 0)); jj = j0[k2]
    g2 = good[jj] & good[jj + 1]; k2 = k2[g2]; jj = jj[g2]; b = a[k2]
    out[k2] = (1 - b) * f[jj] + b * f[jj + 1]; iv[k2] = 1.0 / ((1 - b)**2 * var[jj] + b**2 * var[jj + 1])
    # one of the two neighbours masked: the nearer native pixel if it is good (what the integer shift did everywhere). Without
    # this every isolated masked pixel (BOSS sky lines) costs two rest pixels, and band coverage fell below the 70 % of
    # ab_mags for ~2 % of DR17 galaxies (MABS NaN -> no cell).
    k1 = np.flatnonzero(~on & (iv == 0)); jn = j0[k1] + (a[k1] >= 0.5)
    inr = (jn >= 0) & (jn < npix); k1 = k1[inr]; jn = jn[inr]; g1 = good[jn]; k1 = k1[g1]; jn = jn[g1]
    out[k1] = f[jn]; iv[k1] = 1.0 / var[jn]
    return out, iv


def sdss_rest(coeff0, flux, ivar, mask, z, rid_min, nrest, dlog=DLOG, native_dlog=DLOG):
    """spPlate / spec-lite arrays [n, npix] on the native log lattice log10 lambda = coeff0 + i * native_dlog (COEFF1);
    z [n] -> (rflux, rivar) [n, nrest] float32 on rest pixels rid_min .. rid_min + nrest - 1 of the lattice log10 lambda = rid * dlog.
    mask: nonzero = bad (pass the and_mask selection the thread uses). Rows with non-finite z come back empty."""
    flux = np.atleast_2d(flux); ivar = np.atleast_2d(ivar); mask = np.atleast_2d(mask)
    n, npix = flux.shape
    good = (ivar > 0) & (mask == 0)
    var = np.where(good, 1.0 / np.where(ivar > 0, ivar, 1.0), np.inf)
    rflux = np.zeros((n, nrest), F32); rivar = np.zeros((n, nrest), F32)
    z = np.atleast_1d(np.asarray(z, float))
    for i in range(n):
        if not np.isfinite(z[i]):
            continue
        x0 = (rid_min * dlog - coeff0 + np.log10(1 + z[i])) / native_dlog      # native position of rest pixel 0
        o, v = lanczos_shift(flux[i].astype(float), var[i], good[i], x0, nrest, dlog / native_dlog)
        rflux[i] = o; rivar[i] = v
    return rflux, rivar
