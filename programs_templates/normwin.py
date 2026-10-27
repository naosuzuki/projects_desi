"""Line-free normalization windows, chosen per redshift (ms.tex Sec. 3.6 item 4 / Sec. 3.9 Step 6).

Each spectrum (and each stack) is normalized to unit mean F_lambda in one rest-frame continuum window that
(1) contains no strong emission or absorption feature, (2) lies fully inside the spectrum's rest-frame coverage
with a 30 A margin, (3) does not fall on the telluric A band (7580-7700 A observed), and (4) among the
remaining candidates is the one whose observed wavelength is closest to 6400 A (high throughput, below the
OH forest). The choice is deterministic in z, so all members of a Delta z = 0.05 slice share one window.

Galaxies - avoided features: [OII] 3727, [NeIII] 3869, Balmer series (Hd 4102 .. Ha 6563), [OIII] 4959/5007,
[NI] 5199, HeI 5876, [OI] 6300, [NII] 6548/6583, [SII] 6716/6731, [ArIII] 7136, [OII] 7320/7330,
[SIII] 9069/9531; Ca H&K 3934/3968, 4000 A break (3850-4100), G band 4304, Mg b 5175, Na D 5893,
Ca triplet 8498-8662, Balmer jump 3646, Mg II 2800 and Fe II 2344/2382/2586/2600 absorption, C III] 1909.
Quasars - standard continuum windows (Vanden Berk et al. 2001; Shen et al. 2011) plus 1700-1705 and 6005-6035.
All wavelengths are rest-frame vacuum Angstrom.
"""
import numpy as np

OBS_LO, OBS_HI = 3600.0, 9824.0          # DESI; SDSS BOSS 3600-10300, DR8 3800-9200: pass obs_range for those
A_BAND = (7580.0, 7700.0)
BEST_OBS = 6400.0
MARGIN = 30.0

GALAXY_WINDOWS = [(5400, 5500), (6000, 6250), (4600, 4800), (3000, 3100), (2100, 2200), (8100, 8300), (7000, 7100), (1950, 2050)]
QSO_WINDOWS = [(1350, 1360), (1445, 1465), (1700, 1705), (2200, 2250), (2650, 2710), (3010, 3040), (4210, 4230),
               (5080, 5130), (6005, 6035), (7020, 7060)]


def pick_window(z, kind='GALAXY', obs_range=(OBS_LO, OBS_HI)):
    """Return (rest_lo, rest_hi) of the normalization window for redshift z, or None if no candidate fits."""
    cands = GALAXY_WINDOWS if kind.upper().startswith('GAL') else QSO_WINDOWS
    lo, hi = obs_range[0] / (1 + z) + MARGIN, obs_range[1] / (1 + z) - MARGIN
    best = None
    for a, b in cands:
        if a < lo or b > hi:
            continue
        oa, ob = a * (1 + z), b * (1 + z)
        if ob > A_BAND[0] and oa < A_BAND[1]:
            continue
        score = abs(0.5 * (oa + ob) - BEST_OBS)
        if best is None or score < best[0]:
            best = (score, (a, b))
    return None if best is None else best[1]


def slice_window(z_lo, z_hi, kind='GALAXY', obs_range=(OBS_LO, OBS_HI)):
    """Window for a redshift slice: chosen at the slice centre, then checked at both edges (None if it fails there)."""
    w = pick_window(0.5 * (z_lo + z_hi), kind, obs_range)
    if w is None:
        return None
    for z in (z_lo, z_hi):
        lo, hi = obs_range[0] / (1 + z) + MARGIN, obs_range[1] / (1 + z) - MARGIN
        if w[0] < lo or w[1] > hi:
            return pick_window(z, kind, obs_range)
    return w


def normalize(rwave, flux, ivar, window, min_good=0.7):
    """Scale factor s so that the mean F_lambda over good pixels in `window` is 1; returns (s, ok).
    flux/ivar may be 1-D or 2-D [n, npix]. ok is False when < min_good of the window pixels are usable."""
    sel = (rwave >= window[0]) & (rwave <= window[1])
    f = np.atleast_2d(flux)[:, sel]; w = np.atleast_2d(ivar)[:, sel] > 0
    n = w.sum(axis=1)
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, np.sum(np.where(w, f, 0.0), axis=1) / np.where(n > 0, n, 1), np.nan)
        s = np.where(mean > 0, 1.0 / mean, np.nan)
    ok = (n >= min_good * sel.sum()) & np.isfinite(s)
    if np.ndim(flux) == 1:
        return float(s[0]), bool(ok[0])
    return s, ok


if __name__ == '__main__':
    print('galaxies'); [print(f'  {z:.2f}-{z+.05:.2f}  {slice_window(z, z+.05)}') for z in np.arange(0, 1.6, 0.05)]
    print('quasars');  [print(f'  {z:.2f}-{z+.25:.2f}  {slice_window(z, z+.25, "QSO")}') for z in np.arange(0, 5.0, 0.25)]
