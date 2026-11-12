#!/usr/bin/env python3
"""
02_speczcalib.py — igmdust/DES: calibrate the photo-z absorption of
line-of-sight reddening in the standard-crayon estimator.

The crayon estimator (01) measures color residuals about the red-sequence
locus evaluated at zredmagic.  zredmagic is itself fit from the same
colors, so reddening E along the line of sight shifts the photo-z by
beta = dz_RM/dE, and the locus subtraction then removes L'_c(z)*beta*E
of the signal.  The measured response per color is

    M_c = dR_c - L'_c(z) * beta        (dust alone would give dR_c),

so the 01 profile is suppressed by  f = sum_c w_c dR_c M_c / sum_c w_c
dR_c^2  and the de-absorbed reddening is  E_true = E_measured / f.

beta is measured empirically from crayons positionally matched (1") to
DESI DR1 LRG spectra: at fixed z_spec, the per-galaxy pseudo-reddening
E_hat (color residual about the z_spec locus, projected on the dust
vector) is uncorrelated with redshift errors UNLESS the photo-z responds
to reddening-like color shifts; the slope of (z_RM - z_spec) vs E_hat is
exactly that response.  L'_c(z) comes from the full-sample locus of 01.

WHERE TO RUN: astro4 or astrofs4/SLURM (~2 min).
    python3 02_speczcalib.py
    python3 02_speczcalib.py --selftest
Output: table + $OUT/speczcalib.npz
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "crayon01", os.path.join(_here, "01_crayon_paircount.py"))
m01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m01)

DESI_BASE = "/data/astrofs2_1/suzuki/data/DESI/catalog/dr1"
OUT_DEFAULT = os.path.expanduser("~/data/igmdust/paircount")
MATCH_ARCSEC = 1.0
DZ_CLIP = 0.15                       # |z_RM - z_spec| outlier cut
EHAT_CLIP = 0.3                      # |E_hat| outlier cut [mag]
DR = np.array([3.214 - 2.165, 2.165 - 1.592, 1.592 - 1.211, 3.214 - 1.592])
CNAMES = [f"{a}-{b}" for a, b in m01.COLORS_DES]


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    c = np.cos(dec)
    return np.column_stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)])


def load_lrg_specz(base):
    ra, dec, z = [], [], []
    for cap in ("NGC", "SGC"):
        with fits.open(f"{base}/LRG_{cap}_clustering.dat.fits") as h:
            d = h[1].data
            ra.append(np.asarray(d["RA"], float))
            dec.append(np.asarray(d["DEC"], float))
            z.append(np.asarray(d["Z"], float))
    ra, dec, z = map(np.concatenate, (ra, dec, z))
    print(f"  DESI DR1 LRG spectra: {len(z):,}")
    return ra, dec, z


def running_median(x, y, nbin=25):
    """Median trend of y vs x, interpolated back onto x."""
    q = np.quantile(x, np.linspace(0, 1, nbin + 1))
    q[-1] += 1e-9
    xm, ym = [], []
    for i in range(nbin):
        m = (x >= q[i]) & (x < q[i + 1])
        if m.sum() >= 50:
            xm.append(np.median(x[m]))
            ym.append(np.median(y[m]))
    return np.interp(x, np.array(xm), np.array(ym))


def run(a):
    print("== crayons (raw colors)")
    cra, cdec, zrm, colors, valid = m01.load_crayons(
        a.des_dir, raw=True, phot="mof" if a.mof else "sof")
    print("== DESI LRG spectroscopic redshifts")
    lra, ldec, lz = load_lrg_specz(a.base)

    chord = 2 * np.sin(np.radians(MATCH_ARCSEC / 3600) / 2)
    tree = cKDTree(unitvec(lra, ldec))
    dist, idx = tree.query(unitvec(cra, cdec), k=1,
                           distance_upper_bound=chord)
    m = np.isfinite(dist)
    zs = lz[np.where(m, idx, 0)]
    dz = zrm - zs
    good = m & (np.abs(dz) < DZ_CLIP)
    print(f"  matched within {MATCH_ARCSEC}\": {m.sum():,} crayons "
          f"({m.mean():.1%}); {good.sum():,} after |dz|<{DZ_CLIP}")
    nmad = 1.4826 * np.median(np.abs(dz[good] - np.median(dz[good])))
    print(f"  z_RM - z_spec: median {np.median(dz[good]):+.4f}, "
          f"sigma_NMAD {nmad:.4f}, outliers(|dz|>{DZ_CLIP}) "
          f"{1 - good.sum() / max(m.sum(), 1):.2%}")

    # --- per-galaxy pseudo-reddening at FIXED z_spec -------------------
    zsg, cg, vg = zs[good], colors[good], valid[good]
    zc, med_s = m01.locus_table(zsg, cg, vg, nmin=100)
    u = np.full_like(cg, np.nan)
    for c in range(cg.shape[1]):
        ok = np.isfinite(med_s[:, c])
        u[:, c] = cg[:, c] - np.interp(zsg, zc[ok], med_s[ok, c])
    w = np.zeros(cg.shape[1])
    for c in range(cg.shape[1]):
        r = u[vg[:, c], c]
        w[c] = 1.0 / max(1.4826 * np.median(np.abs(r - np.median(r))),
                         1e-3) ** 2
    wv = vg * w                                       # per-galaxy weights
    ehat = (wv * DR * np.where(vg, u, 0)).sum(1) / (wv * DR ** 2).sum(1)

    # --- beta = dz_RM / dE ---------------------------------------------
    dzg = dz[good] - running_median(zsg, dz[good])    # de-trend vs z
    sel = np.abs(ehat) < EHAT_CLIP
    x, y = ehat[sel], dzg[sel]
    beta = np.sum((x - x.mean()) * (y - y.mean())) / np.sum(
        (x - x.mean()) ** 2)
    rng = np.random.default_rng(1)
    bs = [np.sum((xb - xb.mean()) * (yb - yb.mean()))
          / np.sum((xb - xb.mean()) ** 2)
          for _ in range(200)
          for pick in [rng.integers(0, len(x), len(x))]
          for xb, yb in [(x[pick], y[pick])]]
    beta_err = np.std(bs)
    print(f"\n  beta = dz_RM/dE(B-V) = {beta:+.3f} +- {beta_err:.3f} "
          f"(N={sel.sum():,})")

    # --- locus slopes on the FULL sample vs z_RM (what 01 subtracts) ---
    zc_f, med_f = m01.locus_table(zrm, colors, valid)
    Mc, sup = np.zeros(len(CNAMES)), np.zeros(len(CNAMES))
    print(f"\n  {'color':>6s} {'dR_c':>7s} {'L_c(z)':>8s} "
          f"{'M_c':>7s} {'M_c/dR_c':>9s}")
    Lbar = np.zeros(len(CNAMES))
    for c in range(len(CNAMES)):
        ok = np.isfinite(med_f[:, c])
        slope = np.gradient(med_f[ok, c], zc_f[ok])
        Lbar[c] = np.interp(zsg, zc_f[ok], slope).mean()
        Mc[c] = DR[c] - Lbar[c] * beta
        sup[c] = Mc[c] / DR[c]
        print(f"  {CNAMES[c]:>6s} {DR[c]:7.3f} {Lbar[c]:8.3f} "
              f"{Mc[c]:7.3f} {sup[c]:9.3f}")

    # --- overall absorption factor f ------------------------------------
    f = (w * DR * Mc).sum() / (w * DR ** 2).sum()
    ferr = abs((w * DR * Lbar).sum() / (w * DR ** 2).sum()) * beta_err
    print(f"\n  absorption factor f = {f:.3f} +- {ferr:.3f}"
          f"   (E_true = E_measured / f)")

    npz = f"{a.out}/paircount_descrayon_qso.npz"
    if os.path.exists(npz):
        d = np.load(npz)
        rp, ebv, err = d["rp_mid"], d["ebv"], d["ebv_err"]
        print(f"\n  de-absorbed crayon profile ({npz}):")
        print(f"   rp[kpc]   E_meas      E_true = E/f  +- err")
        for j in range(2, 7):
            print(f"  {rp[j]*1e3:8.0f}  {ebv[j]:+9.2e}  {ebv[j]/f:+9.2e}"
                  f"  {err[j]/abs(f):8.1e}")

    os.makedirs(a.out, exist_ok=True)
    fo = f"{a.out}/speczcalib.npz"
    np.savez(fo, beta=beta, beta_err=beta_err, f=f, f_err=ferr, Mc=Mc,
             Lbar=Lbar, dR=DR, w=w, nmatch=good.sum(), nmad=nmad,
             meta=json.dumps(dict(match_arcsec=MATCH_ARCSEC,
                                  date=time.strftime("%F %T"))))
    print(f"saved -> {fo}")


# --------------------------------------------------------------- selftest ---
def selftest():
    """Synthetic crayons with a KNOWN photo-z response: recover beta.

    z_RM is built to respond to the dust-projected part of the color
    noise with beta_true; the pipeline's regression must recover it and
    the implied per-color suppression."""
    rng = np.random.default_rng(11)
    n = 300_000
    beta_true = 0.6
    zs = rng.uniform(0.4, 0.9, n)
    slopes = np.array([1.6, 0.8, 0.4, 2.4])          # locus dL/dz per color
    locus = 0.5 + slopes * (zs[:, None] - 0.4)
    noise = rng.normal(0, 0.04, (n, 4))
    colors = locus + noise
    valid = np.ones((n, 4), bool)
    w = 1.0 / 0.04 ** 2
    ehat_n = (w * DR * noise).sum(1) / (w * DR ** 2).sum()
    zrm = zs + beta_true * ehat_n + rng.normal(0, 0.015, n)

    # regression exactly as in run(): residuals at fixed z_spec
    zc, med = m01.locus_table(zs, colors, valid, nmin=100)
    u = np.column_stack([
        colors[:, c] - np.interp(zs, zc[np.isfinite(med[:, c])],
                                 med[np.isfinite(med[:, c]), c])
        for c in range(4)])
    ehat = (w * DR * u).sum(1) / (w * DR ** 2).sum()
    dz = (zrm - zs) - running_median(zs, zrm - zs)
    sel = np.abs(ehat) < EHAT_CLIP
    x, y = ehat[sel], dz[sel]
    beta = np.sum((x - x.mean()) * (y - y.mean())) / np.sum(
        (x - x.mean()) ** 2)
    print(f"  beta recovered {beta:+.3f} (true {beta_true:+.3f})")
    if abs(beta - beta_true) < 0.1 * beta_true:
        print("SELFTEST PASS")
        return 0
    print("SELFTEST FAIL")
    return 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--mof", action="store_true",
                    help="calibrate on MOF photometry (iteration 3)")
    ap.add_argument("--des-dir", default=m01.DES_DIR_DEFAULT)
    ap.add_argument("--base", default=DESI_BASE)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    sys.exit(selftest() if a.selftest else run(a))
