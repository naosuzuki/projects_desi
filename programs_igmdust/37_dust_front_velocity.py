"""37_dust_front_velocity.py -- does the tilt of the profiles measure a
dust propagation speed?  (Nao 2026-11-19: "it takes time for the dust to
reach that distance; from these tilts we can deduce the velocity".)

For every lw profile of this work (anchored at >= 10 Mpc, which the
models reproduce by subtracting their own anchor-zone mean):
  * a free power-law slope alpha, overall, inside 1 Mpc and over the
    common megaparsec range 0.6-10 Mpc;
  * an expanding-front model: a steady constant-speed outflow has
    rho ~ r^-2 out to the front R = v (t - t_on) and nothing beyond, so
        Sigma(r_p) ~ (1/r_p) arctan( sqrt(R^2 - r_p^2) / r_p ),
    slope -1 well inside the front and steepening toward it.  R is
    profiled on a grid with the amplitude solved linearly.
  * speeds: the minimum mean speed R / age(z) (dust launched at t = 0)
    and R / (age(z) - age(z=2)) (launched at cosmic noon).

    python3 37_dust_front_velocity.py [NPZDIR]
"""
import os
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
SC = (sys.argv[1] if len(sys.argv) > 1
      else os.path.expanduser("~/data/igmdust/paircount"))
KMS = 977.8                                   # 1 Mpc/Gyr in km/s
S = {  # name: (npz stem, mean z probed, r_min Mpc, drop the 125 kpc bin)
    "BGS all": ("bgs_mqcsfdlwml", 0.15, 0.02, True),
    "BGS z0": ("bgsz0_mqcsfdlwml", 0.14, 0.02, True),
    "BGS z1": ("bgsz1_mqcsfdlwml", 0.22, 0.02, True),
    "BGS z2": ("bgsz2_mqcsfdlwml", 0.33, 0.02, True),
    "PathB": ("bgsmcc_mqcsfdlwml", 0.25, 0.02, True),
    "replica": ("photgal_mqcsfdlwml", 0.30, 0.02, False),
    "LRG all": ("lrg_mqcsfdlwml", 0.69, 0.6, False),
    "LRG z0": ("lrgz0_mqcsfdlwml", 0.53, 0.6, False),
    "LRG z1": ("lrgz1_mqcsfdlwml", 0.74, 0.6, False),
    "LRG z2": ("lrgz2_mqcsfdlwml", 0.93, 0.6, False)}
ALPHAS = np.linspace(-2.5, 0.0, 251)
RS = np.logspace(np.log10(0.05), np.log10(300), 160)


def load(stem, rmin, ex, rmax=10.0):
    d = np.load(f"{SC}/paircount_{stem}_qso.npz", allow_pickle=True)
    rp, e, s = d["rp_mid"], d["ebv"] * 1e3, d["ebv_err"] * 1e3
    k = (rp >= rmin) & (rp < rmax) & (s > 0)
    if ex:
        k &= ~np.isclose(rp, 0.125, atol=0.02)
    return rp, e, s, k, rp >= 10


def chi2(rp, e, s, k, an, shape):
    m = shape(rp)
    w = 1.0 / s[an] ** 2
    m = m - np.sum(m[an] * w) / w.sum()          # the 10 Mpc anchor
    den = np.sum(m[k] ** 2 / s[k] ** 2)
    if not np.isfinite(den) or den <= 0:
        return np.sum((e[k] / s[k]) ** 2), 0.0, np.inf
    n = np.sum(e[k] * m[k] / s[k] ** 2) / den
    return np.sum(((e[k] - n * m[k]) / s[k]) ** 2), n, 1 / np.sqrt(den)


def front(R):
    def f(rp):
        out = np.zeros_like(rp)
        ins = rp < R
        out[ins] = np.arctan(np.sqrt(R * R - rp[ins] ** 2) / rp[ins]) / rp[ins]
        return out
    return f


def slope(name, rmin=None, rmax=10.0, drop=None):
    stem, z, r0, ex = S[name]
    rp, e, s, k, an = load(stem, rmin or r0, ex if drop is None else drop,
                           rmax)
    c = np.array([chi2(rp, e, s, k, an, lambda r, a=a: r ** a)[0]
                  for a in ALPHAS])
    ok = ALPHAS[c <= c.min() + 1]
    return ALPHAS[c.argmin()], ok.min(), ok.max()


print("== power-law slopes (1 sigma range; the grid stops at -2.5 and 0)")
print("   sample    <z>   all radii            inside 1 Mpc         "
      "0.6-10 Mpc")
for name, (stem, z, r0, ex) in S.items():
    a = slope(name)
    b = slope(name, rmax=1.0) if r0 < 0.5 else None
    c = slope(name, rmin=0.6, drop=False)
    fmt = lambda t: (f"{t[0]:+.2f} ({t[1]:+.2f},{t[2]:+.2f})"
                     if t else "       --          ")
    print(f"   {name:8s} {z:.2f}  {fmt(a)}  {fmt(b)}  {fmt(c)}")

print("\n== expanding front (rho ~ r^-2 inside R): front radius and speed")
print("   sample    <z>  age   R best  (1 sig)         2 sig lower  "
      "dchi2(no front)   v_min since t=0 / since z=2 [km/s], 2 sig")
for name, (stem, z, r0, ex) in S.items():
    rp, e, s, k, an = load(stem, r0, ex)
    c = np.array([chi2(rp, e, s, k, an, front(R))[0] for R in RS])
    cinf = chi2(rp, e, s, k, an, lambda r: 1 / r)[0]
    ref = min(c.min(), cinf)
    r1, r2 = RS[c <= ref + 1], RS[c <= ref + 4]
    hi = f"{r1.max():6.2f}" if r1.max() < RS.max() * 0.99 else "   inf"
    age = COSMO.age(z).value
    t2 = age - COSMO.age(2.0).value
    print(f"   {name:8s} {z:.2f} {age:5.1f}  {RS[c.argmin()]:6.2f} "
          f"({r1.min():6.2f},{hi})   > {r2.min():5.2f} Mpc   "
          f"{cinf - c.min():+5.1f}            > {r2.min()/age*KMS:5.0f} / "
          f"> {r2.min()/t2*KMS:5.0f}")
print("\n   to reach 0.3 / 1 / 3 / 7 Mpc in 5 Gyr: "
      + " / ".join(f"{R/5*KMS:.0f}" for R in (0.3, 1, 3, 7)) + " km/s")
