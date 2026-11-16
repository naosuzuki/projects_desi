#!/usr/bin/env python3
"""30_bins_analysis.py — igmdust: analyze the 2026-11-15 wave suite.

Loads the cov-ML npz of the five binned studies (brightness split, SFR
split, BGS/LRG z bins, W1-luminosity mass bins), prints profiles, fits
a Menard-template amplitude per sample over its clean radial range, and
writes bins_analysis.json for the figure scripts.  The DES 2+2 bins
(BGS in the DES footprint x redMaGiC crayons, des/08 cov-ML) are fitted
the same way as the cross-survey check; their 125 kpc bin is kept (the
contamination there is a Legacy-photometry effect).

A is defined by E(rp) = A * T(rp), T = 3.7e-3 (rp/0.027 Mpc)^-0.86
(the Menard+10 extrapolation), inverse-variance weighted over the fit
bins; A = 1 means Menard-level reddening.  BGS-type samples fit over
23-383 kpc (the CGM range; the 125 kpc bin is EXCLUDED from split
fits, it is the contamination discriminator), LRG z bins over
0.67-3.6 Mpc (inside 500 kpc LRG is foreground-light dominated).
"""
import json
import os

import numpy as np

SC = os.path.expanduser("~/data/igmdust/paircount")
OUT = os.path.join(SC, "bins_analysis.json")
T0, RP0, SLOPE = 3.7e-3, 0.027, -0.86

META = {  # name -> (label, fit range Mpc, median z, extra)
    "bgsbright":  ("BGS bright half (r<19.12)", (0.02, 0.45), 0.22, ""),
    "bgsfaint":   ("BGS faint half (r>19.12)",  (0.02, 0.45), 0.22, ""),
    "bgssf":      ("BGS star-forming (blue) half", (0.02, 0.45), 0.227, ""),
    "bgspassive": ("BGS passive (red) half",    (0.02, 0.45), 0.228, ""),
    "bgsz0": ("BGS 0.05<z<0.18", (0.02, 0.45), 0.14, "zbin"),
    "bgsz1": ("BGS 0.18<z<0.27", (0.02, 0.45), 0.22, "zbin"),
    "bgsz2": ("BGS 0.27<z<0.45", (0.02, 0.45), 0.33, "zbin"),
    "lrgz0": ("LRG 0.40<z<0.66", (0.6, 4.0), 0.53, "zbin"),
    "lrgz1": ("LRG 0.66<z<0.83", (0.6, 4.0), 0.74, "zbin"),
    "lrgz2": ("LRG 0.83<z<1.10", (0.6, 4.0), 0.93, "zbin"),
    "bgsm0": ("BGS M_W1<-22.4 (massive)", (0.02, 0.45), 0.325, "massbin"),
    "bgsm1": ("BGS -22.4<M_W1<-21.2", (0.02, 0.45), 0.238, "massbin"),
    "bgsm2": ("BGS M_W1>-21.2 (low mass)", (0.02, 0.45), 0.149, "massbin"),
    # DES cross-survey check: BGS inside the DES footprint (350k) x
    # redMaGiC crayons, MOF colors, 2 bins each (des/01 --nzbin/--nmassbin)
    "desz0": ("DES: BGS 0.05<z<0.23", (0.02, 0.45), 0.159, "deszbin"),
    "desz1": ("DES: BGS 0.23<z<0.45", (0.02, 0.45), 0.298, "deszbin"),
    "desm0": ("DES: BGS M_W1<-21.9 (massive)", (0.02, 0.45), 0.297,
              "desmassbin"),
    "desm1": ("DES: BGS M_W1>-21.9 (low mass)", (0.02, 0.45), 0.176,
              "desmassbin"),
}
# npz stem per sample; the DES counts carry their own naming
STEM = {k: f"{k}_mqcsfdlwml" for k in META if not k.startswith("des")}
STEM.update({f"des{b}": f"descrayon{b}mofcsfdml"
             for b in ("z0", "z1", "m0", "m1")})
# the 125 kpc bin is the diagnosed bright-foreground light-contamination
# bin at BGS brightness (confirmed by the split: -3.6 sigma bright half,
# +0.1 faint) -- excluded from every BGS-type fit
EXCL_125 = {k for k in META if k.startswith("bgs")}


def template(rp):
    return T0 * (rp / RP0) ** SLOPE


res = {}
for name, (label, (rlo, rhi), zmed, kind) in META.items():
    d = np.load(f"{SC}/paircount_{STEM[name]}_qso.npz")
    rp, e, s = d["rp_mid"], d["ebv"], d["ebv_err"]
    m = (rp >= rlo) & (rp <= rhi) & (s > 0)
    if name in EXCL_125:
        m &= ~np.isclose(rp, 0.125, atol=0.02)
    T = template(rp[m])
    w = T ** 2 / s[m] ** 2
    A = float(np.sum(e[m] * T / s[m] ** 2) / np.sum(w))
    sA = float(1.0 / np.sqrt(np.sum(w)))
    res[name] = dict(label=label, zmed=zmed, kind=kind, A=A, sA=sA,
                     rp=rp.tolist(), ebv=e.tolist(), err=s.tolist())
    print(f"== {name:11s} {label}")
    print(f"   A = {A:+.3f} +- {sA:.3f} x Menard "
          f"(fit {rlo*1e3:.0f}-{rhi*1e3:.0f} kpc, "
          f"{A/sA:+.1f} sigma), median z = {zmed}")
    for r, ee, ss in zip(rp, e, s):
        if 0.02 <= r <= 1.3:
            print(f"   {r*1e3:7.0f} kpc  E = {ee:+.2e} +- {ss:.1e} "
                  f"({ee/ss:+5.1f} sig)")

with open(OUT, "w") as fh:
    json.dump(res, fh)
print(f"\nsaved -> {OUT}")
