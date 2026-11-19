"""36_zevolution_fig.py -- one plot of the halo-reddening amplitude versus
foreground redshift: Menard+10, McCleary+25 and the DESI x Legacy
measurements of this work (Nao 2026-11-19: the three experiments probe
different redshifts).

Common quantity: observer-frame A_V at r_p = 100 h^-1 kpc (143 kpc
physical for h = 0.7), the pivot at which both published laws are
normalized, so the literature points need no slope extrapolation.

  Menard+10     <A_V> = (4.14 +- 0.19)e-3 (r_p/100 h^-1 kpc)^-0.84,
                observer-frame V, i < 21 photometric galaxies,
                <z> = 0.36 (their Eqs. 21, 29).
  McCleary+25   A_V = C (r/100 h^-1 kpc)^alpha (arXiv:2503.04098):
                WISExSCOS  C = (3.1 +- 0.06)e-3, alpha = -0.87, zbar = 0.14
                redMaGiC   C = (4.6 +- 0.2)e-3,  alpha = -0.77, zbar = 0.35
                plotted as published; the paper does not say whether its
                radii are physical or comoving (comoving would lower the
                points by 11% and 21% at this pivot).
  this work     A_V = 3.1 A T(143 kpc), with A the template amplitude of
                30_bins_analysis.py and T = 1.2e-3 (r_p/100 kpc)^-0.86 in
                E(B-V): BGS z bins (fit 20-450 kpc), LRG z bins (fit
                0.6-4 Mpc, i.e. EXTRAPOLATED inward to the pivot: open
                symbols), the photometric replica (the Menard analog,
                r < 21 galaxies, z_eff = 0.30) and Path B (the McCleary
                analog: BGS passive, M_r < -20.2; her slope).

NOTE on the template: T(143 kpc) x 3.1 = 2.73e-3, i.e. 0.66 of the
published Menard+10 value at the same radius.  "A = 1" in the other
figures is therefore 0.66 x Menard+10, not Menard+10.

    python3 36_zevolution_fig.py [NPZDIR]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SP = (sys.argv[1] if len(sys.argv) > 1
      else os.path.expanduser("~/data/igmdust/paircount"))
FIG = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
RV, H = 3.1, 0.7
RPIV = 100.0 / H                               # kpc physical
T_PIV = 1.2e-3 * (RPIV / 100.0) ** -0.86       # template E(B-V) at pivot
K = RV * T_PIV * 1e3                           # A -> A_V [1e-3 mag]

R = json.load(open(f"{SP}/bins_analysis.json"))
M = json.load(open(f"{SP}/mccleary_legacy.json"))


def series(keys, edges):
    z = np.array([R[k]["zmed"] for k in keys])
    a = np.array([R[k]["A"] for k in keys]) * K
    s = np.array([R[k]["sA"] for k in keys]) * K
    lo = z - np.array(edges[:-1])
    hi = np.array(edges[1:]) - z
    return z, a, s, [lo, hi]


# the replica: same fit as 30_ over 20-450 kpc (no bin excluded)
d = np.load(f"{SP}/paircount_photgal_mqcsfdlwml_qso.npz", allow_pickle=True)
rp, e, s = d["rp_mid"], d["ebv"], d["ebv_err"]
m = (rp >= 0.02) & (rp <= 0.45)
T = 1.2e-3 * (rp[m] / 0.1) ** -0.86
A_rep = np.sum(e[m] * T / s[m] ** 2) / np.sum(T ** 2 / s[m] ** 2)
sA_rep = 1.0 / np.sqrt(np.sum(T ** 2 / s[m] ** 2))

# Path B: best-fit profile C * (her law in the comoving reading used by
# 34_) evaluated at the pivot -- independent of that reading
law_piv = 4.6e-3 * (0.91 * RPIV / 100.0) ** -0.77 * 1e3
B, sB = M["bgsmcc"]["C"] * law_piv, M["bgsmcc"]["sC"] * law_piv

print(f"pivot r_p = {RPIV:.0f} kpc physical;  A -> A_V factor {K:.2f}e-3"
      f"  (template/Menard+10 at pivot = {K/4.14:.2f})")

fig, ax = plt.subplots(figsize=(8.0, 5.2))
zz = np.linspace(0.0, 1.05, 100)
ax.plot(zz, 4.14 * ((1 + zz) / 1.36) ** 1.1, color="0.55", ls=":", lw=1.2,
        label=r"no evolution: fixed dust column seen in the observed "
              r"bands, $\propto(1+z)^{1.1}$")
ax.axhline(0, color="0.7", lw=0.8)

ax.errorbar([0.36], [4.14], yerr=[0.19], fmt="*", ms=17, color="k",
            capsize=3, zorder=6,
            label=r"M$\acute{\rm e}$nard+10: SDSS $i<21$ galaxies "
                  r"$\times$ QSOs")
ax.errorbar([0.14, 0.35], [3.1, 4.6], yerr=[0.06, 0.2],
            xerr=[[0.14, 0.20], [0.06, 0.10]], fmt="^", ms=10, color="C4",
            capsize=3, zorder=5,
            label=r"McCleary+25: WISE$\times$SCOS ($\bar z=0.14$), "
                  r"redMaGiC LRG ($\bar z=0.35$) $\times$ DES crayons")

z, a, s_, xe = series(["bgsz0", "bgsz1", "bgsz2"],
                      [0.05, 0.179, 0.272, 0.45])
ax.errorbar(z, a, yerr=s_, xerr=xe, fmt="o", ms=7, color="C3", capsize=3,
            zorder=4, label=r"this work: DESI BGS $\times$ QSOs, 3 $z$ "
                            "bins (fit 20–450 kpc)")
print("BGS :", " ".join(f"z={x:.2f} {y:+.2f}+-{w:.2f}"
                        for x, y, w in zip(z, a, s_)))
z, a, s_, xe = series(["lrgz0", "lrgz1", "lrgz2"],
                      [0.40, 0.661, 0.830, 1.10])
ax.errorbar(z, a, yerr=s_, xerr=xe, fmt="D", ms=7, color="C1", mfc="w",
            capsize=3, zorder=4,
            label=r"this work: DESI LRG $\times$ QSOs, 3 $z$ bins "
                  "(fit 0.6–4 Mpc, extrapolated to the pivot)")
print("LRG :", " ".join(f"z={x:.2f} {y:+.2f}+-{w:.2f}"
                        for x, y, w in zip(z, a, s_)))
ax.errorbar([0.30], [A_rep * K], yerr=[sA_rep * K], fmt="s", ms=8,
            color="C0", capsize=3, zorder=5,
            label=r"this work: Legacy $r<21$ galaxies $\times$ QSOs "
                  r"(M$\acute{\rm e}$nard analog)")
ax.errorbar([0.27], [B], yerr=[sB], xerr=[[0.12], [0.18]], fmt="P", ms=9,
            color="C2", capsize=3, zorder=5,
            label=r"this work: BGS passive, $M_r<-20.2$ $\times$ QSOs "
                  "(McCleary analog)")
print(f"replica: z=0.30 {A_rep*K:+.2f}+-{sA_rep*K:.2f}   "
      f"(A = {A_rep:.2f}+-{sA_rep:.2f})")
print(f"Path B : z=0.27 {B:+.2f}+-{sB:.2f}")

ax.set_xlim(0, 1.12)
ax.set_ylim(-1.5, 15.5)
ax.set_xlabel("foreground redshift", fontsize=12)
ax.set_ylabel(r"$A_V$ at $r_p=100\,h^{-1}$ kpc  [$10^{-3}$ mag, "
              "observer frame]", fontsize=11.5)
ax.set_title("Halo reddening amplitude versus foreground redshift: "
             "three experiments", fontsize=10.5)
ax.legend(fontsize=7.6, frameon=False, loc="upper left")
ax.tick_params(labelsize=11)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/av_of_z_experiments.{ext}", dpi=150)
print(f"saved -> {FIG}/av_of_z_experiments.pdf")
