import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 27_hsc_profile_fig.py — av_profile_hsc: the HSC-SSP DR4 first-wave
# A_V profiles (DESI BGS / LRG foregrounds x quiescent crayons,
# i<24, photo-z 0.7-1.7), anchored cov-ML (10_ml_legacy --hsc, grizy
# responses, all four colors, 10 Mpc anchor).  Style follows
# 18_profile_fig.py (figure of record).  f_abs = 1: the photo-z
# absorption factor is NOT yet calibrated for the HSC crayons (the
# DES value 0.806 does not transfer), so amplitudes are lower
# limits at the ~20% level.  The HSC star nulls (point sources x
# the same foregrounds) are NOT usable as subtraction templates:
# PSF-vs-extended systematics dominate their profiles, so no mlnull
# variant is plotted.
SC = (sys.argv[1] if len(sys.argv) > 1
      else os.path.expanduser("~/data/igmdust/paircount"))
OUT = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
os.makedirs(OUT, exist_ok=True)
RV = 3.1
YMIN, YMAX = 0.01, 100.0
LRG_RPMIN_HSC = 2.0  # Mpc: foreground-light contamination reaches
                     # ~1.2 Mpc in the deep HSC photometry (-13/-20
                     # sigma at 23/41 kpc, -2.8 sigma at 671 kpc)

bgs = np.load(f"{SC}/paircount_hscbgscsfdml_qso.npz")
lrg = np.load(f"{SC}/paircount_hsclrgcsfdml_qso.npz")

fig, ax = plt.subplots(figsize=(10.0, 5.6))


def plot_profile(ax, x, av, s, color, marker, label, xoff=1.0):
    """Detections (>1 sigma) as points; others as 2 sigma upper limits."""
    x = np.asarray(x) * xoff
    det = av > s
    lo = np.clip(np.minimum(s[det], av[det] - YMIN * 1.05), 0, None)
    ax.errorbar(x[det], av[det], yerr=[lo, s[det]], fmt=marker, color=color,
                ms=5, capsize=2, lw=1.2, label=label)
    ul = av[~det] + 2 * s[~det]
    ok = ul > YMIN * 1.3
    ax.errorbar(x[~det][ok], ul[ok], yerr=ul[ok] * 0.22, uplims=True,
                fmt=marker, color=color, ms=3.5, mfc="none", lw=1.0,
                alpha=0.85)


rr = np.logspace(np.log10(0.02), np.log10(30.0), 100)
# Menard+10 Eq. 29, h = 0.7, physical radii (A_V in 1e-3 mag)
ax.plot(rr * 1e3, 4.14 * (rr * 0.7 / 0.1) ** -0.84, "k--", lw=1,
        label="Menard+10 law")

plot_profile(ax, bgs["rp_mid"] * 1e3, RV * bgs["ebv"] * 1e3,
             RV * bgs["ebv_err"] * 1e3, "C3", "o",
             "BGS $\\times$ HSC quiescent crayons (CSFD, cov-ML)")
g = lrg["rp_mid"] >= LRG_RPMIN_HSC
plot_profile(ax, lrg["rp_mid"][g] * 1e3, RV * lrg["ebv"][g] * 1e3,
             RV * lrg["ebv_err"][g] * 1e3, "C1", "D",
             "LRG $\\times$ HSC crayons ($r_p>2$ Mpc, CSFD, cov-ML)",
             xoff=1.07)
ax.set_title("HSC-SSP DR4 imaging (quiescent crayons $i<24$, "
             "$z_{\\rm phot}$ 0.7–1.7, $f_{\\rm abs}=1$, 10 Mpc anchor)",
             fontsize=11)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10, 50000)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=12)
ax.set_ylabel("$A_V$  [$10^{-3}$ mag]", fontsize=12)
ax.legend(fontsize=9.5, loc="lower left", frameon=False)
ax.tick_params(labelsize=11)

# ---- foreground n(z) inset (17_fgz_stats.py) -----------------------------
with open(f"{SC}/fgz_stats.json") as fh:
    FGZ = json.load(fh)
ins = ax.inset_axes([0.615, 0.615, 0.365, 0.345])
for key, color, label in (("bgs", "C3", "BGS"), ("lrg", "C1", "LRG")):
    d = FGZ[key]
    e = np.asarray(d["edges"], float)
    h = np.asarray(d["hist"], float)
    zc = 0.5 * (e[:-1] + e[1:])
    ins.fill_between(zc, h, color=color, alpha=0.16, lw=0)
    ins.plot(zc, h, color=color, lw=1.3,
             label=f"{label} $\\langle z\\rangle{{=}}{d['mean']:.2f}$, "
                   f"$\\tilde z{{=}}{d['median']:.2f}$")
ins.axvspan(0.7, 1.15, color="C2", alpha=0.10, lw=0)
ins.plot([], [], color="C2", lw=5, alpha=0.3,
         label="crayons $z_{\\rm phot}$ 0.7–1.7")
ins.set_xlim(0, 1.15)
ins.set_ylim(bottom=0)
ins.set_yticks([])
ins.set_xlabel("$z$", fontsize=8, labelpad=1)
ins.tick_params(labelsize=7)
ins.legend(fontsize=6.4, frameon=False, loc="upper right",
           handlelength=1.3, borderaxespad=0.1, labelspacing=0.35)
ins.set_title("foreground $n(z)$ (weighted)", fontsize=7.5, pad=2)

fig.tight_layout()
fig.savefig(f"{OUT}/av_profile_hsc.pdf")
fig.savefig(f"{OUT}/av_profile_hsc.png", dpi=150)
print(f"saved -> {OUT}/av_profile_hsc.(pdf,png)")
