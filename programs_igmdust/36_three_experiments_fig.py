"""36_three_experiments_fig.py -- Menard+10, McCleary+25 and the DESI x
Legacy measurements of this work in ONE plot, as published and as
measured, each labelled with the mean foreground redshift it probes
(Nao 2026-11-19: no renormalization -- the three experiments sit at
different redshifts and the plot is there to show that).

Nothing is rescaled.  Published laws are drawn over their fitted
ranges; our points are A_V = 3.1 E(B-V) straight from the cov-ML npz.
All are observer-frame.  The only conversion is h^-1 kpc -> kpc with
h = 0.7 so that the three share a radius axis (the published radii are
taken as physical; McCleary+25 does not say).

  Menard+10     <A_V> = 4.14e-3 (r_p/100 h^-1 kpc)^-0.84, 20 kpc - few
                Mpc; SDSS i<21 galaxies, <z> = 0.36 (number-weighted,
                angular stack) x quasars at z > 1
  McCleary+25   A_V = C (r/100 h^-1 kpc)^alpha  (arXiv:2503.04098)
                WISExSCOS  C = 3.1e-3, alpha = -0.87, zbar = 0.14,
                           fitted 50 h^-1 kpc - 1 h^-1 Mpc
                redMaGiC   C = 4.6e-3, alpha = -0.77, zbar = 0.35,
                           fitted 50 h^-1 kpc - 30 h^-1 Mpc
  this work     BGS x QSOs, photometric r<21 galaxies x QSOs (replica),
                LRG x QSOs (r_p > 0.5 Mpc); lw recounts.

MEAN REDSHIFTS OF OUR MEASUREMENTS.  A stack at fixed PHYSICAL radius
weights each foreground by its pair count, w x N_bg(z_bg > z_fg + 0.2)
/ d_A(z_fg)^2, so the redshift a profile probes is the pair-weighted
mean, lower than the catalog mean (BGS 0.15 vs 0.23).  Both are
computed here from the weighted n(z) of 17_fgz_stats / 19_bgz_stats and
printed; the figure quotes the pair-weighted value.  The replica has no
redshifts: z_eff = 0.30 is assumed by design.

    python3 36_three_experiments_fig.py [NPZDIR]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM

HERE = os.path.dirname(os.path.abspath(__file__))
SP = (sys.argv[1] if len(sys.argv) > 1
      else os.path.expanduser("~/data/igmdust/paircount"))
FIG = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
RV, H, ZBUF = 3.1, 0.7, 0.2
YMIN, YMAX = 0.01, 100.0                 # A_V [1e-3 mag], as 18_profile_fig
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
FGZ = json.load(open(f"{HERE}/fgz_stats_2026-11-12.json"))
BGZ = json.load(open(f"{HERE}/bgz_stats_2026-11-12.json"))


def pair_weighted(fg, bg):
    """(pair-weighted mean z, pair-weighted n(z) on the fg grid)"""
    e = np.asarray(fg["edges"], float)
    zc = 0.5 * (e[:-1] + e[1:])
    n = np.asarray(fg["hist"], float)
    be = np.asarray(bg["edges"], float)
    bh = np.asarray(bg["hist"], float) * np.diff(be)
    above = 1.0 - np.interp(zc + ZBUF, be,
                            np.concatenate([[0], np.cumsum(bh / bh.sum())]))
    w = n * above / COSMO.angular_diameter_distance(zc).value ** 2
    return float(np.sum(w * zc) / w.sum()), zc, w


Z = {}
for k in ("bgs", "lrg"):
    Z[k], zc, w = pair_weighted(FGZ[k], BGZ["qso"])
    print(f"{k.upper()} x QSOs: catalog-weighted <z> = {FGZ[k]['mean']:.3f}, "
          f"pair-weighted <z> = {Z[k]:.3f}")
Z["rep"] = FGZ["photgal"]["zeff"]
print(f"replica: z_eff = {Z['rep']:.2f} (assumed; photometric foregrounds)")
print(f"QSO backgrounds: <z> = {BGZ['qso']['mean']:.2f}")

# one color per measurement, ordered blue -> red with mean redshift
ZCOL = [(0.14, "#1f3fbf"), (0.152, "#1aa0d8"), (0.30, "#2ca02c"),
        (0.35, "#d49100"), (0.36, "#e8590c"), (0.694, "#a50f15")]


def zcol(z):
    return min(ZCOL, key=lambda t: abs(t[0] - z))[1]


def plot_profile(ax, x, av, s, color, marker, label, xoff=1.0):
    """18_profile_fig look: >1 sigma bins as points, the rest as
    2 sigma upper limits (arrows)."""
    x = np.asarray(x) * xoff
    det = av > s
    lo = np.clip(np.minimum(s[det], av[det] - YMIN * 1.05), 0, None)
    ax.errorbar(x[det], av[det], yerr=[lo, s[det]], fmt=marker, color=color,
                ms=6, capsize=2, lw=1.2, label=label, zorder=4)
    ul = av[~det] + 2 * s[~det]
    ok = ul > YMIN * 1.3
    ax.errorbar(x[~det][ok], ul[ok], yerr=ul[ok] * 0.22, uplims=True,
                fmt=marker, color=color, ms=3.5, mfc="none", lw=1.0,
                alpha=0.85, zorder=3)


fig, ax = plt.subplots(figsize=(10.0, 5.6))


def law(C, alpha, rlo, rhi, z, ls, lw, label):
    rr = np.logspace(np.log10(rlo), np.log10(rhi), 100)      # kpc
    ax.plot(rr, C * (rr * H / 100.0) ** alpha, ls=ls, lw=lw, color=zcol(z),
            label=label, zorder=5)


law(4.14, -0.84, 20.0, 5000.0, 0.36, "-", 2.4,
    r"Ménard+10: SDSS $i<21$ galaxies $\times$ QSOs, "
    r"$\langle z\rangle=0.36$")
law(3.1, -0.87, 50.0 / H, 1000.0 / H, 0.14, "--", 2.0,
    r"McCleary+25: WISE$\times$SCOS $\times$ DES crayons, $\bar z=0.14$")
law(4.6, -0.77, 50.0 / H, 30000.0 / H, 0.35, "-.", 2.0,
    r"McCleary+25: redMaGiC LRG $\times$ DES crayons, $\bar z=0.35$")

for key, stem, marker, xoff, rmin, label in (
        ("bgs", "bgs_mqcsfdlwml", "o", 1.0, 0.0,
         r"this work: DESI BGS $\times$ QSOs, "
         r"$\langle z\rangle_{\rm pair}=%.2f$"),
        ("rep", "photgal_mqcsfdlwml", "s", 0.93, 0.0,
         r"this work: Legacy $r<21$ galaxies $\times$ QSOs, "
         r"$z_{\rm eff}=%.2f$ (assumed)"),
        ("lrg", "lrg_mqcsfdlwml", "D", 1.07, 500.0,
         r"this work: DESI LRG $\times$ QSOs ($r_p>0.5$ Mpc), "
         r"$\langle z\rangle_{\rm pair}=%.2f$")):
    d = np.load(f"{SP}/paircount_{stem}_qso.npz", allow_pickle=True)
    rp = d["rp_mid"] * 1e3
    g = rp >= rmin
    plot_profile(ax, rp[g], RV * d["ebv"][g] * 1e3, RV * d["ebv_err"][g] * 1e3,
                 zcol(Z[key]), marker, label % Z[key], xoff=xoff)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(15, 50000)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=13)
ax.set_ylabel("$A_V$  [$10^{-3}$ mag]", fontsize=13)
ax.set_title("Three experiments at three redshifts, as published and as "
             "measured (observer frame, no rescaling)", fontsize=11)
ax.legend(fontsize=8.6, frameon=False, loc="lower left")
ax.tick_params(labelsize=11)

# inset: where each measurement sits in redshift
ins = ax.inset_axes([0.60, 0.60, 0.385, 0.365])
for k, lab, dx, ha in (("bgs", "BGS", 0.02, "left"),
                       ("lrg", "LRG", 0.0, "center")):
    zm, zc, w = pair_weighted(FGZ[k], BGZ["qso"])
    ins.fill_between(zc, w / w.max(), color=zcol(zm), alpha=0.22, lw=0)
    ins.plot(zc, w / w.max(), color=zcol(zm), lw=1.2)
    ins.plot([zm], [1.09], marker="v", ms=6, color=zcol(zm))
    ins.text(zm + dx, 1.50 if k == "bgs" else 1.22, f"{lab} {zm:.2f}",
             ha=ha, va="center", fontsize=6.8, color=zcol(zm))
for z, ls, lab, y, ha in ((0.14, "--", "McC 0.14", 1.50, "right"),
                          (Z["rep"], ":", "replica 0.30", 0.80, "right"),
                          (0.35, "-.", "McC 0.35", 1.22, "left"),
                          (0.36, "-", "M10 0.36", 1.50, "left")):
    ins.axvline(z, ls=ls, lw=1.4, color=zcol(z))
    ins.text(z + (0.015 if ha == "left" else -0.015) + (0.01 if lab.startswith("M10") else 0.0),
             y, lab, ha=ha, va="center", fontsize=6.8, color=zcol(z))
ins.set_xlim(0, 1.12)
ins.set_ylim(0, 1.72)
ins.set_yticks([])
ins.set_xlabel("foreground redshift", fontsize=8, labelpad=1)
ins.tick_params(labelsize=7)
ins.set_title("mean redshift of each measurement "
              "(curves: our pair-weighted $n(z)$)", fontsize=7.5, pad=2)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/av_three_experiments.{ext}", dpi=150)
print(f"saved -> {FIG}/av_three_experiments.pdf")
