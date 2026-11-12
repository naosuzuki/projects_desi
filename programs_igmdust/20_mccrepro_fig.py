"""20_mccrepro_fig.py -- figure script of record for the McCleary
reproduction exhibit (ms.tex fig:mccrepro): DES iteration-4 profile
(MOF + CSFD + covariance-ML) against her published LRG power law.

v2 (default, 2026-11-12, -> mccleary_repro_v2.pdf/png): adds the
redshift inset -- weighted n(z) of the redMaGiC high-density
foregrounds (17_fgz_stats.py) and the redMaGiC highlum_highz crayon
backgrounds (19_bgz_stats.py), each annotated with its mean and
median redshift.  --v1 reproduces the original 2026-11-08 figure
byte-for-byte (recovered from the session transcript after the
scratchpad wipe; it was never a repo file before).

Inputs: paircount_rmgmofcsfdml_qso.npz in NPZDIR (argv[1], default
script dir; master copy astro4:~/data/igmdust/paircount/), plus
fgz_stats_2026-11-12.json / bgz_stats_2026-11-12.json in the repo.

    python3 20_mccrepro_fig.py [NPZDIR]        # v2 figure
    python3 20_mccrepro_fig.py [NPZDIR] --v1   # original figure
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("npzdir", nargs="?", default=None)
ap.add_argument("--v1", action="store_true",
                help="reproduce the original 2026-11-08 figure")
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))
SP = a.npzdir or HERE
FIG = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
F_ABS, RV = 0.806, 3.1


def law_av(r_kpc):                      # McCleary+25 LRG power law
    return 4.6e-3 * (0.91 * r_kpc / 100.0) ** -0.77


d = np.load(f"{SP}/paircount_rmgmofcsfdml_qso.npz", allow_pickle=True)
rp = d["rp_mid"] * 1e3                  # kpc
av = RV * d["ebv"] / F_ABS
err = RV * d["ebv_err"] / F_ABS
keep = rp >= 150

fig, ax = plt.subplots(figsize=(10.0, 5.6))
rr = np.logspace(np.log10(140), np.log10(40000), 200)
ax.plot(rr, law_av(rr), "k--", lw=1.6, zorder=3,
        label=r"McCleary et al. (2026) LRG law:"
              r" $A_V = 4.6\times10^{-3}\,(r/100\,h^{-1}{\rm kpc})^{-0.77}$")
ax.fill_between(rr, law_av(rr) * 3.1 / 4.6, law_av(rr), color="0.82",
                zorder=1, label=r"their redMaGiC range $C=(3.1$--$4.6)"
                                r"\times10^{-3}$")
ax.errorbar(rp[keep], av[keep], yerr=err[keep], fmt="o", ms=7, color="C3",
            capsize=3, zorder=4,
            label=r"this work: redMaGiC $\times$ DES crayons"
                  r" (MOF + CSFD + cov-ML), $A_V = 3.1\,{\rm E}/f_{\rm abs}$")
neg = av[keep] < 0
if neg.any():
    ax.errorbar(rp[keep][neg], -av[keep][neg], yerr=err[keep][neg],
                fmt="o", ms=7, mfc="none", color="C3", capsize=3, zorder=4)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(140, 45000)
ax.set_ylim(2e-5, 8e-3)
ax.set_xlabel(r"$r_p$ [physical kpc]", fontsize=13)
ax.set_ylabel(r"$A_V$ [mag]", fontsize=13)
ax.text(0.03 if a.v1 else 0.345, 0.08,
        "mid zone (0.15–1.5 Mpc): mean pull +0.24, $\\chi^2/n=1.00$\n"
        "tail zone (1.5–34 Mpc): mean pull +0.78, $\\chi^2/n=0.71$\n"
        "MCCLEARY_REPRODUCTION_PASS (2026-11-08)",
        transform=ax.transAxes, fontsize=10.5, va="bottom",
        bbox=dict(fc="white", ec="0.6", alpha=0.9))
ax.legend(loc="upper right", fontsize=9.5)
ax.tick_params(labelsize=11)

suffix = ""
if not a.v1:                            # ---- v2: fg + bg n(z) inset ----
    suffix = "_v2"
    with open(f"{HERE}/fgz_stats_2026-11-12.json") as fh:
        FGZ = json.load(fh)
    with open(f"{HERE}/bgz_stats_2026-11-12.json") as fh:
        BGZ = json.load(fh)
    ins = ax.inset_axes([0.04, 0.07, 0.27, 0.33])
    for d_, color, label in ((FGZ["redmagic"], "C4",
                              "fg: redMaGiC highdens"),
                             (BGZ["descrayon"], "C0",
                              "bg: highlum crayons")):
        e = np.asarray(d_["edges"], float)
        h = np.asarray(d_["hist"], float)
        h = h / h.max()                 # peak-normalized shapes
        zc = 0.5 * (e[:-1] + e[1:])
        ins.fill_between(zc, h, color=color, alpha=0.16, lw=0)
        ins.plot(zc, h, color=color, lw=1.3,
                 label=f"{label}\n$\\langle z\\rangle{{=}}{d_['mean']:.2f}$,"
                       f" $\\tilde z{{=}}{d_['median']:.2f}$")
    ins.set_xlim(0, 1.0)
    ins.set_ylim(0, 1.55)
    ins.set_yticks([])
    ins.set_xlabel("$z$", fontsize=8, labelpad=1)
    ins.tick_params(labelsize=7)
    ins.legend(fontsize=6.6, frameon=False, loc="upper left",
               handlelength=1.2, borderaxespad=0.1, labelspacing=0.45)
    ins.set_title("foreground / background $n(z)$", fontsize=7.5, pad=2)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_repro{suffix}.{ext}", dpi=150)
print("pulls:", np.round((av[keep] - law_av(rp[keep])) / err[keep], 2))
print(f"saved -> {FIG}/mccleary_repro{suffix}.pdf")
