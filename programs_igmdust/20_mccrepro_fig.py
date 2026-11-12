"""20_mccrepro_fig.py -- figure script of record for the McCleary
reproduction exhibit (ms.tex fig:mccrepro): DES iteration-4 profile
(MOF + CSFD + covariance-ML) against her published LRG power law.

Versions (Nao, 2026-11-12: try insets and separate n(z) as different
version numbers):
  v2  mccleary_repro_v2.pdf/png -- profile with the redshift inset:
      weighted n(z) of the redMaGiC high-density foregrounds
      (17_fgz_stats.py) and the redMaGiC highlum_highz crayon
      backgrounds (19_bgz_stats.py), annotated with mean and median
      z.  No +/- sign coding: every kept bin is one filled symbol
      at |A_V| (the 2026-11-08..11 open-symbol convention for
      negative bins is dropped; no kept bin is negative anyway).
  v3  mccleary_repro_v3.pdf/png (the same profile without the
      inset) plus mccleary_repro_v3_nz.pdf/png, the redshift
      distribution as its own standalone figure.
  --v1 reproduces the original 2026-11-08 figure byte-for-byte
      (recovered from the session transcript after the scratchpad
      wipe).

Inputs: paircount_rmgmofcsfdml_qso.npz in NPZDIR (argv[1], default
~/data/igmdust/paircount, local mirror of
astro4:~/data/igmdust/paircount/), plus fgz_stats_2026-11-12.json /
bgz_stats_2026-11-12.json in the repo.

    python3 20_mccrepro_fig.py [NPZDIR]        # v2 + v3 figures
    python3 20_mccrepro_fig.py [NPZDIR] --v1   # original figure
"""
import argparse
import json
import os

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
SP = a.npzdir or os.path.expanduser("~/data/igmdust/paircount")
FIG = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
F_ABS, RV = 0.806, 3.1


def law_av(r_kpc):                      # McCleary+25 LRG power law
    return 4.6e-3 * (0.91 * r_kpc / 100.0) ** -0.77


d = np.load(f"{SP}/paircount_rmgmofcsfdml_qso.npz", allow_pickle=True)
rp = d["rp_mid"] * 1e3                  # kpc
av = RV * d["ebv"] / F_ABS
err = RV * d["ebv_err"] / F_ABS
keep = rp >= 150


def draw_profile(ax, signed, text_x):
    rr = np.logspace(np.log10(140), np.log10(40000), 200)
    ax.plot(rr, law_av(rr), "k--", lw=1.6, zorder=3,
            label=r"McCleary et al. (2026) LRG law:"
                  r" $A_V = 4.6\times10^{-3}\,(r/100\,h^{-1}{\rm kpc})"
                  r"^{-0.77}$")
    ax.fill_between(rr, law_av(rr) * 3.1 / 4.6, law_av(rr), color="0.82",
                    zorder=1, label=r"their redMaGiC range $C=(3.1$--$4.6)"
                                    r"\times10^{-3}$")
    y = av[keep] if signed else np.abs(av[keep])
    ax.errorbar(rp[keep], y, yerr=err[keep], fmt="o", ms=7, color="C3",
                capsize=3, zorder=4,
                label=r"this work: redMaGiC $\times$ DES crayons"
                      r" (MOF + CSFD + cov-ML),"
                      r" $A_V = 3.1\,{\rm E}/f_{\rm abs}$")
    if signed:                          # v1: open symbols at |A_V| for E<0
        neg = av[keep] < 0
        if neg.any():
            ax.errorbar(rp[keep][neg], -av[keep][neg], yerr=err[keep][neg],
                        fmt="o", ms=7, mfc="none", color="C3", capsize=3,
                        zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(140, 45000)
    ax.set_ylim(2e-5, 8e-3)
    ax.set_xlabel(r"$r_p$ [physical kpc]", fontsize=13)
    ax.set_ylabel(r"$A_V$ [mag]", fontsize=13)
    ax.text(text_x, 0.08,
            "mid zone (0.15–1.5 Mpc): mean pull +0.24, $\\chi^2/n=1.00$\n"
            "tail zone (1.5–34 Mpc): mean pull +0.78, $\\chi^2/n=0.71$\n"
            "MCCLEARY_REPRODUCTION_PASS (2026-11-08)",
            transform=ax.transAxes, fontsize=10.5, va="bottom",
            bbox=dict(fc="white", ec="0.6", alpha=0.9))
    ax.legend(loc="upper right", fontsize=9.5)
    ax.tick_params(labelsize=11)


if a.v1:
    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    draw_profile(ax, signed=True, text_x=0.03)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/mccleary_repro.{ext}", dpi=150)
    print("pulls:", np.round((av[keep] - law_av(rp[keep])) / err[keep], 2))
    print(f"saved -> {FIG}/mccleary_repro.pdf")
    raise SystemExit

# ---------------- v2 + v3 ----------------
with open(f"{HERE}/fgz_stats_2026-11-12.json") as fh:
    FGZ = json.load(fh)
with open(f"{HERE}/bgz_stats_2026-11-12.json") as fh:
    BGZ = json.load(fh)
NZ = ((FGZ["redmagic"], "C4", "fg: redMaGiC highdens"),
      (BGZ["descrayon"], "C0", "bg: highlum crayons"))


def nz_curves(ax, lw, fontsize, two_line):
    for d_, color, label in NZ:
        e = np.asarray(d_["edges"], float)
        h = np.asarray(d_["hist"], float)
        h = h / h.max()                 # peak-normalized shapes
        zc = 0.5 * (e[:-1] + e[1:])
        sep = "\n" if two_line else " "
        ax.fill_between(zc, h, color=color, alpha=0.16, lw=0)
        ax.plot(zc, h, color=color, lw=lw,
                label=f"{label}{sep}$\\langle z\\rangle{{=}}"
                      f"{d_['mean']:.2f}$, $\\tilde z{{=}}"
                      f"{d_['median']:.2f}$")
    ax.set_xlim(0, 1.0)
    ax.set_yticks([])


# ---- v2: profile with the n(z) inset ----
fig, ax = plt.subplots(figsize=(10.0, 5.6))
draw_profile(ax, signed=False, text_x=0.345)
ins = ax.inset_axes([0.04, 0.07, 0.27, 0.33])
nz_curves(ins, lw=1.3, fontsize=6.6, two_line=True)
ins.set_ylim(0, 1.55)
ins.set_xlabel("$z$", fontsize=8, labelpad=1)
ins.tick_params(labelsize=7)
ins.legend(fontsize=6.6, frameon=False, loc="upper left",
           handlelength=1.2, borderaxespad=0.1, labelspacing=0.45)
ins.set_title("foreground / background $n(z)$", fontsize=7.5, pad=2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_repro_v2.{ext}", dpi=150)
print(f"saved -> {FIG}/mccleary_repro_v2.pdf")
plt.close(fig)

# ---- v3: profile without the inset ----
fig, ax = plt.subplots(figsize=(10.0, 5.6))
draw_profile(ax, signed=False, text_x=0.03)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_repro_v3.{ext}", dpi=150)
print(f"saved -> {FIG}/mccleary_repro_v3.pdf")
plt.close(fig)

# ---- v3 companion: standalone n(z) figure ----
fig, ax = plt.subplots(figsize=(6.2, 4.0))
nz_curves(ax, lw=1.6, fontsize=9.5, two_line=False)
ax.set_ylim(0, 1.45)
ax.set_xlabel("$z$", fontsize=12)
ax.set_ylabel("$n(z)$  [peak-normalized]", fontsize=12)
ax.legend(fontsize=9.5, frameon=False, loc="upper left")
ax.tick_params(labelsize=11)
ax.set_title("Redshift distributions of the McCleary-reproduction "
             "samples", fontsize=10.5)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_repro_v3_nz.{ext}", dpi=150)
print("pulls:", np.round((av[keep] - law_av(rp[keep])) / err[keep], 2))
print(f"saved -> {FIG}/mccleary_repro_v3_nz.pdf")
