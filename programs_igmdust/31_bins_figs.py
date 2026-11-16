#!/usr/bin/env python3
"""31_bins_figs.py — figures for the 2026-11-15 binned-wave results.

Reads bins_analysis.json (30_bins_analysis.py).  One plot per figure
(Nao 2026-11-12): ebv_split_bgs, ebv_sfr_bgs, a_of_z, a_of_m.
Profiles follow the house style of 21_crayon_deficit_fig.py: log-log
|E(B-V)|, filled symbols for >=1 sigma bins, 2 sigma upper-limit arrows
otherwise, n(z) inset; signs go in the caption.  The two A-plots are
linear and carry the DES 2+2 bins as the cross-survey check."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SC = os.path.expanduser("~/data/igmdust/paircount")
OUT = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
os.makedirs(OUT, exist_ok=True)
R = json.load(open(f"{SC}/bins_analysis.json"))
T0, RP0, SLOPE = 3.7e-3, 0.027, -0.86


YMIN, YMAX = 8e-3, 40.0            # |E(B-V)| in 1e-3 mag (21_ house style)
HERE = os.path.dirname(os.path.abspath(__file__))
FGZ = json.load(open(f"{HERE}/fgz_stats_2026-11-12.json"))
BGZ = json.load(open(f"{HERE}/bgz_stats_2026-11-12.json"))


def plot_abs(ax, x, e, s, color, marker, label):
    """v2/v3 look (21_crayon_deficit_fig): one filled symbol per
    >=1 sigma bin at |E|, no sign coding (signs stay in the caption);
    <1 sigma bins are 2 sigma upper limits (arrows)."""
    x, e, s = np.asarray(x), np.asarray(e), np.asarray(s)
    det = np.abs(e) > s
    a_ = np.abs(e[det])
    lo = np.clip(np.minimum(s[det], a_ - YMIN * 1.05), 0, None)
    ax.errorbar(x[det], a_, yerr=[lo, s[det]], fmt=marker, color=color,
                ms=5.5, capsize=2, lw=1.2, label=label)
    ul = np.abs(e[~det]) + 2 * s[~det]
    ok = ul > YMIN * 1.3
    ax.errorbar(x[~det][ok], ul[ok], yerr=ul[ok] * 0.22, uplims=True,
                fmt=marker, color=color, ms=3.0, mfc="none", lw=0.9,
                alpha=0.6)


def nz_curve(ax, d_, color, label):
    e = np.asarray(d_["edges"], float)
    h = np.asarray(d_["hist"], float)
    h = h / h.max()
    zc = 0.5 * (e[:-1] + e[1:])
    ax.fill_between(zc, h, color=color, alpha=0.16, lw=0)
    ax.plot(zc, h, color=color, lw=1.2,
            label=label + (f" $\\langle z\\rangle{{=}}{d_['mean']:.2f}$,"
                           f" $\\tilde z{{=}}{d_['median']:.2f}$"))


def nz_inset(ax, entries, zmax):
    ins = ax.inset_axes([0.615, 0.615, 0.365, 0.345])
    for d_, color, label in entries:
        nz_curve(ins, d_, color, label)
    ins.set_xlim(0, zmax)
    ins.set_ylim(0, 1.6)
    ins.set_yticks([])
    ins.set_xlabel("$z$", fontsize=8, labelpad=1)
    ins.tick_params(labelsize=7)
    ins.legend(fontsize=6.2, frameon=False, loc="upper right",
               handlelength=1.1, borderaxespad=0.1, labelspacing=0.35)
    ins.set_title("fg / bg $n(z)$", fontsize=7.5, pad=2)


def profile_fig(pairs, fname, title):
    fig, ax = plt.subplots(figsize=(6.2, 4.9))
    rr = np.logspace(np.log10(0.012), np.log10(40.0), 150)
    ax.plot(rr * 1e3, T0 * (rr / RP0) ** SLOPE * 1e3, "k--", lw=1,
            label="Ménard+10 extrapolation")
    for name, color, marker, xoff in pairs:
        d = R[name]
        rp = np.asarray(d["rp"]) * 1e3 * xoff
        e = np.asarray(d["ebv"]) * 1e3
        s = np.asarray(d["err"]) * 1e3
        plot_abs(ax, rp, e, s, color, marker,
                 f"{d['label']}  ($A={d['A']:.2f}\\pm{d['sA']:.2f}$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(12, 46000)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=12)
    ax.set_ylabel("$|{\\rm E(B\\!-\\!V)}|$  [$10^{-3}$ mag]", fontsize=12)
    ax.set_title(title, fontsize=10.5)
    ax.legend(fontsize=8.5, frameon=False, loc="lower left")
    ax.tick_params(labelsize=11)
    nz_inset(ax, [(FGZ["bgs"], "C2", "fg: BGS"),
                  (BGZ["qso"], "C3", "bg: QSOs")], 4.0)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{fname}.pdf")
    fig.savefig(f"{OUT}/{fname}.png", dpi=150)
    print(f"saved -> {OUT}/{fname}.(pdf,png)")


profile_fig([("bgsbright", "C3", "o", 1.0), ("bgsfaint", "C0", "s", 1.07)],
            "ebv_split_bgs",
            "BGS brightness split $\\times$ merged QSOs "
            "(lw, CSFD, cov-ML, 10 Mpc anchor)")
profile_fig([("bgspassive", "C3", "o", 1.0), ("bgssf", "C0", "s", 1.07)],
            "ebv_sfr_bgs",
            "BGS star-formation split (z-conditional $g-r$) "
            "$\\times$ merged QSOs")

# ---- A(z): six points, BGS and LRG --------------------------------------
fig, ax = plt.subplots(figsize=(8.0, 4.8))
for keys, color, marker, label in (
        (["bgsz0", "bgsz1", "bgsz2"], "C3", "o",
         "BGS $\\times$ QSOs (fit 20–450 kpc)"),
        (["lrgz0", "lrgz1", "lrgz2"], "C1", "D",
         "LRG $\\times$ QSOs (fit 0.6–4 Mpc)"),
        (["desz0", "desz1"], "C0", "s",
         "DES: BGS $\\times$ redMaGiC crayons (fit 20–450 kpc)")):
    z = [R[k]["zmed"] for k in keys]
    A = [R[k]["A"] for k in keys]
    s = [R[k]["sA"] for k in keys]
    des = keys[0].startswith("des")
    ax.errorbar(z, A, yerr=s, fmt=marker + ("--" if des else "-"),
                color=color, ms=6, capsize=3, lw=1.2, label=label,
                mfc="w" if des else color)
ax.axhline(1, color="k", ls="--", lw=0.9,
           label="Menard+10 amplitude ($A=1$)")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlabel("foreground median redshift", fontsize=12)
ax.set_ylabel("$A$  [$\\times$ Menard+10 template]", fontsize=12)
ax.set_title("Template amplitude versus foreground redshift "
             "(equal-count $z$ bins: 3 per Legacy sample, 2 for DES)",
             fontsize=10.5)
ax.legend(fontsize=9.5, frameon=False, loc="upper left")
ax.tick_params(labelsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/a_of_z.pdf")
fig.savefig(f"{OUT}/a_of_z.png", dpi=150)
print(f"saved -> {OUT}/a_of_z.(pdf,png)")

# ---- A(M): three luminosity bins ----------------------------------------
fig, ax = plt.subplots(figsize=(8.0, 4.8))
Mmid = [-23.0, -21.8, -20.6]      # representative bin centers (M_W1)
keys = ["bgsm0", "bgsm1", "bgsm2"]
A = [R[k]["A"] for k in keys]
s = [R[k]["sA"] for k in keys]
ax.errorbar(Mmid, A, yerr=s, fmt="o-", color="C3", ms=6, capsize=3,
            lw=1.2, label="BGS $\\times$ QSOs (fit 20–450 kpc)")
for x, k in zip(Mmid, keys):
    ax.annotate(f"$\\tilde z={R[k]['zmed']:.2f}$",
                xy=(x, R[k]["A"]), xytext=(x + 0.08, R[k]["A"] + 0.25),
                fontsize=9)
Mdes = [-22.7, -20.9]             # halves split at M_W1 = -21.889
kdes = ["desm0", "desm1"]
ax.errorbar(Mdes, [R[k]["A"] for k in kdes], yerr=[R[k]["sA"] for k in kdes],
            fmt="s--", color="C0", mfc="w", ms=6, capsize=3, lw=1.2,
            label="DES: BGS $\\times$ redMaGiC crayons (2 bins)")
for x, k in zip(Mdes, kdes):
    ax.annotate(f"$\\tilde z={R[k]['zmed']:.2f}$",
                xy=(x, R[k]["A"]), xytext=(x + 0.08, R[k]["A"] - 0.42),
                fontsize=9, color="C0")
ax.axhline(1, color="k", ls="--", lw=0.9,
           label="Menard+10 amplitude ($A=1$)")
ax.axhline(0, color="0.6", lw=0.8)
ax.invert_xaxis()                  # luminous (massive) on the right
ax.set_xlabel("absolute $W1$ magnitude (stellar-mass proxy; "
              "massive $\\rightarrow$ right)", fontsize=12)
ax.set_ylabel("$A$  [$\\times$ Menard+10 template]", fontsize=12)
ax.set_title("Template amplitude versus foreground luminosity "
             "(equal-count $M_{W1}$ bins, $0.05<z<0.45$)",
             fontsize=10.5)
ax.legend(fontsize=9.5, frameon=False, loc="upper left")
ax.tick_params(labelsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/a_of_m.pdf")
fig.savefig(f"{OUT}/a_of_m.png", dpi=150)
print(f"saved -> {OUT}/a_of_m.(pdf,png)")
