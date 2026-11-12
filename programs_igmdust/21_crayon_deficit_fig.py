"""21_crayon_deficit_fig.py -- figure script of record for the
two-recipe exhibit (ms.tex fig:crayondeficit): LRG spec-z crayons vs
merged-quasar colors around the SAME DR1 BGS foregrounds, |E(B-V)| on
log-log axes.

v2 (default, 2026-11-12, -> crayon_deficit_legacy_v2.pdf/png), per
Nao: no +/- sign coding (every >=1 sigma bin is one filled symbol at
|E|; the sign story stays in the caption), and each panel gains a
redshift inset -- the BGS foreground n(z) (17_fgz_stats.py) with that
panel's background n(z) (19_bgz_stats.py: LRG spec-z crayons left,
merged QSOs right), annotated with mean and median redshifts.
Legends move to lower left to clear the insets.  --v1 reproduces the
original 2026-11-11 sign-coded figure byte-for-byte (recovered from
the session transcript after the scratchpad wipe; it was never a repo
file before).

Inputs in NPZDIR (argv[1], default script dir; master copies on
astro4:~/data/igmdust/paircount/): fixwave/paircount_bgslrgcrayonml_
qso.npz and fixwave/paircount_bgs_mqcsfdml_qso.npz, plus
fgz_stats_2026-11-12.json / bgz_stats_2026-11-12.json in the repo.

    python3 21_crayon_deficit_fig.py [NPZDIR]        # v2 figure
    python3 21_crayon_deficit_fig.py [NPZDIR] --v1   # original
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
                help="reproduce the original 2026-11-11 figure")
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))
SC = a.npzdir or HERE
OUT = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
YMIN, YMAX = 8e-3, 40.0            # |E(B-V)| in 1e-3 mag
RV = 3.1

cr = np.load(f"{SC}/fixwave/paircount_bgslrgcrayonml_qso.npz")
qs = np.load(f"{SC}/fixwave/paircount_bgs_mqcsfdml_qso.npz")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)


def plot_signed(ax, x, e, s, color, marker, label):
    """v1 look: |E| on log axes; filled = positive, open = negative;
    <1 sigma bins become 2 sigma upper limits on |E| (small arrows)."""
    x, e, s = np.asarray(x), np.asarray(e), np.asarray(s)
    det = np.abs(e) > s
    for sign, mfc, lab in ((e > 0, color, label + " ($+$)"),
                           (e < 0, "white", label + " ($-$, at $|$E$|$)")):
        m = det & sign
        if not m.any():
            continue
        a_ = np.abs(e[m])
        lo = np.clip(np.minimum(s[m], a_ - YMIN * 1.05), 0, None)
        ax.errorbar(x[m], a_, yerr=[lo, s[m]], fmt=marker, color=color,
                    mfc=mfc, ms=5.5, capsize=2, lw=1.2, label=lab)
    ul = np.abs(e[~det]) + 2 * s[~det]
    ok = ul > YMIN * 1.3
    ax.errorbar(x[~det][ok], ul[ok], yerr=ul[ok] * 0.22, uplims=True,
                fmt=marker, color=color, ms=3.0, mfc="none", lw=0.9,
                alpha=0.6)


def plot_abs(ax, x, e, s, color, marker, label):
    """v2 look: one filled symbol per >=1 sigma bin at |E|, no sign
    coding; <1 sigma bins stay 2 sigma upper limits (arrows)."""
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


plot = plot_signed if a.v1 else plot_abs
rr = np.logspace(np.log10(0.012), np.log10(40.0), 150)   # Mpc

# ---- left: LRG spec-z crayons, McCleary+26 law as guide -----------------
axL.plot(rr * 1e3, (4.6 / RV) * (0.91 * rr / 0.1) ** -0.77, "k:",
         lw=1.2, label="McCleary+26 LRG power law")
plot(axL, cr["rp_mid"] * 1e3, cr["ebv"] * 1e3, cr["ebv_err"] * 1e3,
     "C0", "s", "LRG spec-$z$ crayons")
axL.set_title("LRG spec-$z$ standard crayons ($f_{\\rm abs}\\equiv1$)",
              fontsize=10.5)

# ---- right: merged-quasar colors, Menard+10 law as guide ----------------
axR.plot(rr * 1e3, 3.7 * (rr / 0.027) ** -0.86, "k--", lw=1.1,
         label="Ménard+10 extrapolation")
plot(axR, qs["rp_mid"] * 1e3, qs["ebv"] * 1e3, qs["ebv_err"] * 1e3,
     "C3", "o", "merged-quasar colors")
axR.set_title("Merged-quasar colors (CSFD, cov-ML)", fontsize=10.5)

for ax in (axL, axR):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 45000)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=12)
    ax.legend(fontsize=9.5,
              loc="upper right" if a.v1 else "lower left", frameon=False)
    ax.tick_params(labelsize=11)
axL.set_ylabel("$|{\\rm E(B}-{\\rm V)}|$  [$10^{-3}$ mag]", fontsize=12)

if a.v1:
    suffix = ""
    fig.suptitle("Two recipes, one foreground: backgrounds around DR1 BGS "
                 "(filled $=$ red excess, open $=$ deficit)", fontsize=11,
                 y=0.995)
else:                                   # ---- v2: n(z) insets ----
    suffix = "_v2"
    fig.suptitle("Two recipes, one foreground: backgrounds around DR1 BGS",
                 fontsize=11, y=0.995)
    with open(f"{HERE}/fgz_stats_2026-11-12.json") as fh:
        FGZ = json.load(fh)
    with open(f"{HERE}/bgz_stats_2026-11-12.json") as fh:
        BGZ = json.load(fh)

    def nz_inset(ax, entries, zmax):
        ins = ax.inset_axes([0.615, 0.615, 0.365, 0.345])
        for d_, color, label in entries:
            e = np.asarray(d_["edges"], float)
            h = np.asarray(d_["hist"], float)
            h = h / h.max()             # peak-normalized shapes
            zc = 0.5 * (e[:-1] + e[1:])
            ins.fill_between(zc, h, color=color, alpha=0.16, lw=0)
            ins.plot(zc, h, color=color, lw=1.2,
                     label=f"{label} $\\langle z\\rangle{{=}}"
                           f"{d_['mean']:.2f}$,"
                           f" $\\tilde z{{=}}{d_['median']:.2f}$")
        ins.set_xlim(0, zmax)
        ins.set_ylim(0, 1.6)
        ins.set_yticks([])
        ins.set_xlabel("$z$", fontsize=8, labelpad=1)
        ins.tick_params(labelsize=7)
        ins.legend(fontsize=6.2, frameon=False, loc="upper right",
                   handlelength=1.1, borderaxespad=0.1, labelspacing=0.35)
        ins.set_title("fg / bg $n(z)$", fontsize=7.5, pad=2)

    nz_inset(axL, [(FGZ["bgs"], "C3", "fg: BGS"),
                   (BGZ["lrgcrayon"], "C0", "bg: LRG crayons")], zmax=1.15)
    nz_inset(axR, [(FGZ["bgs"], "C0", "fg: BGS"),
                   (BGZ["qso"], "C3", "bg: QSOs")], zmax=4.0)

fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{OUT}/crayon_deficit_legacy{suffix}.pdf")
fig.savefig(f"{OUT}/crayon_deficit_legacy{suffix}.png", dpi=150)
print("saved", f"{OUT}/crayon_deficit_legacy{suffix}.pdf")
