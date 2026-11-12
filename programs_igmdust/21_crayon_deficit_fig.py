"""21_crayon_deficit_fig.py -- figure script of record for the
two-recipe exhibit (ms.tex fig:crayondeficit / fig:crayondeficitqso):
LRG spec-z crayons vs merged-quasar colors around the SAME DR1 BGS
foregrounds, |E(B-V)| on log-log axes.

Versions (Nao, 2026-11-12: one plot per file -- never two panels in
one figure; try insets and separate n(z) as different versions):
  v2  crayon_deficit_legacy_v2_crayon.pdf/png and _v2_qso.pdf/png --
      one single-panel plot per recipe, no +/- sign coding (every
      >=1 sigma bin is one filled symbol at |E|; the sign story
      stays in the captions), each with a redshift inset: the BGS
      foreground n(z) (17_fgz_stats.py) with that plot's background
      n(z) (19_bgz_stats.py), annotated with mean and median z.
      (The 2026-11-12 two-panel combined v2 is superseded.)
  v3  crayon_deficit_legacy_v3_crayon/_v3_qso (profiles without the
      inset) plus crayon_deficit_legacy_v3_nz -- ONE standalone
      redshift-distribution figure shared by both profiles (fg BGS
      + bg LRG crayons + bg QSOs; both recipes share the same
      foreground).
  --v1 reproduces the original 2026-11-11 sign-coded two-panel
      figure byte-for-byte (recovered from the session transcript
      after the scratchpad wipe).

Inputs in NPZDIR (argv[1], default ~/data/igmdust/paircount, local
mirror of astro4:~/data/igmdust/paircount/):
fixwave/paircount_bgslrgcrayonml_qso.npz and
fixwave/paircount_bgs_mqcsfdml_qso.npz, plus
fgz_stats_2026-11-12.json / bgz_stats_2026-11-12.json in the repo.

    python3 21_crayon_deficit_fig.py [NPZDIR]        # v2 + v3
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
SC = a.npzdir or os.path.expanduser("~/data/igmdust/paircount")
OUT = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
YMIN, YMAX = 8e-3, 40.0            # |E(B-V)| in 1e-3 mag
RV = 3.1

cr = np.load(f"{SC}/fixwave/paircount_bgslrgcrayonml_qso.npz")
qs = np.load(f"{SC}/fixwave/paircount_bgs_mqcsfdml_qso.npz")


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
    """v2/v3 look: one filled symbol per >=1 sigma bin at |E|, no sign
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


rr = np.logspace(np.log10(0.012), np.log10(40.0), 150)   # Mpc


def draw_crayon(ax, plot):
    """Left recipe: LRG spec-z crayons, McCleary+26 law as guide."""
    ax.plot(rr * 1e3, (4.6 / RV) * (0.91 * rr / 0.1) ** -0.77, "k:",
            lw=1.2, label="McCleary+26 LRG power law")
    plot(ax, cr["rp_mid"] * 1e3, cr["ebv"] * 1e3, cr["ebv_err"] * 1e3,
         "C0", "s", "LRG spec-$z$ crayons")
    ax.set_title("LRG spec-$z$ standard crayons ($f_{\\rm abs}\\equiv1$)",
                 fontsize=10.5)


def draw_qso(ax, plot):
    """Right recipe: merged-quasar colors, Menard+10 law as guide."""
    ax.plot(rr * 1e3, 3.7 * (rr / 0.027) ** -0.86, "k--", lw=1.1,
            label="Ménard+10 extrapolation")
    plot(ax, qs["rp_mid"] * 1e3, qs["ebv"] * 1e3, qs["ebv_err"] * 1e3,
         "C3", "o", "merged-quasar colors")
    ax.set_title("Merged-quasar colors (CSFD, cov-ML)", fontsize=10.5)


def style(ax, legend_loc):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 45000)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=12)
    ax.legend(fontsize=9.5, loc=legend_loc, frameon=False)
    ax.tick_params(labelsize=11)


YLAB = "$|{\\rm E(B}-{\\rm V)}|$  [$10^{-3}$ mag]"

if a.v1:
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)
    draw_crayon(axL, plot_signed)
    draw_qso(axR, plot_signed)
    for ax in (axL, axR):
        style(ax, "upper right")
    axL.set_ylabel(YLAB, fontsize=12)
    fig.suptitle("Two recipes, one foreground: backgrounds around DR1 BGS "
                 "(filled $=$ red excess, open $=$ deficit)", fontsize=11,
                 y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f"{OUT}/crayon_deficit_legacy.pdf")
    fig.savefig(f"{OUT}/crayon_deficit_legacy.png", dpi=150)
    print("saved", f"{OUT}/crayon_deficit_legacy.pdf")
    raise SystemExit

# ---------------- v2 + v3: one plot per file ----------------
with open(f"{HERE}/fgz_stats_2026-11-12.json") as fh:
    FGZ = json.load(fh)
with open(f"{HERE}/bgz_stats_2026-11-12.json") as fh:
    BGZ = json.load(fh)


def nz_curve(ax, d_, color, label, lw=1.2, legend_z=True):
    e = np.asarray(d_["edges"], float)
    h = np.asarray(d_["hist"], float)
    h = h / h.max()                     # peak-normalized shapes
    zc = 0.5 * (e[:-1] + e[1:])
    ax.fill_between(zc, h, color=color, alpha=0.16, lw=0)
    lab = label
    if legend_z:
        lab += (f" $\\langle z\\rangle{{=}}{d_['mean']:.2f}$,"
                f" $\\tilde z{{=}}{d_['median']:.2f}$")
    ax.plot(zc, h, color=color, lw=lw, label=lab)


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


PANELS = (
    ("crayon", draw_crayon,
     [(FGZ["bgs"], "C3", "fg: BGS"),
      (BGZ["lrgcrayon"], "C0", "bg: LRG crayons")], 1.15),
    ("qso", draw_qso,
     [(FGZ["bgs"], "C0", "fg: BGS"),
      (BGZ["qso"], "C3", "bg: QSOs")], 4.0),
)

for key, draw, entries, zmax in PANELS:
    # ---- v2: single-panel profile with the n(z) inset ----
    fig, ax = plt.subplots(figsize=(6.2, 4.9))
    draw(ax, plot_abs)
    style(ax, "lower left")
    ax.set_ylabel(YLAB, fontsize=12)
    nz_inset(ax, entries, zmax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/crayon_deficit_legacy_v2_{key}.pdf")
    fig.savefig(f"{OUT}/crayon_deficit_legacy_v2_{key}.png", dpi=150)
    print("saved", f"{OUT}/crayon_deficit_legacy_v2_{key}.pdf")
    plt.close(fig)

    # ---- v3: the same profile without the inset ----
    fig, ax = plt.subplots(figsize=(6.2, 4.9))
    draw(ax, plot_abs)
    style(ax, "lower left")
    ax.set_ylabel(YLAB, fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT}/crayon_deficit_legacy_v3_{key}.pdf")
    fig.savefig(f"{OUT}/crayon_deficit_legacy_v3_{key}.png", dpi=150)
    print("saved", f"{OUT}/crayon_deficit_legacy_v3_{key}.pdf")
    plt.close(fig)

# ---- v3 companion: one standalone n(z) figure for both profiles ----
fig, ax = plt.subplots(figsize=(6.2, 4.0))
nz_curve(ax, FGZ["bgs"], "C2", "foreground: BGS")
nz_curve(ax, BGZ["lrgcrayon"], "C0", "background: LRG spec-$z$ crayons")
nz_curve(ax, BGZ["qso"], "C3", "background: merged QSOs")
ax.set_xlim(0, 4.0)
ax.set_ylim(0, 1.45)
ax.set_yticks([])
ax.set_xlabel("$z$", fontsize=12)
ax.set_ylabel("$n(z)$  [peak-normalized]", fontsize=12)
ax.legend(fontsize=9.5, frameon=False, loc="upper right")
ax.tick_params(labelsize=11)
ax.set_title("Redshift distributions: shared BGS foreground, two "
             "background recipes", fontsize=10.5)
fig.tight_layout()
fig.savefig(f"{OUT}/crayon_deficit_legacy_v3_nz.pdf")
fig.savefig(f"{OUT}/crayon_deficit_legacy_v3_nz.png", dpi=150)
print("saved", f"{OUT}/crayon_deficit_legacy_v3_nz.pdf")
