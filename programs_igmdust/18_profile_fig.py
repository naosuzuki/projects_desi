import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# NOTE 2026-11-11: scratchpad was wiped overnight; the npz inputs must be
# re-fetched from astro4:~/data/igmdust/paircount/ (scp via cm-astro4gw)
# into $SC and $SC/fixwave before running.
import sys
SC = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
os.makedirs(OUT, exist_ok=True)
F_ABS = 0.806                      # speczcalib absorption factor (crayons)
RV = 3.1                           # A_V = R_V * E(B-V)
YMIN, YMAX = 0.01, 100.0           # A_V [1e-3 mag]; McCleary-wide ranges

# fix wave 2026-11-08: CSFD + 15 Mpc anchor + covariance-ML, 45 Mpc
bgs = np.load(f"{SC}/fixwave/paircount_bgs_mqcsfdml_qso.npz")
lrg = np.load(f"{SC}/fixwave/paircount_lrg_mqcsfdml_qso.npz")
# Menard replica 2026-11-10: 135.9M r<21 phot. galaxies x QSOs, 4 shards
pg = np.load(f"{SC}/fixwave/paircount_photgal_mqcsfdml_qso.npz")
mqd = np.load(f"{SC}/paircount_menarddes_qso.npz")  # redMaGiC x QSO colors
cr = np.load(f"{SC}/paircount_descrayon_qso.npz")
mc = np.load(f"{SC}/paircount_rmgmofcsfdml_qso.npz")   # reproduction (PASS)
F_ABS_MOF = 0.803                  # MOF-recalibrated absorption factor
# redMaPPer EXTENDED run (23 Mpc) from igmdust_rme_63.log
rm = dict(
    rp_mid=np.array([26, 47, 82, 143, 251, 441, 775, 1360, 2386, 4188,
                     7351, 12902, 22645]) / 1e3,
    ebv=np.array([-8.54e-2, -2.61e-2, -1.81e-2, -9.61e-4, +1.57e-4,
                  -1.52e-3, +1.25e-3, +8.10e-5, +8.25e-5, +6.06e-5,
                  -2.29e-5, +3.75e-5, -2.86e-5]),
    ebv_err=np.array([1.0e-2, 7.0e-3, 2.6e-3, 1.8e-3, 1.0e-3, 6.3e-4,
                      3.1e-4, 1.9e-4, 1.0e-4, 7.1e-5, 4.7e-5, 3.3e-5,
                      2.4e-5]))

figL, axL = plt.subplots(figsize=(10.0, 5.6))
figR, axR = plt.subplots(figsize=(10.0, 5.6))


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


def menard(ax, rmax=2.0):
    rr = np.logspace(np.log10(0.02), np.log10(rmax), 100)
    ax.plot(rr * 1e3, RV * 3.7 * (rr / 0.027) ** -0.86, "k--", lw=1,
            label="Menard+10 extrapolation")


# ---- left: DESI spectroscopy x DESI Legacy Surveys imaging (QSO colors) --
menard(axL, rmax=30.0)
plot_profile(axL, bgs["rp_mid"] * 1e3, RV * bgs["ebv"] * 1e3,
             RV * bgs["ebv_err"] * 1e3, "C3", "o",
             "BGS $\\times$ QSO colors (CSFD, cov-ML)")
plot_profile(axL, pg["rp_mid"] * 1e3, RV * pg["ebv"] * 1e3,
             RV * pg["ebv_err"] * 1e3, "C0", "s",
             "Phot. galaxies $\\times$ QSO colors (Ménard replica)",
             xoff=0.93)
g = lrg["rp_mid"] * 1e3 >= 500      # fiducial LRG mask (contamination)
plot_profile(axL, lrg["rp_mid"][g] * 1e3, RV * lrg["ebv"][g] * 1e3,
             RV * lrg["ebv_err"][g] * 1e3, "C1", "D",
             "LRG $\\times$ QSO colors ($r_p>0.5$ Mpc, CSFD, cov-ML)",
             xoff=1.07)
axL.set_title("DESI $\\times$ Legacy Surveys imaging "
              "(15 Mpc anchor, 45 Mpc reach)", fontsize=11)

# ---- right: DES imaging (redMaGiC standard crayons) ----------------------
menard(axR)
rr = np.logspace(np.log10(0.05), np.log10(30.0), 120)
axR.plot(rr * 1e3, 4.6 * (0.91 * rr / 0.1) ** -0.77, "k:", lw=1.2,
         label="McCleary+26 LRG power law")
g = cr["rp_mid"] * 1e3 >= 60
plot_profile(axR, cr["rp_mid"][g] * 1e3, RV * cr["ebv"][g] * 1e3 / F_ABS,
             RV * cr["ebv_err"][g] * 1e3 / F_ABS, "C0", "s",
             "DESI BGS $\\times$ crayons", xoff=1.07)
g = rm["rp_mid"] * 1e3 >= 150
plot_profile(axR, rm["rp_mid"][g] * 1e3, RV * rm["ebv"][g] * 1e3 / F_ABS,
             RV * rm["ebv_err"][g] * 1e3 / F_ABS, "C2", "D",
             "redMaPPer $\\lambda>20$ $\\times$ crayons", xoff=0.93)
g = mc["rp_mid"] * 1e3 >= 120
plot_profile(axR, mc["rp_mid"][g] * 1e3,
             RV * mc["ebv"][g] * 1e3 / F_ABS_MOF,
             RV * mc["ebv_err"][g] * 1e3 / F_ABS_MOF, "C4", "^",
             "redMaGiC $\\times$ crayons (MOF+CSFD, cov-ML)")
g = mqd["rp_mid"] * 1e3 >= 550      # QSO extension (f = 1)
plot_profile(axR, mqd["rp_mid"][g] * 1e3, RV * mqd["ebv"][g] * 1e3,
             RV * mqd["ebv_err"][g] * 1e3, "C1", "v",
             "redMaGiC $\\times$ merged QSOs", xoff=1.14)
axR.set_title("DES imaging (standard crayons, $/f_{\\rm abs}$)",
              fontsize=11)

for ax in (axL, axR):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(8, 50000)
    ax.set_ylim(YMIN, YMAX)
axL.set_xlim(15, 50000)     # wide rerun reaches 45 Mpc
axL.set_ylim(0.01, 100)
for ax in (axL, axR):
    ax.set_xlabel("$r_p$  [kpc, physical]", fontsize=12)
    ax.set_ylabel("$A_V$  [$10^{-3}$ mag]", fontsize=12)
    ax.legend(fontsize=9.5, loc="lower left", frameon=False)
    ax.tick_params(labelsize=11)
# ---- foreground n(z) insets (17_fgz_stats.py, run 2026-11-12) ------------
# DES-vs-Legacy amplitude comparisons carry a redshift-range difference:
# Legacy LRG <z>=0.74 vs DES redMaGiC 0.34 / redMaPPer 0.45.
import json
with open(f"{SC}/fgz_stats.json") as fh:
    FGZ = json.load(fh)


def nz_inset(ax, entries, zeff=None, zmax=1.15):
    ins = ax.inset_axes([0.615, 0.615, 0.365, 0.345])
    for key, color, label in entries:
        d = FGZ[key]
        e = np.asarray(d["edges"], float)
        h = np.asarray(d["hist"], float)
        zc = 0.5 * (e[:-1] + e[1:])
        ins.fill_between(zc, h, color=color, alpha=0.16, lw=0)
        ins.plot(zc, h, color=color, lw=1.3,
                 label=f"{label} $\\langle z\\rangle{{=}}{d['mean']:.2f}$, "
                       f"$\\tilde z{{=}}{d['median']:.2f}$")
    if zeff is not None:
        ins.axvline(zeff[0], color=zeff[1], lw=1.3, ls="--",
                    label=f"{zeff[2]} $z_{{\\rm eff}}{{=}}{zeff[0]:.2f}$")
    ins.set_xlim(0, zmax)
    ins.set_ylim(bottom=0)
    ins.set_yticks([])
    ins.set_xlabel("$z_{\\rm fg}$", fontsize=8, labelpad=1)
    ins.tick_params(labelsize=7)
    ins.legend(fontsize=6.4, frameon=False, loc="upper right",
               handlelength=1.3, borderaxespad=0.1, labelspacing=0.35)
    ins.set_title("foreground $n(z)$ (weighted)", fontsize=7.5, pad=2)


nz_inset(axL, [("bgs", "C3", "BGS"), ("lrg", "C1", "LRG")],
         zeff=(0.30, "C0", "phot. gal."), zmax=1.15)
nz_inset(axR, [("bgs", "C0", "DESI BGS"), ("redmagic", "C4", "redMaGiC"),
               ("redmapper", "C2", "redMaPPer")], zmax=0.75)

for f, name in ((figL, "av_profile_legacy"), (figR, "av_profile_des")):
    f.tight_layout()
    f.savefig(f"{OUT}/{name}.pdf")
    f.savefig(f"{OUT}/{name}.png", dpi=150)
    print("saved", f"{OUT}/{name}.pdf")
