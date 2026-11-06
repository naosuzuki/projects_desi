#!/usr/bin/env python3
"""
01_feasibility.py — igmdust phase 0: S/N forecast for the galaxy-dust reddening
measurement (Menard+10 redone with DESI DR1), per estimator and foreground split.

Forecasts the per-bin and cumulative S/N of <E(B-V)>(r_p) for:
  1a. QSO colors around BGS   (Menard+10's estimator; DESI spectroscopic QSOs)
  1a'. QSO colors around LRG
  1b. LRG "standard crayons" around BGS (Peek+15 / McCleary+25 style, spec-z)

Signal model: E(B-V)(r_p) = EBV100 * (r_p/100 kpc)^-0.86   (Menard+10 / McCleary+25)
Noise: per-source color scatter sigma_c per color, sqrt(N_pairs) averaging,
with a systematics floor added in quadrature per bin.

Catalog counts are the measured DR1 values (FITS headers, 2026-11-06).
Areas are provisional until computed from the randoms in 02.

Output: table to stdout + figure $HOME/data/igmdust/feasibility/snr_forecast.png
Pure forecast - needs numpy/matplotlib/astropy only, no catalogs.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

# ----------------------------------------------------------------- config ---
# measured DR1 row counts (00_download + FITS NAXIS2, 2026-11-06)
N_BGS = 4_081_227 + 1_441_126        # BGS_ANY clustering NGC+SGC
N_LRG = 1_476_135 + 662_492          # LRG clustering NGC+SGC
N_QSO_CAT = 2_182_309                # QSO_cat_iron_cumulative_v0 (all QSOs)

# provisional footprint areas [deg^2] — replace with randoms-based areas in 02
AREA_BGS = 7_400.0
AREA_DARK = 9_100.0

# usable background QSOs: main survey, z_q > z_fg + 0.2, BAL cut (~65% of cat)
F_QSO_USABLE = 0.65
Z_BGS, Z_LRG = 0.20, 0.65            # median foreground redshifts
QSO_ZMED = 1.7

# signal normalization: E(B-V) at 100 kpc [mag] (M10 ~1e-3; McCleary+25 A_V/3.1)
EBV100, SLOPE = 1.2e-3, -0.86

# extinction coefficients R_X for DECam/WISE bands (SFD-style, R_V=3.1 MW law;
# SMC at z_fg modifies these at the ~20% level — fine for a forecast)
R = {"g": 3.214, "r": 2.165, "z": 1.211, "W1": 0.184}
COLORS = [("g", "r"), ("r", "z"), ("g", "z"), ("z", "W1")]

# per-source color scatter [mag] (intrinsic + photometric, per color)
SIGMA_C = {"QSO": 0.15, "LRG": 0.05}
SYS_FLOOR = 1.0e-4                   # per-bin systematics floor on E(B-V) [mag]

# binning: r_p = 20 kpc - 10 Mpc, 4 bins/dex (physical at z_fg)
RP_EDGES = np.logspace(np.log10(0.02), np.log10(10.0), int(4 * 2.7) + 1)  # Mpc

# ------------------------------------------------------------- machinery ---
def forecast(name, n_fg, area_fg, n_src, area_src, z_fg, sigma_c):
    """Per-bin S/N on E(B-V) for one foreground x background combination."""
    sigma_src = n_src / area_src                       # sources per deg^2
    da = cosmo.angular_diameter_distance(z_fg).to(u.Mpc).value
    deg_per_mpc = np.degrees(1.0 / da)
    lo, hi = RP_EDGES[:-1], RP_EDGES[1:]
    rp_mid = np.sqrt(lo * hi)
    ann_deg2 = np.pi * (hi**2 - lo**2) * deg_per_mpc**2
    n_pairs = n_fg * sigma_src * ann_deg2
    ebv = EBV100 * (rp_mid / 0.1) ** SLOPE
    # E(B-V) sensitivity from all colors: sigma_EBV^-2 = N_pairs * sum (dR_i/sigma_c)^2
    w = sum(((R[a] - R[b]) / sigma_c) ** 2 for a, b in COLORS)
    sigma_ebv = np.sqrt(1.0 / (n_pairs * w) + SYS_FLOOR**2)
    return dict(name=name, rp=rp_mid, ebv=ebv, sigma=sigma_ebv,
                snr=ebv / sigma_ebv, n_pairs=n_pairs)

cases = [
    forecast("QSO x BGS  (Menard est.)", N_BGS, AREA_BGS,
             F_QSO_USABLE * N_QSO_CAT, AREA_DARK, Z_BGS, SIGMA_C["QSO"]),
    forecast("QSO x LRG", N_LRG, AREA_DARK,
             F_QSO_USABLE * N_QSO_CAT, AREA_DARK, Z_LRG, SIGMA_C["QSO"]),
    forecast("LRG crayons x BGS (McCleary est.)", N_BGS, AREA_BGS,
             N_LRG, AREA_DARK, Z_BGS, SIGMA_C["LRG"]),
]

# ---------------------------------------------------------------- report ---
print(f"igmdust S/N forecast  |  EBV(rp)={EBV100:.1e}(rp/100kpc)^{SLOPE}"
      f"  |  sys floor {SYS_FLOOR:.0e} mag/bin")
print(f"BGS={N_BGS:,}  LRG={N_LRG:,}  QSO_usable={F_QSO_USABLE*N_QSO_CAT:,.0f}"
      f"  (areas {AREA_BGS:.0f}/{AREA_DARK:.0f} deg2, provisional)\n")
for c in cases:
    inner = slice(0, np.searchsorted(RP_EDGES, 1.0))   # 20 kpc - 1 Mpc
    tot = np.sqrt((c["snr"][inner] ** 2).sum())
    print(f"--- {c['name']}:  cumulative S/N(20 kpc-1 Mpc) = {tot:.0f}")
    for rp, ebv, sg, sn, npair in zip(c["rp"], c["ebv"], c["sigma"],
                                      c["snr"], c["n_pairs"]):
        print(f"  rp={rp*1000:8.0f} kpc  pairs={npair:12.3e}"
              f"  EBV={ebv:9.2e}  sig={sg:9.2e}  S/N={sn:7.1f}")
    print()

# ---------------------------------------------------------------- figure ---
out = os.path.expanduser("~/data/igmdust/feasibility")
os.makedirs(out, exist_ok=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
for c in cases:
    ax1.loglog(c["rp"] * 1e3, c["snr"], "o-", label=c["name"])
    ax2.loglog(c["rp"] * 1e3, c["sigma"], "o--", label=c["name"])
ax1.axhline(1, color="k", lw=0.5)
ax1.set(xlabel=r"$r_p$ [kpc]", ylabel="S/N per bin", title="per-bin detection S/N")
rp = cases[0]["rp"]
ax2.loglog(rp * 1e3, cases[0]["ebv"], "k-", lw=2, label="signal E(B-V)")
ax2.set(xlabel=r"$r_p$ [kpc]", ylabel=r"E(B$-$V) [mag]",
        title="signal vs per-bin error")
ax1.legend(fontsize=8); ax2.legend(fontsize=8)
fig.suptitle("igmdust DR1 pilot forecast (4 bins/dex, DR1 counts, provisional areas)")
fig.savefig(f"{out}/snr_forecast.png", dpi=150)
print(f"figure -> {out}/snr_forecast.png")
