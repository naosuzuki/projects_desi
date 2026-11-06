#!/usr/bin/env python3
"""
01_feasibility.py — igmdust phase 0: S/N forecast for the galaxy-dust reddening
measurement (Menard+10 redone with DESI), per estimator, foreground split,
and imaging survey (Legacy Surveys / DES / HSC).

Estimator combinations forecast:
  LS pilot (DESI DR1 spectroscopy x Legacy Surveys photometry):
    1a.  QSO colors around BGS   (Menard+10's estimator)
    1a'. QSO colors around LRG
    1b.  LRG "standard crayons" around BGS (Peek+15 / McCleary+25, spec-z)
  DES  (DESI BGS foregrounds in the DESI x DES overlap; DES red-sequence
        "crayon" galaxies as backgrounds, redMaGiC-like density)
  HSC  (DESI BGS foregrounds in the DESI x HSC-Wide overlap; HSC quiescent
        photometric galaxies to i~24 as backgrounds -- depth axis)

Signal model: E(B-V)(r_p) = EBV100 * (r_p/100 kpc)^-0.86   (Menard+10 / McCleary+25)
Noise: per-source color scatter per color, sqrt(N_pairs) averaging, and a
per-bin systematics floor in quadrature. E(B-V) sensitivity combines the
survey's colors with extinction-law weights:
   sigma_EBV^-2 = N_pairs * sum_colors [(R_a - R_b)/sigma_c]^2.

DESI DR1 counts are measured (FITS headers, 2026-11-06). Areas, overlap
areas, and the DES/HSC background densities and scatters are assumptions,
flagged in the config below; replace with measured values as they land.

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
AREA_BGS = 7_400.0                   # DESI DR1 BGS footprint
AREA_DARK = 9_100.0                  # DESI DR1 dark-time footprint
AREA_DESI_X_DES = 1_500.0            # [ASSUMED] DESI spec coverage inside DES
AREA_DESI_X_HSC = 1_000.0            # [ASSUMED] DESI spec coverage inside HSC-Wide

# usable background QSOs: main survey, z_q > z_fg + 0.2, BAL cut (~65% of cat)
F_QSO_USABLE = 0.65
Z_BGS, Z_LRG = 0.20, 0.65            # median foreground redshifts

# [ASSUMED] background "crayon" samples in the imaging surveys
DES_CRAYON_DENS = 650.0              # red-sequence gal / deg^2 (redMaGiC-like)
DES_CRAYON_SIG = 0.04                # per-color scatter [mag]
HSC_CRAYON_DENS = 10_800.0           # quiescent gal to i~24 / deg^2 (~3/arcmin^2)
HSC_CRAYON_SIG = 0.06                # per-color scatter [mag] (fainter, photo-z)

# signal normalization: E(B-V) at 100 kpc [mag] (M10 ~1e-3; McCleary+25 A_V/3.1)
EBV100, SLOPE = 1.2e-3, -0.86

# extinction coefficients R_X (SFD-style, R_V=3.1 MW law; SMC at z_fg modifies
# these at the ~20% level — fine for a forecast)
R = {"g": 3.214, "r": 2.165, "i": 1.592, "z": 1.211, "W1": 0.184}
COLORS_LS = [("g", "r"), ("r", "z"), ("g", "z"), ("z", "W1")]   # DECam+WISE
COLORS_GRIZ = [("g", "r"), ("r", "i"), ("i", "z"), ("g", "z")]  # DES / HSC

# per-source color scatter [mag] for the LS pilot sources
SIGMA_QSO, SIGMA_LRG = 0.15, 0.05
SYS_FLOOR = 1.0e-4                   # per-bin systematics floor on E(B-V) [mag]

# binning: r_p = 20 kpc - 10 Mpc, 4 bins/dex (physical at z_fg)
RP_EDGES = np.logspace(np.log10(0.02), np.log10(10.0), int(4 * 2.7) + 1)  # Mpc

# ------------------------------------------------------------- machinery ---
def forecast(name, n_fg, n_src, area_src, z_fg, sigma_c, colors):
    """Per-bin S/N on E(B-V) for one foreground x background combination."""
    sigma_src = n_src / area_src                       # sources per deg^2
    da = cosmo.angular_diameter_distance(z_fg).to(u.Mpc).value
    deg_per_mpc = np.degrees(1.0 / da)
    lo, hi = RP_EDGES[:-1], RP_EDGES[1:]
    rp_mid = np.sqrt(lo * hi)
    ann_deg2 = np.pi * (hi**2 - lo**2) * deg_per_mpc**2
    n_pairs = n_fg * sigma_src * ann_deg2
    ebv = EBV100 * (rp_mid / 0.1) ** SLOPE
    w = sum(((R[a] - R[b]) / sigma_c) ** 2 for a, b in colors)
    sigma_ebv = np.sqrt(1.0 / (n_pairs * w) + SYS_FLOOR**2)
    return dict(name=name, rp=rp_mid, ebv=ebv, sigma=sigma_ebv,
                snr=ebv / sigma_ebv, n_pairs=n_pairs)

bgs_dens = N_BGS / AREA_BGS                            # BGS per deg^2

cases = [
    # --- Legacy Surveys pilot (DESI DR1 spectroscopy, LS photometry) ---
    forecast("LS: QSO x BGS  (Menard est.)", N_BGS,
             F_QSO_USABLE * N_QSO_CAT, AREA_DARK, Z_BGS, SIGMA_QSO, COLORS_LS),
    forecast("LS: QSO x LRG", N_LRG,
             F_QSO_USABLE * N_QSO_CAT, AREA_DARK, Z_LRG, SIGMA_QSO, COLORS_LS),
    forecast("LS: LRG crayons x BGS (McCleary est.)", N_BGS,
             N_LRG, AREA_DARK, Z_BGS, SIGMA_LRG, COLORS_LS),
    # --- DES: DESI BGS foregrounds in overlap, DES red-seq crayons behind ---
    forecast("DES: crayons x BGS (overlap 1500 deg2)",
             bgs_dens * AREA_DESI_X_DES,
             DES_CRAYON_DENS * AREA_DESI_X_DES, AREA_DESI_X_DES,
             Z_BGS, DES_CRAYON_SIG, COLORS_GRIZ),
    # --- HSC: DESI BGS foregrounds in overlap, deep quiescent gals behind ---
    forecast("HSC: crayons x BGS (overlap 1000 deg2)",
             bgs_dens * AREA_DESI_X_HSC,
             HSC_CRAYON_DENS * AREA_DESI_X_HSC, AREA_DESI_X_HSC,
             Z_BGS, HSC_CRAYON_SIG, COLORS_GRIZ),
]

# ---------------------------------------------------------------- report ---
print(f"igmdust S/N forecast  |  EBV(rp)={EBV100:.1e}(rp/100kpc)^{SLOPE}"
      f"  |  sys floor {SYS_FLOOR:.0e} mag/bin")
print(f"BGS={N_BGS:,}  LRG={N_LRG:,}  QSO_usable={F_QSO_USABLE*N_QSO_CAT:,.0f}"
      f"  (areas {AREA_BGS:.0f}/{AREA_DARK:.0f} deg2 provisional;"
      f" DESIxDES={AREA_DESI_X_DES:.0f}, DESIxHSC={AREA_DESI_X_HSC:.0f} assumed)\n")
inner = slice(0, int(np.searchsorted(RP_EDGES, 1.0)))   # 20 kpc - ~1 Mpc
for c in cases:
    tot = np.sqrt((c["snr"][inner] ** 2).sum())
    cum_in = c["n_pairs"][inner].sum()
    print(f"--- {c['name']}:  S/N(<~1 Mpc) = {tot:.0f}"
          f"   pairs(<~1.5 Mpc) = {cum_in:.2e}   pairs(<10 Mpc) = {c['n_pairs'].sum():.2e}")
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
    ax1.loglog(c["rp"] * 1e3, c["snr"], "o-", ms=3.5, label=c["name"])
    ax2.loglog(c["rp"] * 1e3, c["sigma"], "o--", ms=3.5, label=c["name"])
ax1.axhline(1, color="k", lw=0.5)
ax1.set(xlabel=r"$r_p$ [kpc]", ylabel="S/N per bin", title="per-bin detection S/N")
rp = cases[0]["rp"]
ax2.loglog(rp * 1e3, cases[0]["ebv"], "k-", lw=2, label="signal E(B-V)")
ax2.set(xlabel=r"$r_p$ [kpc]", ylabel=r"E(B$-$V) [mag]",
        title="signal vs per-bin error")
ax1.legend(fontsize=7); ax2.legend(fontsize=7)
fig.suptitle("igmdust forecast: LS pilot + DES/HSC extensions "
             "(4 bins/dex; DR1 counts; overlap areas assumed)")
fig.savefig(f"{out}/snr_forecast.png", dpi=150)
print(f"figure -> {out}/snr_forecast.png")
