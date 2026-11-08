#!/usr/bin/env python3
"""
14_forecast.py — igmdust: can Legacy beat Menard+10 and McCleary+26?

Empirical noise model built from the MEASURED BGSxQSO pair counts
(no assumed photometric errors -- the 2026-11-08 power budget showed
the per-pair scatter is intrinsic QSO diversity):

    sigma_c(b)^2 = s_c^2 / N_pair(b)  +  sigma_sat,c(b)^2

with s_c the per-pair scatter from the shot-limited inner bins and
sigma_sat the measured excess (quasar-multiplicity saturation plus
correlated systematics; scales as 1/sqrt(N_QSO)).  Scenarios:

  bgs      the measured profile (sanity: model == data)
  replica  photgal foregrounds: N_pair x boost, same QSOs
  m10      Menard+10: their N_pair (24M gal, 85k QSO, 3,800 deg^2,
           z=0.36), our per-pair scatter, shot noise ONLY (no
           systematics floor charged to them -- flatters M10, so the
           inner-radius comparison is conservative for us)

Compares per-bin S/N against the Menard law E = EBV100 (rp/100kpc)^s
and Fisher precision of the power-law amplitude against McCleary's
+-20% (C = 3.1-4.6e-3).  Run anywhere with the npz:
    python3 14_forecast.py --npz paircount_bgs_mq_qso.npz
"""
import argparse
import os

import numpy as np

EBV100, SLOPE = 1.2e-3, -0.86                 # 01_feasibility signal model
R_COEF = {"g": 3.214, "r": 2.165, "z": 1.211, "W1": 0.184}
M10 = dict(ngal=24e6, nqso=85e3, area=3800.0, z=0.36, da=1065.0)


def law(rp_mpc):
    return EBV100 * (rp_mpc / 0.1) ** SLOPE


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--npz", required=True,
                    help="measured BGSxQSO pair counts (merged QSOs)")
    ap.add_argument("--boost", type=float, default=17.0,
                    help="replica pair-count boost (photgal/BGS density; "
                         "refine once the sweep scan reports N_gal)")
    ap.add_argument("--nqso", type=float, default=2.37e6)
    ap.add_argument("--des-npz", default=None,
                    help="measured DES crayon npz (rmgmofcsfdml): tail "
                         "amplitude precision, scaled to Legacy area")
    a = ap.parse_args()

    d = np.load(a.npz, allow_pickle=True)
    rp = np.asarray(d["rp_mid"], float)
    sig = np.asarray(d["sigma"], float)           # (ncol, nbin)
    npair = np.asarray(d["npair_data"], float)
    cnames = [str(c) for c in d["colors"]]
    edges = np.asarray(d["rp_edges"], float)
    ncol, nbin = sig.shape

    # per-pair scatter from the 4 innermost shot-limited bins
    s = np.nanmedian(sig[:, :4] * np.sqrt(npair[None, :4]), axis=1)
    shot = s[:, None] / np.sqrt(npair[None, :])
    sat2 = np.clip(sig ** 2 - shot ** 2, 0.0, None)
    print("per-pair scatter s_c:",
          {c: round(v, 3) for c, v in zip(cnames, s)})

    # M10 pair counts: annulus area (Mpc^2) -> deg^2 at their z
    mpc2deg = 1.0 / (M10["da"] * np.pi / 180.0)
    ann = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)          # Mpc^2
    ann_deg = ann * mpc2deg ** 2
    np_m10 = M10["ngal"] * (M10["nqso"] / M10["area"]) * ann_deg

    def sigma_ebv(sig_c):
        dR = np.array([R_COEF[c.split("-")[0]] - R_COEF[c.split("-")[1]]
                       for c in cnames])[:, None]
        with np.errstate(invalid="ignore", divide="ignore"):
            den = np.nansum(dR ** 2 / sig_c ** 2, 0)
            return 1.0 / np.sqrt(den)

    scen = {
        "bgs (measured)": sigma_ebv(np.sqrt(shot ** 2 + sat2)),
        f"replica (x{a.boost:.0f} pairs)":
            sigma_ebv(np.sqrt(shot ** 2 / a.boost + sat2)),
        "M10 (shot only)": sigma_ebv(s[:, None] / np.sqrt(np_m10[None, :])),
    }

    print(f"\n{'rp[kpc]':>9} {'law E':>9}"
          + "".join(f"{k:>22}" for k in scen) + "   (sigma_E | S/N)")
    totals = {k: 0.0 for k in scen}
    for b in range(nbin):
        row = f"{rp[b]*1e3:9.0f} {law(rp[b]):9.2e}"
        for k, v in scen.items():
            snr = law(rp[b]) / v[b] if np.isfinite(v[b]) else 0.0
            totals[k] += snr ** 2
            row += f"  {v[b]:9.2e}|{snr:7.1f}"
        print(row)
    print("\nprofile S/N (quadrature over bins, Menard-law signal):")
    for k, t in totals.items():
        print(f"  {k:>24}: {np.sqrt(t):6.1f} sigma")

    # Fisher precision of a power-law amplitude over the tail (>1.5 Mpc)
    m = rp > 1.5
    print("\ntail (>1.5 Mpc) amplitude precision dC/C "
          "(McCleary+26: C known to ~20%):")
    for k, v in scen.items():
        w = np.where(np.isfinite(v[m]), (law(rp[m]) / v[m]) ** 2, 0.0)
        print(f"  {k:>24}: {100.0/np.sqrt(w.sum()):5.1f} %")

    if a.des_npz:                       # crayon route: measured, not modeled
        dd = np.load(a.des_npz, allow_pickle=True)
        rpd = np.asarray(dd["rp_mid"], float)
        ebv, err = np.asarray(dd["ebv"]), np.asarray(dd["ebv_err"])
        law_d = 4.6e-3 * (0.91 * rpd * 1e3 / 100.0) ** -0.77 / 3.1
        t = rpd > 1.5
        wd = (law_d[t] / err[t]) ** 2
        dcc = 100.0 / np.sqrt(np.nansum(wd))
        print(f"\ncrayon route (McCleary law template, tail >1.5 Mpc):"
              f"\n  {'DES crayons (measured)':>24}: {dcc:5.1f} %"
              f"\n  {'Legacy crayons (x2.8 area)':>26}: "
              f"{dcc/np.sqrt(14000.0/5000.0):5.1f} %"
              "   (spec-z: no f_abs systematic)")


if __name__ == "__main__":
    main()
