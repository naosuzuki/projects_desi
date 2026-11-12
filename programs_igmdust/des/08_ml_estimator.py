#!/usr/bin/env python3
"""
08_ml_estimator.py — igmdust/DES: covariance-weighted maximum-likelihood
E(B-V), the final rung of the McCleary-reproduction ladder.

WHY: McCleary (2026) estimates reddening with a Gaussian ML estimator whose
covariance is "calculated directly from the data vector", encapsulating
measurement noise AND cosmic variance, including color-color correlations.
Our baseline solve weights each color independently (w_c = 1/sigma_c^2),
which gives the g-band colors (largest dust response) the largest weight.
Iteration 3 showed the MOF profile inflated at ALL scales in g-r/g-i while
i-z stayed put — correlated g systematics that full-covariance weighting
suppresses by construction.

HOW: the paircount npz (02_paircount.py with jk_* arrays saved) carries the
per-(color, region, bin) accumulations. We rebuild the leave-one-region-out
delta vectors d_k, form the jackknife color-color covariance per radial bin,

    C_i = (K-1)/K  sum_k (d_k - mean)(d_k - mean)^T    (Hartlap-corrected)

and solve
    E_i     = dR^T C^-1 delta / dR^T C^-1 dR
    sigma_i = (dR^T C^-1 dR)^-1/2 .

Writes paircount_<name>ml_qso.npz (same schema), assessable by
07_assess_mccleary.py --npz.

Usage:  python3 08_ml_estimator.py --npz .../paircount_rmgmofcsfd_qso.npz
"""
import argparse
import json
import os
import time

import numpy as np

R_DES = {"g": 3.186, "r": 2.140, "i": 1.569, "z": 1.196}


def loo_deltas(d):
    """Leave-one-region-out delta vectors, (K, ncol, nbin)."""
    sw_d, swc_d = d["jk_sw_d"], d["jk_swc_d"]
    sw_r, swc_r = d["jk_sw_r"], d["jk_swc_r"]
    K = sw_d.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        dk = np.array([
            (swc_d.sum(1) - swc_d[:, k]) / (sw_d.sum(1) - sw_d[:, k])
            - (swc_r.sum(1) - swc_r[:, k]) / (sw_r.sum(1) - sw_r[:, k])
            for k in range(K)])
    return dk


def ml_solve(d):
    cols = [c.split("-") for c in d["colors"]]
    dR = np.array([R_DES[a] - R_DES[b] for a, b in cols])
    delta = d["delta"]                                  # (ncol, nbin)
    dk = loo_deltas(d)                                  # (K, ncol, nbin)
    K, ncol, nbin = dk.shape
    hartlap = (K - ncol - 2) / (K - 1)
    E = np.full(nbin, np.nan)
    Eerr = np.full(nbin, np.nan)
    for i in range(nbin):
        x = dk[:, :, i]
        good = np.all(np.isfinite(x), 1)
        if good.sum() < ncol + 2 or not np.all(np.isfinite(delta[:, i])):
            continue
        x = x[good]
        k = len(x)
        m = x.mean(0)
        C = (k - 1) / k * (x - m).T @ (x - m)           # jackknife covariance
        Cinv = hartlap * np.linalg.pinv(C)
        den = dR @ Cinv @ dR
        E[i] = dR @ Cinv @ delta[:, i] / den
        Eerr[i] = 1.0 / np.sqrt(den)
    return E, Eerr


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--npz", required=True,
                    help="paircount npz with jk_* arrays")
    a = ap.parse_args()
    d = np.load(a.npz, allow_pickle=True)
    if "jk_sw_d" not in d:
        raise SystemExit(f"{a.npz} lacks jk_* arrays -- rerun the pair "
                         "count with the updated 02_paircount.py")
    E, Eerr = ml_solve(d)
    ebv0, err0 = d["ebv"], d["ebv_err"]
    print(f"== ML (full color covariance) vs diagonal solve: "
          f"{os.path.basename(a.npz)}")
    print("   rp[kpc]    E_ML +- err          E_diag +- err      ratio")
    for i, rp in enumerate(d["rp_mid"]):
        r = E[i] / ebv0[i] if abs(ebv0[i]) > 0 else np.nan
        print(f"  {rp*1e3:8.0f}  {E[i]:+9.2e} {Eerr[i]:8.1e}"
              f"   {ebv0[i]:+9.2e} {err0[i]:8.1e}   {r:6.2f}")

    meta = json.loads(str(d["meta"]))
    meta.update(estimator="ml_fullcov", date=time.strftime("%F %T"))
    out = a.npz.replace("_qso.npz", "ml_qso.npz")
    np.savez(out, rp_edges=d["rp_edges"], rp_mid=d["rp_mid"],
             delta=d["delta"], sigma=d["sigma"], ebv=E, ebv_err=Eerr,
             npair_data=d["npair_data"], npair_ran=d["npair_ran"],
             colors=d["colors"], meta=json.dumps(meta))
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
