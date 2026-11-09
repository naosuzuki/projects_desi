#!/usr/bin/env python3
"""
10_ml_legacy.py — igmdust: covariance-weighted ML E(B-V) for Legacy runs.

The third DES ingredient (projects_des 08_ml_estimator.py), adapted to
the Legacy QSO-color npz.  Per radial bin,

    E = dR^T C^-1 delta / dR^T C^-1 dR

with C the color-color covariance of delta built from the per-region
jackknife accumulations that 02_paircount saves (Hartlap-corrected),
which downweights the spatially correlated g-band directions a diagonal
1/sigma^2 solve amplifies.  Differences from the DES script:
  - colors {g-r, r-z, z-W1} only (g-z is their sum: a singular C);
    z-W1 dropped inside WISE_RP_MIN (6" blending);
  - QSO backgrounds carry per-color constant imaging offsets, so O_c is
    estimated from the bins beyond --anchor (default 10 Mpc, where the
    Menard law is < 2e-5) and subtracted — per jackknife sample, so the
    anchor noise propagates into C;
  - --null <npz> instead subtracts a measured star-null profile
    (interpolated in log r_p) before the anchor step.

    python3 10_ml_legacy.py --npz $OUT/paircount_bgs_mqcsfd_qso.npz
writes paircount_<name>ml_qso.npz (ebv/ebv_err replaced; same schema).
"""
import argparse
import json
import os

import numpy as np

R_COEF = {"g": 3.214, "r": 2.165, "z": 1.211, "W1": 0.184}
COLS_USE = ["g-r", "r-z", "z-W1"]      # independent subset
WISE_RP_MIN = 0.15                     # Mpc; z-W1 dropped inside


def loo_deltas(d):
    """Full-sample delta and leave-one-region-out deltas (K, ncol, nbin)
    from the saved accumulations, with the exact 02 jackknife formula."""
    sw_d, swc_d = d["jk_sw_d"], d["jk_swc_d"]
    sw_r, swc_r = d["jk_sw_r"], d["jk_swc_r"]
    K = sw_d.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        delta = swc_d.sum(1) / sw_d.sum(1) - swc_r.sum(1) / sw_r.sum(1)
        dk = np.array([
            (swc_d.sum(1) - swc_d[:, k]) / (sw_d.sum(1) - sw_d[:, k])
            - (swc_r.sum(1) - swc_r[:, k]) / (sw_r.sum(1) - sw_r[:, k])
            for k in range(K)])
    return delta, dk


def offsets(delta, sig, rp, anchor):
    """Per-color constants from the anchor bins (inverse-variance)."""
    m = rp >= anchor
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(np.isfinite(sig[:, m]) & (sig[:, m] > 0),
                     1.0 / sig[:, m] ** 2, 0.0)
        return (np.nansum(w * np.nan_to_num(delta[:, m]), 1)
                / np.maximum(w.sum(1), 1e-30))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--npz", required=True)
    ap.add_argument("--anchor", type=float, default=10.0,
                    help="O_c anchor radius, Mpc (ignored with --null)")
    ap.add_argument("--null", default=None,
                    help="star-null paircount npz: subtract its per-color "
                         "profile instead of anchoring constants")
    a = ap.parse_args()

    d = np.load(a.npz, allow_pickle=True)
    rp = np.asarray(d["rp_mid"], float)
    cnames = [str(c) for c in d["colors"]]
    sig = np.asarray(d["sigma"], float)
    delta, dk = loo_deltas(d)
    K = dk.shape[0]

    if a.null:
        n = np.load(a.null, allow_pickle=True)
        nrp, ndel = np.asarray(n["rp_mid"], float), np.asarray(n["delta"])
        prof = np.array([np.interp(np.log(rp), np.log(nrp), ndel[c])
                         for c in range(ndel.shape[0])])
        delta = delta - prof
        dk = dk - prof[None]
        print(f"  null profile subtracted from {a.null}")
    else:
        O = offsets(delta, sig, rp, a.anchor)
        delta = delta - O[:, None]
        for k in range(K):
            dk[k] = dk[k] - offsets(dk[k], sig, rp, a.anchor)[:, None]
        print("  offsets O_c (anchor >= %.1f Mpc): %s" % (
            a.anchor, np.array2string(O, precision=6)))

    use = [cnames.index(c) for c in COLS_USE if c in cnames]
    dR_all = {c: R_COEF[c.split("-")[0]] - R_COEF[c.split("-")[1]]
              for c in cnames}
    nbin = len(rp)
    ebv = np.full(nbin, np.nan)
    err = np.full(nbin, np.nan)
    for b in range(nbin):
        cols = [c for c in use
                if (cnames[c] != "z-W1" or rp[b] >= WISE_RP_MIN)]
        x = dk[:, cols, b]
        ok = np.all(np.isfinite(x), 1)
        k, nc = ok.sum(), len(cols)
        if k < nc + 3 or not np.all(np.isfinite(delta[cols, b])):
            continue
        x = x[ok]
        m = x.mean(0)
        C = (k - 1) / k * (x - m).T @ (x - m)
        hartlap = (k - nc - 2) / (k - 1)
        Cinv = hartlap * np.linalg.pinv(C)
        dR = np.array([dR_all[cnames[c]] for c in cols])
        den = dR @ Cinv @ dR
        ebv[b] = dR @ Cinv @ delta[cols, b] / den
        err[b] = 1.0 / np.sqrt(den)

    print(f"\n== cov-ML result ({os.path.basename(a.npz)}, K={K}):")
    for b in range(nbin):
        print(f"  {rp[b]*1e3:8.0f} kpc  E(B-V) = {ebv[b]:+9.2e} "
              f"+- {err[b]:8.1e}")

    meta = json.loads(str(d["meta"]))
    meta.update(estimator="ml_fullcov",
                offset_model=("null:" + a.null) if a.null
                else f"anchor>={a.anchor}Mpc")
    fo = a.npz.replace("_qso.npz", "ml_qso.npz")
    np.savez(fo, rp_edges=d["rp_edges"], rp_mid=rp, delta=delta, sigma=sig,
             ebv=ebv, ebv_err=err, npair_data=d["npair_data"],
             npair_ran=d["npair_ran"], colors=d["colors"],
             meta=json.dumps(meta))
    print(f"saved -> {fo}")


if __name__ == "__main__":
    main()
