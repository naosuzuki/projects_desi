"""Investigate the -3 sigma 125 kpc bin of the BGS lw recount.

Per-color deltas there vs the photgal replica and the earlier (pre-lw)
BGS run; leave-one-region-out distribution of the ML estimate to see
whether a few sky regions drive it.
"""
import numpy as np

P = "/Users/suzuki/data/igmdust/paircount"
R = {"g": 3.214, "r": 2.165, "z": 1.211, "W1": 0.184}

def loo(d):
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

def report(name, path):
    d = np.load(path, allow_pickle=True)
    rp = np.asarray(d["rp_mid"], float)
    b = int(np.argmin(np.abs(rp - 0.125)))
    cn = [str(c) for c in d["colors"]]
    delta, dk = loo(d)
    sig = np.asarray(d["sigma"], float)
    print(f"== {name}  (bin {b}: {rp[b]*1e3:.0f} kpc)")
    for i, c in enumerate(cn):
        a_, bnd = c.split("-")
        dr = R.get(a_, np.nan) - R.get(bnd, np.nan)
        s = sig[i, b]
        print(f"  {c:5s} delta={delta[i,b]:+9.2e} +-{s:8.1e} "
              f"({delta[i,b]/s:+5.1f} sig)  dR={dr:+5.2f} "
              f"implied E={delta[i,b]/dr:+9.2e}")
    # LOO spread of the diagonal estimate at this bin (which regions move it)
    cols = [i for i, c in enumerate(cn) if c in ("g-r", "r-z", "z-W1")]
    dR = np.array([R[cn[i].split("-")[0]] - R[cn[i].split("-")[1]]
                   for i in cols])
    w = 1.0 / sig[cols, b] ** 2
    e_full = np.nansum(w * dR * delta[cols, b]) / np.nansum(w * dR ** 2)
    ek = np.array([np.nansum(w * dR * dk[k][cols, b])
                   / np.nansum(w * dR ** 2) for k in range(dk.shape[0])])
    jump = (e_full - ek) * (dk.shape[0] - 1)   # per-region pull on E
    top = np.argsort(jump)[:5]
    print(f"  diag E={e_full:+9.2e}; top-5 negative-pulling regions "
          f"(of {dk.shape[0]}): " +
          ", ".join(f"{k}:{jump[k]:+6.1e}" for k in top))
    print(f"  LOO spread: median jump {np.median(np.abs(jump)):6.1e}, "
          f"max |jump| {np.max(np.abs(jump)):6.1e}")

report("BGS lw recount", f"{P}/paircount_bgs_mqcsfdlw_qso.npz")
report("photgal replica", f"{P}/paircount_photgal_mqcsfdlw_qso.npz")
report("BGS pre-lw (fixwave)", f"{P}/fixwave/paircount_bgs_mqcsfdml_qso.npz")
