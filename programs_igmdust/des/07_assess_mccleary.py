#!/usr/bin/env python3
"""
07_assess_mccleary.py — igmdust/DES: automatic pass/fail assessment of
the McCleary et al. reproduction, run as a chained SLURM job so the
verdict is computed and logged even with nobody connected.

Compares a crayon profile (default: the CSFD-recorrected wide run)
against the published LRG power law
    A_V = 4.6e-3 (r / 100 h^-1 kpc)^-0.77
(r comoving h^-1; at the sample's z~0.3 and h=0.7, r_h = 0.91 r_phys),
after the photo-z absorption correction (E_true = E/0.806).  Bins
inside 150 kpc are foreground-light contaminated and excluded.  Zones:
mid-scale 0.15-1.5 Mpc (the open deficit), tail 1.5-34 Mpc (previously
reproduced).  A zone passes when |mean pull| < 2 and chi2/n < 2.
"""
import argparse
import os

import numpy as np

OUT = os.path.expanduser("~/data/igmdust/paircount")
F_ABS = 0.806
RV = 3.1


def law_E(r_kpc):
    return 4.6e-3 * (0.91 * r_kpc / 100.0) ** -0.77 / RV


def assess(npz):
    d = np.load(npz)
    rp = d["rp_mid"] * 1e3
    E = d["ebv"] / F_ABS
    S = d["ebv_err"] / F_ABS
    keep = rp >= 150
    print(f"== assessment of {os.path.basename(npz)} vs McCleary law")
    print("   rp[kpc]   E/f [1e-3]   law [1e-3]   pull")
    pulls = {}
    for r, e, s in zip(rp[keep], E[keep], S[keep]):
        p = (e - law_E(r)) / s
        zone = "mid" if r < 1500 else "tail"
        pulls.setdefault(zone, []).append(p)
        print(f"  {r:8.0f}   {e*1e3:+8.3f}     {law_E(r)*1e3:8.3f}"
              f"   {p:+6.2f}  [{zone}]")
    ok = True
    for zone in ("mid", "tail"):
        p = np.array(pulls.get(zone, []))
        if len(p) == 0:
            continue
        mp, chi2 = p.mean(), (p ** 2).mean()
        zok = abs(mp) < 2 and chi2 < 2
        ok &= zok
        print(f"  zone {zone:>4s}: mean pull {mp:+.2f}, chi2/n {chi2:.2f}"
              f" over {len(p)} bins -> {'PASS' if zok else 'FAIL'}")
    print("MCCLEARY_REPRODUCTION_" + ("PASS" if ok else "FAIL"))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--npz", default=f"{OUT}/paircount_rmgcsfd_qso.npz")
    a = ap.parse_args()
    if not os.path.exists(a.npz):
        alt = f"{OUT}/paircount_rmg_qso.npz"
        print(f"  {a.npz} missing; assessing {alt}")
        a.npz = alt
    raise SystemExit(assess(a.npz))
