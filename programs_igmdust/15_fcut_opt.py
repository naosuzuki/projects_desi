"""fcut_opt.py — is there an optimal faint-quasar cut, or do ivar
weights beat every cut?  (2026-11-11, power-budget follow-up)

Loads the merged QSO catalog with the NEW 02_paircount.load_qso_merged
(returns per-color photometric sigma), applies locus_weight (z-locus
residual colors + per-QSO weights), then compares the error on a mean
color as a function of an r-band faint cut:
    unweighted:  err ~ std(resid) / sqrt(N)      (current estimator)
    weighted:    err ~ 1 / sqrt(sum w)           (new estimator)
normalized to weighted-no-cut = 1.  Also the E(B-V)-combined figure of
merit using the MW-law color responses.  CSFD recorrection skipped
(rms 0.0027 mag, irrelevant for scatter).
"""
import importlib.util
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "pc", f"{HERE}/02_paircount.py")
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)

CAT = ("/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/"
       "qso_merged_dr11.fits")

ra, dec, z, col, valid, cerr = pc.load_qso_merged(CAT)
rcol, wt = pc.locus_weight(z, col, valid, cerr)

# r magnitude for the cut axis (same selection as the loader)
with fits.open(CAT) as h:
    d = h[1].data
    zz = np.asarray(d["Z"], float)
    s = (np.asarray(d["MATCHED"], bool)
         & (np.char.strip(d["LS_TYPE"].astype(str)) == "PSF")
         & np.isfinite(zz) & (zz > 0.3) & (zz < 4.0))
    fr = np.asarray(d["FLUX_R"], float)[s] / np.asarray(d["MWT_R"],
                                                        float)[s]
mr = np.where(fr > 0, 22.5 - 2.5 * np.log10(np.maximum(fr, 1e-12)),
              np.inf)
assert len(mr) == len(z), (len(mr), len(z))

dR = np.array([pc.R_COEF[a] - pc.R_COEF[b] for a, b in pc.COLORS])
cuts = [20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, np.inf]
names = [f"{a}-{b}" for a, b in pc.COLORS]

# reference: weighted, no cut
ref = np.array([wt[valid[:, c], c].sum() for c in range(4)])
err_ref = 1.0 / np.sqrt(ref)
fom_ref = np.sqrt(np.sum(dR ** 2 / err_ref ** 2))

print("\n   r_cut    N_kept   " + "".join(f"  unw({n:>4s}) wt({n:>4s})"
                                          for n in names)
      + "   EBV_unw  EBV_wt   (err rel. to weighted no-cut)")
best = (np.inf, None)
for cut in cuts:
    row = f"  {cut:6.1f}" if np.isfinite(cut) else "    none"
    keep = mr < cut
    row += f"  {keep.sum():9,d}"
    eu, ew = np.zeros(4), np.zeros(4)
    for c in range(4):
        m = keep & valid[:, c]
        if m.sum() < 1000:
            eu[c] = ew[c] = np.inf
            continue
        eu[c] = rcol[m, c].std() / np.sqrt(m.sum())
        ew[c] = 1.0 / np.sqrt(wt[m, c].sum())
        row += f"   {eu[c]/err_ref[c]:7.3f} {ew[c]/err_ref[c]:8.3f}"
    fu = np.sqrt(np.sum(dR ** 2 / eu ** 2))
    fw = np.sqrt(np.sum(dR ** 2 / ew ** 2))
    row += f"   {fom_ref/fu:7.3f} {fom_ref/fw:7.3f}"
    print(row)
    if fom_ref / fu < best[0]:
        best = (fom_ref / fu, cut)

# weight share of the faint end: do faint QSOs matter once weighted?
print("\nweight carried by r > 22 quasars (weighted estimator):")
for c in range(4):
    v = valid[:, c]
    frac = wt[v & (mr > 22), c].sum() / wt[v, c].sum()
    nfrac = (v & (mr > 22)).sum() / v.sum()
    print(f"  {names[c]:>5s}: {frac:6.1%} of weight from {nfrac:5.1%} "
          f"of objects")

print(f"\nbest UNWEIGHTED cut: r < {best[1]}  "
      f"(EBV err {best[0]:.3f}x weighted no-cut)")
print("DONE fcut_opt")
