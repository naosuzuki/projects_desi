"""35_mofdr11_diff.py -- the differential test of Path A: where do DES
MOF and Legacy DR11 (Tractor) photometry of the SAME crayons part?

Input: paircount_rmgmofmdr11csfd_qso.npz from
    des/01_crayon_paircount.py --fg redmagic --mof --diff-colors ...
whose eight data columns are per-crayon (MOF minus DR11) quantities
stacked around the redMaGiC foregrounds minus around randoms:
    g-r, r-i, i-z, g-i   locus-residual differences
    g, r, i, z           dereddened magnitude differences

Outputs
  * Delta A_V(r_p) = 3.1 Delta E / f_abs: the cov-ML dust projection of
    the color-difference vector (des/08 estimator), i.e. the part of
    the MOF signal that the Tractor colors of the same galaxies do not
    show; as counted and with the constants anchored at >= 10 Mpc; its
    amplitude in units of her law.
  * per-color and per-band differential profiles with jackknife errors
    (which band, which radii), anchored.
  * figures/mccleary_legacy_diff.pdf (20_mccrepro_fig v2 style).

    python3 35_mofdr11_diff.py [NPZDIR]
"""
import importlib.util
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SP = (sys.argv[1] if len(sys.argv) > 1
      else os.path.expanduser("~/data/igmdust/paircount"))
FIG = "/Users/suzuki/github/papers/26_igmreddening_study/figures"
RV, F_ABS, RMIN, ANCHOR = 3.1, 0.806, 150.0, 10.0
NAMES = ["g-r", "r-i", "i-z", "g-i", "g", "r", "i", "z"]

spec = importlib.util.spec_from_file_location(
    "ml08", f"{HERE}/des/08_ml_estimator.py")
ml = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ml)


def law_av(r_kpc):
    return 4.6e-3 * (0.91 * r_kpc / 100.0) ** -0.77


d = np.load(f"{SP}/paircount_rmgmofmdr11csfd_qso.npz", allow_pickle=True)
rp = d["rp_mid"] * 1e3
an = d["rp_mid"] >= ANCHOR
dk0 = ml.loo_deltas(d)                          # (K, 8, nbin)
K = dk0.shape[0]
w = 1.0 / d["sigma"][:, an] ** 2


def anchored(x):
    """subtract the per-column ivar mean over the anchor bins"""
    return x - ((x[..., an] * w).sum(-1) / w.sum(-1))[..., None]


def jk_err(x):
    return np.sqrt((K - 1) / K * ((x - x.mean(0)) ** 2).sum(0))


def ml_dust(delta, dk):
    """cov-ML projection of the four color columns on the dust vector"""
    dR = np.array([ml.R_DES[a] - ml.R_DES[b]
                   for a, b in (n.split("-") for n in NAMES[:4])])
    nb = delta.shape[1]
    E, S = np.zeros(nb), np.zeros(nb)
    for i in range(nb):
        x = dk[:, :4, i]
        m = x.mean(0)
        Ci = ((K - 4 - 2) / (K - 1)
              * np.linalg.pinv((K - 1) / K * (x - m).T @ (x - m)))
        den = dR @ Ci @ dR
        E[i], S[i] = dR @ Ci @ delta[:4, i] / den, 1.0 / np.sqrt(den)
    return RV * E / F_ABS, RV * S / F_ABS


def amplitude(av, err, rmax, sub_anchor):
    k = (rp >= RMIN) & (rp <= rmax)
    L = law_av(rp[k])
    if sub_anchor:
        wa = 1.0 / err[an] ** 2
        L = L - np.sum(law_av(rp[an]) * wa) / wa.sum()
    C = np.sum(av[k] * L / err[k] ** 2) / np.sum(L ** 2 / err[k] ** 2)
    sC = 1.0 / np.sqrt(np.sum(L ** 2 / err[k] ** 2))
    return float(C), float(sC), float(np.sum((av[k] / err[k]) ** 2)), int(
        k.sum())


RES = {}
print("== Delta A_V = A_V(MOF) - A_V(DR11), same crayons and pairs "
      "[1e-3 mag]")
for tag, delta, dk, rmax, sub in (
        ("as counted", d["delta"], dk0, 4.0e4, False),
        ("anchored", anchored(d["delta"]), anchored(dk0), ANCHOR * 1e3,
         True)):
    av, err = ml_dust(delta, dk)
    C, sC, chi0, n = amplitude(av, err, rmax, sub)
    RES[tag] = dict(rp=rp.tolist(), av=av.tolist(), err=err.tolist(),
                    C=C, sC=sC, chi2_zero=chi0, nbin=n)
    print(f"-- {tag}: amplitude {C:+.2f} +- {sC:.2f} x her law; chi2 vs "
          f"zero {chi0:.1f}/{n}")
    for r, a, e in zip(rp, av, err):
        if r >= 10:
            print(f"   {r:8.0f} kpc  {a*1e3:+8.3f} +- {e*1e3:6.3f}  "
                  f"({a/e:+5.1f} sig)   law {law_av(r)*1e3:6.3f}")

print("\n== per-column differential profiles, anchored at >= 10 Mpc "
      "[1e-4 mag] (MOF minus DR11; fg minus random)")
da, dka = anchored(d["delta"]), anchored(dk0)
ea = jk_err(dka)
RES["columns"] = dict(names=NAMES, delta=da.tolist(), err=ea.tolist(),
                      const=((d["delta"][:, an] * w).sum(1)
                             / w.sum(1)).tolist())
print("   rp[kpc] " + "".join(f"{n:>16s}" for n in NAMES))
for i, r in enumerate(rp):
    if r >= 20:
        print(f"  {r:8.0f} " + "".join(
            f"  {da[c, i]*1e4:+6.2f}({da[c, i]/ea[c, i]:+4.1f})"
            for c in range(8)))
print("   anchor constants [1e-4 mag]: " + "  ".join(
    f"{n} {c*1e4:+.2f}" for n, c in zip(NAMES, RES["columns"]["const"])))
with open(f"{SP}/mofdr11_diff.json", "w") as fh:
    json.dump(RES, fh)

# ------------------------------------------------------------- figure ---
YMIN, YMAX = 2e-5, 2e-2
r = RES["anchored"]
av, err = np.array(r["av"]), np.array(r["err"])
k = (rp >= RMIN) & (rp <= ANCHOR * 1e3)
fig, ax = plt.subplots(figsize=(10.0, 5.6))
rr = np.logspace(np.log10(140), np.log10(40000), 200)
ax.plot(rr, law_av(rr), "k--", lw=1.6, zorder=3,
        label=r"McCleary et al. (2026) LRG law:"
              r" $A_V = 4.6\times10^{-3}\,(r/100\,h^{-1}{\rm kpc})^{-0.77}$")
ax.fill_between(rr, law_av(rr) * 3.1 / 4.6, law_av(rr), color="0.82",
                zorder=1, label=r"their redMaGiC range $C=(3.1$--$4.6)"
                                r"\times10^{-3}$")
det = np.abs(av[k]) > err[k]
y = np.abs(av[k][det])
lo = np.clip(np.minimum(err[k][det], y - YMIN * 1.05), 0, None)
ax.errorbar(rp[k][det], y, yerr=[lo, err[k][det]], fmt="D", ms=7,
            color="C4", capsize=3, zorder=4,
            label=r"$\Delta A_V$: DES MOF minus Legacy DR11 (Tractor), "
                  "per crayon, same pairs")
ul = np.abs(av[k][~det]) + 2 * err[k][~det]
ax.errorbar(rp[k][~det], ul, yerr=ul * 0.25, uplims=True, fmt="D", ms=4,
            mfc="none", color="C4", lw=1, alpha=0.7, zorder=4)
ax.axvspan(ANCHOR * 1e3, 45000, color="0.93", zorder=0)
ax.text(2.1e4, 2.2e-3, "anchor zone\n($\\geq$10 Mpc)", ha="center",
        va="top", fontsize=8.5, color="0.35")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(140, 45000)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel(r"$r_p$ [physical kpc]", fontsize=13)
ax.set_ylabel(r"$|\Delta A_V|$ [mag]", fontsize=13)
ax.tick_params(labelsize=11)
ax.legend(loc="upper right", fontsize=9.5)
r0 = RES["as counted"]
ax.text(0.04, 0.08,
        "amplitude of the difference, in units of her law:\n"
        f"${r['C']:+.2f}\\pm{r['sC']:.2f}$ anchored "
        f"($\\chi^2$ vs. zero {r['chi2_zero']:.0f}/{r['nbin']}),  "
        f"${r0['C']:+.2f}\\pm{r0['sC']:.2f}$ as counted",
        transform=ax.transAxes, fontsize=9.6, va="bottom",
        bbox=dict(fc="white", ec="0.6", alpha=0.9))
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_legacy_diff.{ext}", dpi=150)
print(f"saved -> {FIG}/mccleary_legacy_diff.pdf")
