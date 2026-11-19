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
  * her-law amplitude of each run solved one color at a time (a dust
    signal gives the same amplitude in every color).
  * figures/mccleary_legacy_diff.pdf: the per-band magnitude difference
    versus r_p (20_mccrepro_fig v2 style).

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


# ---- single-color solves of the two runs themselves -------------------
# dust reddens every color by the same E(B-V); a photometric term does
# not.  Solve each run color by color (anchored, f_abs applied) and fit
# her law: where do the two reductions agree?
def run_color(stem, cols):
    r_ = np.load(f"{SP}/paircount_{stem}_qso.npz", allow_pickle=True)
    names = [str(c) for c in r_["colors"]]
    ix = [names.index(c) for c in cols]
    dR = np.array([ml.R_DES[a] - ml.R_DES[b]
                   for a, b in (c.split("-") for c in cols)])
    w_ = 1.0 / r_["sigma"][ix][:, an] ** 2

    def anch(x):
        return x - ((x[..., an] * w_).sum(-1) / w_.sum(-1))[..., None]
    delta = anch(r_["delta"][ix])
    dk = anch(ml.loo_deltas(r_)[:, ix, :])
    Kc, nc, nb = dk.shape
    E, S = np.zeros(nb), np.zeros(nb)
    for i in range(nb):
        x = dk[:, :, i]
        m = x.mean(0)
        Ci = ((Kc - nc - 2) / (Kc - 1) * np.linalg.pinv(
            np.atleast_2d((Kc - 1) / Kc * (x - m).T @ (x - m))))
        den = dR @ Ci @ dR
        E[i], S[i] = dR @ Ci @ delta[:, i] / den, 1.0 / np.sqrt(den)
    return amplitude(RV * E / F_ABS, RV * S / F_ABS, ANCHOR * 1e3, True)


RES["single_color"] = {}
print("\n== her-law amplitude color by color (anchored, 219 kpc - 6.3 Mpc)")
print("   colors             DES MOF            Legacy DR11")
for cols in (["i-z"], ["r-i"], ["g-r"], ["g-r", "r-i", "i-z", "g-i"]):
    a = run_color("rmgmofcsfd", cols)
    b = run_color("rmgdr11csfd", cols)
    RES["single_color"]["+".join(cols)] = dict(mof=a[:2], dr11=b[:2])
    print(f"   {'+'.join(cols):16s} {a[0]:+6.2f} +- {a[1]:.2f}     "
          f"{b[0]:+6.2f} +- {b[1]:.2f}")

with open(f"{SP}/mofdr11_diff.json", "w") as fh:
    json.dump(RES, fh)

# ------------------------------------------------------------- figure ---
# per-band magnitude difference: which band, which radii.  Beyond 71 kpc
# every bin short of the anchor zone is positive (MOF fainter than
# Tractor near foregrounds), so the house log-log |.| axes lose nothing.
YMIN, YMAX = 0.01, 60.0                          # mmag
k = (rp >= RMIN) & (rp <= ANCHOR * 1e3)
fig, ax = plt.subplots(figsize=(10.0, 5.6))
rr = np.logspace(np.log10(140), np.log10(40000), 200)
ax.plot(rr, law_av(rr) * 1e3, "k--", lw=1.6, zorder=3,
        label=r"for scale: the McCleary et al. (2026) LRG law, $A_V(r_p)$")
for c, band, color, marker, xoff in ((4, "g", "C2", "o", 0.94),
                                     (5, "r", "C3", "s", 0.98),
                                     (6, "i", "C1", "D", 1.02),
                                     (7, "z", "C5", "^", 1.06)):
    y, e = da[c][k] * 1e3, ea[c][k] * 1e3
    det = np.abs(y) > e
    lo = np.clip(np.minimum(e[det], np.abs(y[det]) - YMIN * 1.05), 0, None)
    ax.errorbar(rp[k][det] * xoff, np.abs(y[det]), yerr=[lo, e[det]],
                fmt=marker, ms=6.5, color=color, capsize=3, zorder=4,
                label=f"${band}$ band")
    ul = np.abs(y[~det]) + 2 * e[~det]
    ax.errorbar(rp[k][~det] * xoff, ul, yerr=ul * 0.25, uplims=True,
                fmt=marker, ms=4, mfc="none", color=color, lw=1, alpha=0.7)
ax.axvspan(ANCHOR * 1e3, 45000, color="0.93", zorder=0)
ax.text(2.1e4, 2.0, "anchor zone\n($\\geq$10 Mpc)", ha="center",
        va="top", fontsize=8.5, color="0.35")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(140, 45000)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel(r"$r_p$ [physical kpc]", fontsize=13)
ax.set_ylabel(r"$|\Delta m|$ around foregrounds minus around randoms "
              "[mmag]", fontsize=12)
ax.tick_params(labelsize=11)
ax.legend(loc="upper right", fontsize=9.5,
          title="DES MOF minus Legacy DR11 (Tractor) magnitude of the "
                "same crayons", title_fontsize=9.5)
sc = RES["single_color"]
ax.text(0.04, 0.07,
        "her-law amplitude from one color at a time (MOF | Tractor):\n"
        + "\n".join(f"${c}$:  ${sc[c]['mof'][0]:+.2f}\\pm{sc[c]['mof'][1]:.2f}$"
                    f"  |  ${sc[c]['dr11'][0]:+.2f}\\pm{sc[c]['dr11'][1]:.2f}$"
                    for c in ("i-z", "r-i", "g-r")),
        transform=ax.transAxes, fontsize=9.4, va="bottom",
        bbox=dict(fc="white", ec="0.6", alpha=0.9))
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/mccleary_legacy_diff.{ext}", dpi=150)
print(f"saved -> {FIG}/mccleary_legacy_diff.pdf")
