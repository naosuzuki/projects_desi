"""34_mccleary_legacy.py -- McCleary with DESI Legacy data: assessment
and figures for Paths A and B (33_pathab_babysit.sh).

  PATH A  her samples, Legacy photometry: redMaGiC highdens x highlum
          crayons with DR11 (Tractor) griz colors, CSFD, cov-ML
          (paircount_rmgdr11csfdml), against the completed MOF
          reproduction on the SAME galaxies (paircount_rmgmofcsfdml).
          A_V = 3.1 E / f_abs, f_abs = 0.806 (the photo-z absorption
          belongs to the redMaGiC redshifts, not to the photometry).
          Both are ALSO re-solved here with the per-color constants
          anchored at >= 10 Mpc (keys mof10, dr1110) from the raw
          jackknife npz: with DR11 colors the fg-minus-random colors
          carry a coherent negative constant out to 34 Mpc, and the
          anchored solve is the like-for-like test of the profile.
  PATH B  her selection, DESI spectroscopy: BGS passive & 0.15<z<0.45
          & M_r<-20.2 x merged QSOs, lw config, cov-ML
          (paircount_bgsmcc_mqcsfdlwml).  f = 1 (spectroscopic QSOs).
          The lw solve anchors the per-color constants at >= 10 Mpc,
          so its comparison stops at 10 Mpc.

Each profile is compared with her LRG law over r >= 150 kpc exactly as
des/07_assess_mccleary.py does (mid 0.15-1.5 Mpc, tail beyond), plus a
one-parameter amplitude C_fit / 4.6e-3 of the law.  Figures follow
20_mccrepro_fig.py v2 (one plot per figure, n(z) inset):
    mccleary_legacy_a.pdf   MOF vs DR11 colors, same galaxies, anchored
    mccleary_legacy_b.pdf   BGS-selected analog x QSOs

    python3 34_mccleary_legacy.py [NPZDIR]
"""
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
RV, RMIN = 3.1, 150.0

SAMPLES = {  # key -> (npz stem, f_abs, r_max kpc, label)
    "mof": ("rmgmofcsfdml", 0.806, 4.0e4,
            r"redMaGiC $\times$ crayons, DES MOF colors"),
    "dr11": ("rmgdr11csfdml", 0.806, 4.0e4,
             r"same galaxies, Legacy DR11 (Tractor) colors"),
    "bgsmcc": ("bgsmcc_mqcsfdlwml", 1.0, 1.0e4,
               r"BGS passive, $M_r<-20.2$, $0.15<z<0.45$ $\times$ QSOs"),
}


def law_av(r_kpc):                      # McCleary+25 LRG power law
    return 4.6e-3 * (0.91 * r_kpc / 100.0) ** -0.77


def load(key):
    stem, f, rmax, label = SAMPLES[key]
    fn = f"{SP}/paircount_{stem}_qso.npz"
    if not os.path.exists(fn):
        return None
    d = np.load(fn, allow_pickle=True)
    rp = d["rp_mid"] * 1e3
    return dict(rp=rp, av=RV * d["ebv"] / f, err=RV * d["ebv_err"] / f,
                keep=(rp >= RMIN) & (rp <= rmax), label=label, f=f)


def assess(key, s):
    k = s["keep"]
    rp, av, err = s["rp"][k], s["av"][k], s["err"][k]
    L = law_av(rp)
    pull = (av - L) / err
    C = np.sum(av * L / err ** 2) / np.sum(L ** 2 / err ** 2)
    sC = 1.0 / np.sqrt(np.sum(L ** 2 / err ** 2))
    out = dict(C=float(C), sC=float(sC), zones={})
    print(f"\n== {key}: {SAMPLES[key][0]}  (f_abs = {s['f']})")
    print("   rp[kpc]   A_V [1e-3]   +-      law     pull")
    for r, a, e, l, p in zip(rp, av, err, L, pull):
        print(f"  {r:8.0f}   {a*1e3:+8.3f}  {e*1e3:6.3f}  {l*1e3:7.3f}"
              f"  {p:+6.2f}")
    for zone, m in (("mid", rp < 1500), ("tail", rp >= 1500)):
        p = pull[m]
        out["zones"][zone] = dict(n=int(m.sum()), pull=float(p.mean()),
                                  chi2=float((p ** 2).mean()))
        print(f"  zone {zone:>4s}: mean pull {p.mean():+.2f}, chi2/n "
              f"{(p**2).mean():.2f} over {m.sum()} bins")
    print(f"  amplitude: {C:.2f} +- {sC:.2f} x her law "
          f"(C = {C*4.6:.2f} +- {sC*4.6:.2f} e-3)")
    out["chi2_zero"] = float(np.sum((av / err) ** 2))
    out["chi2_law"] = float(np.sum(pull ** 2))
    out["nbin"] = int(k.sum())
    print(f"  chi2 vs zero {out['chi2_zero']:.1f}, vs law "
          f"{out['chi2_law']:.1f}, {out['nbin']} bins")
    if key in ANCHORED:
        # an anchored profile is (true profile - its anchor-zone mean):
        # fit the law minus the law's own anchor-zone mean
        an = s["rp"] >= 1.0e4
        wa = 1.0 / s["err"][an] ** 2
        La = L - np.sum(law_av(s["rp"][an]) * wa) / wa.sum()
        Ca = np.sum(av * La / err ** 2) / np.sum(La ** 2 / err ** 2)
        sa = 1.0 / np.sqrt(np.sum(La ** 2 / err ** 2))
        out.update(C_anch=float(Ca), sC_anch=float(sa))
        print(f"  amplitude vs (law - its anchor mean): "
              f"{Ca:.2f} +- {sa:.2f}")
    return out


def load_anchored(key, stem, f, label, anchor=10.0):
    """cov-ML (des/08) on the raw jackknife npz after subtracting, per
    color, the inverse-variance mean of delta over rp >= anchor Mpc --
    from the data vector and from every leave-one-out vector."""
    fn = f"{SP}/paircount_{stem}_qso.npz"
    if not os.path.exists(fn):
        return None
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ml08", f"{HERE}/des/08_ml_estimator.py")
    ml = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ml)
    d = np.load(fn, allow_pickle=True)
    names = [str(c) for c in d["colors"]]
    dR = np.array([ml.R_DES[a] - ml.R_DES[b]
                   for a, b in (n.split("-") for n in names)])
    rp = d["rp_mid"]
    an = rp >= anchor
    w = 1.0 / d["sigma"][:, an] ** 2
    delta = d["delta"] - ((d["delta"][:, an] * w).sum(1) / w.sum(1))[:, None]
    dk = ml.loo_deltas(d)
    dk = dk - ((dk[:, :, an] * w[None]).sum(2) / w.sum(1)[None])[:, :, None]
    K, nc, nb = dk.shape
    E, Eerr = np.zeros(nb), np.zeros(nb)
    for i in range(nb):
        x = dk[:, :, i]
        m = x.mean(0)
        Ci = ((K - nc - 2) / (K - 1)
              * np.linalg.pinv((K - 1) / K * (x - m).T @ (x - m)))
        den = dR @ Ci @ dR
        E[i], Eerr[i] = dR @ Ci @ delta[:, i] / den, 1.0 / np.sqrt(den)
    SAMPLES[key] = (stem + " [10 Mpc anchor]", f, anchor * 1e3, label)
    rp = rp * 1e3
    return dict(rp=rp, av=RV * E / f, err=RV * Eerr / f,
                keep=(rp >= RMIN) & (rp <= anchor * 1e3), label=label, f=f)


S = {k: load(k) for k in list(SAMPLES)}
S["mof10"] = load_anchored("mof10", "rmgmofcsfd", 0.806,
                           "DES MOF colors, anchored")
S["dr1110"] = load_anchored("dr1110", "rmgdr11csfd", 0.806,
                            "Legacy DR11 colors, anchored")
ANCHORED = ("mof10", "dr1110", "bgsmcc")
RES = {k: assess(k, s) for k, s in S.items() if s is not None}
if S["dr11"] is not None:
    a, b = S["dr11"], S["mof"]
    k = a["keep"]
    print("\n== DR11 colors vs MOF colors, same galaxies (errors are "
          "correlated: same pairs)")
    for r, x, ex, y, ey in zip(a["rp"][k], a["av"][k], a["err"][k],
                               b["av"][k], b["err"][k]):
        print(f"  {r:8.0f}   DR11 {x*1e3:+7.3f} +- {ex*1e3:5.3f}   "
              f"MOF {y*1e3:+7.3f} +- {ey*1e3:5.3f}   ratio {x/y:5.2f}")
with open(f"{SP}/mccleary_legacy.json", "w") as fh:
    json.dump(RES, fh)

# ------------------------------------------------------------ figures ---
with open(f"{HERE}/fgz_stats_2026-11-12.json") as fh:
    FGZ = json.load(fh)
with open(f"{HERE}/bgz_stats_2026-11-12.json") as fh:
    BGZ = json.load(fh)
with open(f"{HERE}/fgz_stats_bgsmcc.json") as fh:
    FGZ.update(json.load(fh))
YMIN, YMAX = 2e-5, 2e-2          # 20_ uses 8e-3; room for the legend


def points(ax, s, color, marker, xoff, label):
    """20_ v2 look: filled symbol at |A_V| for >= 1 sigma bins; bins
    below 1 sigma become 2 sigma upper limits (signs in the caption)."""
    k = s["keep"]
    x, a, e = s["rp"][k] * xoff, s["av"][k], s["err"][k]
    det = np.abs(a) > e
    y = np.abs(a[det])
    lo = np.clip(np.minimum(e[det], y - YMIN * 1.05), 0, None)
    ax.errorbar(x[det], y, yerr=[lo, e[det]], fmt=marker, ms=7,
                color=color, capsize=3, zorder=4, label=label)
    ul = np.abs(a[~det]) + 2 * e[~det]
    ax.errorbar(x[~det], ul, yerr=ul * 0.25, uplims=True, fmt=marker,
                ms=4, mfc="none", color=color, lw=1, alpha=0.7, zorder=4)


def frame(ax):
    rr = np.logspace(np.log10(140), np.log10(40000), 200)
    ax.plot(rr, law_av(rr), "k--", lw=1.6, zorder=3,
            label=r"McCleary et al. (2026) LRG law:"
                  r" $A_V = 4.6\times10^{-3}\,(r/100\,h^{-1}{\rm kpc})"
                  r"^{-0.77}$")
    ax.fill_between(rr, law_av(rr) * 3.1 / 4.6, law_av(rr), color="0.82",
                    zorder=1, label=r"their redMaGiC range $C=(3.1$--$4.6)"
                                    r"\times10^{-3}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(140, 45000)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xlabel(r"$r_p$ [physical kpc]", fontsize=13)
    ax.set_ylabel(r"$A_V$ [mag]", fontsize=13)
    ax.tick_params(labelsize=11)


def inset(ax, entries, zmax):
    ins = ax.inset_axes([0.04, 0.07, 0.27, 0.33])
    for d_, color, label in entries:
        e = np.asarray(d_["edges"], float)
        h = np.asarray(d_["hist"], float)
        h = h / h.max()
        zc = 0.5 * (e[:-1] + e[1:])
        ins.fill_between(zc, h, color=color, alpha=0.16, lw=0)
        ins.plot(zc, h, color=color, lw=1.3,
                 label=f"{label}\n$\\langle z\\rangle{{=}}{d_['mean']:.2f}$,"
                       f" $\\tilde z{{=}}{d_['median']:.2f}$")
    ins.set_xlim(0, zmax)
    ins.set_ylim(0, 1.55)
    ins.set_yticks([])
    ins.set_xlabel("$z$", fontsize=8, labelpad=1)
    ins.tick_params(labelsize=7)
    ins.legend(fontsize=6.6, frameon=False, loc="upper left",
               handlelength=1.2, borderaxespad=0.1, labelspacing=0.45)
    ins.set_title("foreground / background $n(z)$", fontsize=7.5, pad=2)


def zone_text(key, tag=""):
    z = RES[key]["zones"]
    return (f"{tag}mid {z['mid']['pull']:+.2f} ({z['mid']['chi2']:.2f}),  "
            f"tail {z['tail']['pull']:+.2f} ({z['tail']['chi2']:.2f}),  "
            f"amplitude ${RES[key]['C']:.2f}\\pm{RES[key]['sC']:.2f}$"
            " $\\times$ law")


def box(ax, lines):
    head = "vs. her law, mean pull ($\\chi^2/n$) per zone:"
    ax.text(0.345, 0.08, "\n".join([head] + lines),
            transform=ax.transAxes, fontsize=9.4, va="bottom",
            bbox=dict(fc="white", ec="0.6", alpha=0.9))


def save(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/{name}.{ext}", dpi=150)
    print(f"saved -> {FIG}/{name}.pdf")
    plt.close(fig)


def amp_line(tag, key):
    r, r0 = RES[key], RES[key[:-2]]
    return (f"{tag} ${r['C_anch']:+.2f}\\pm{r['sC_anch']:.2f}$ anchored "
            f"($\\chi^2$ vs. zero {r['chi2_zero']:.0f}/{r['nbin']}),  "
            f"${r0['C']:+.2f}\\pm{r0['sC']:.2f}$ as counted")


def shade_anchor(ax):
    ax.axvspan(1.0e4, 45000, color="0.93", zorder=0)
    ax.text(2.1e4, 2.2e-3, "anchor zone\n($\\geq$10 Mpc)", ha="center",
            va="top", fontsize=8.5, color="0.35")


if S["dr1110"] is not None:
    # like-for-like: same galaxies, pairs, estimator and anchor; only
    # the photometry differs
    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    frame(ax)
    points(ax, S["mof10"], "C3", "o", 1.0,
           r"redMaGiC $\times$ DES crayons, DES MOF colors")
    points(ax, S["dr1110"], "C0", "s", 1.07,
           "same galaxies and pairs, Legacy DR11 (Tractor) colors")
    shade_anchor(ax)
    ax.legend(loc="upper right", fontsize=9.5)
    ax.text(0.345, 0.08, "\n".join(
        ["amplitude in units of her law:",
         amp_line("MOF: ", "mof10"), amp_line("DR11:", "dr1110")]),
        transform=ax.transAxes, fontsize=9.0, va="bottom",
        bbox=dict(fc="white", ec="0.6", alpha=0.9))
    inset(ax, ((FGZ["redmagic"], "C4", "fg: redMaGiC highdens"),
               (BGZ["descrayon"], "C0", "bg: highlum crayons")), 1.0)
    save(fig, "mccleary_legacy_a")

fig, ax = plt.subplots(figsize=(10.0, 5.6))
frame(ax)
points(ax, S["bgsmcc"], "C2", "D", 1.0,
       r"DESI BGS passive, $M_r<-20.2$, $0.15<z<0.45$ $\times$ "
       r"spectroscopic QSOs ($f\equiv1$)")
shade_anchor(ax)
ax.legend(loc="upper right", fontsize=9.5)
box(ax, [zone_text("bgsmcc")])
inset(ax, ((FGZ["bgsmcc"], "C2", "fg: BGS passive, $M_r<-20.2$"),
           (BGZ["qso"], "C3", "bg: merged QSOs")), 4.0)
save(fig, "mccleary_legacy_b")
