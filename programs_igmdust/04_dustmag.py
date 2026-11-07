#!/usr/bin/env python3
"""
04_dustmag.py — igmdust: joint dust + magnification decomposition (route 2).

The mean QSO color excess around foregrounds mixes two signals:
  dust          : delta_c = (R_a - R_b) E(B-V)        (fixed color ratios)
  magnification : an achromatic flux boost mu = 1 + dmu leaves every QSO's
      color unchanged but pulls intrinsically fainter QSOs into the
      flux-limited sample; because QSO color correlates with magnitude the
      SAMPLE mean color shifts by M_c dmu, with M_c = d<c>/d(dmu) measured
      empirically from the catalog by relaxing the r limit by 2.5 log10(mu).

Per r_p bin this fits    delta_c = dR_c E + M_c dmu + O_c    with
  - colors {g-r, r-z, z-W1} only (g-z is their sum: no extra information),
  - z-W1 dropped inside WISE_RP_MIN (6 arcsec blending; ms.tex sec. 6),
  - per-color constants O_c anchored on bins r_p > RP_ANCHOR, where
    E = dmu = 0 is assumed (this suppresses any genuine large-scale tail;
    the halo-scale decomposition is what this script is for).
Inside WISE_RP_MIN only 2 colors remain -> the fit is exactly determined
(0 dof): E and dmu are still measured, there is just no residual chi2 check.

WHERE TO RUN: astro4 (reads the 02 npz + QSO catalog; ~1 min, no SLURM).
    python3 04_dustmag.py --fg lrg
    python3 04_dustmag.py --selftest       # no data needed
Output: table, $OUT/dustmag_<fg>.npz, Delta-chi2 contours dustmag_<fg>.png
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "paircount02", os.path.join(_here, "02_paircount.py"))
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)

BASE_DEFAULT = pc.BASE_DEFAULT
OUT_DEFAULT = pc.OUT_DEFAULT
COLS_USE = ["g-r", "r-z", "z-W1"]   # g-z = (g-r)+(r-z): linearly dependent
WISE_RP_MIN = 0.15                  # Mpc; z-W1 dropped inside (W1 blending)
RP_ANCHOR = 2.0                     # Mpc; offsets O_c from bins beyond this
RLIMS = (22.4, 22.6, 22.8)          # r limits for the response (edge ~23.0)
DMU = 0.2                           # finite-difference magnification step
DR_COEF = np.array([pc.R_COEF[a] - pc.R_COEF[b] for a, b in pc.COLORS])


# ------------------------------------------------------ QSO color response ---
def load_qso_mags(base):
    """QSO colors + r mag; selection IDENTICAL to 02_paircount.load_qso."""
    f = f"{base}/QSO_cat_iron_cumulative_v0.fits"
    with fits.open(f) as h:
        d = h[1]
        z = pc.getcol(d, "Z")
        names = {c.upper() for c in d.columns.names}
        sel = np.isfinite(z) & (z > 0.3) & (z < 4.0)
        if "SURVEY" in names:
            sel &= np.char.strip(d.data["SURVEY"].astype(str)) == "main"
        if "PROGRAM" in names:
            sel &= np.char.strip(d.data["PROGRAM"].astype(str)) == "dark"
        if "BI_CIV" in names:
            bi = np.asarray(d.data["BI_CIV"], dtype=float)
            sel &= ~(np.isfinite(bi) & (bi > 0))
        mags = {}
        for b in ("G", "R", "Z", "W1"):
            flux = np.asarray(d.data[f"FLUX_{b}"], dtype=float)
            if f"MW_TRANSMISSION_{b}" in names:
                flux = flux / np.asarray(d.data[f"MW_TRANSMISSION_{b}"],
                                         dtype=float)
            ok = np.isfinite(flux) & (flux > 0)
            m = np.full(len(flux), np.nan)
            m[ok] = 22.5 - 2.5 * np.log10(flux[ok])
            mags[b if b == "W1" else b.lower()] = m
    colors = np.column_stack([mags[a] - mags[b] for a, b in pc.COLORS])
    valid = np.isfinite(colors)
    colors = np.where(valid, colors, 0.0)
    print(f"  QSO: {sel.sum():,} / {len(z):,} after cuts (same as 02)")
    return colors[sel], valid[sel], mags["r"][sel]


def mag_response(colors, valid, rmag, rlims=RLIMS, dmu=DMU):
    """M_c = d<c>_sample/d(dmu), by relaxing the r limit by 2.5 log10(1+dmu).

    The magnified sample at limit r_lim is the unmagnified catalog cut at
    r_lim + 2.5 log10(1+dmu) (individual colors are unchanged; the counts
    re-weighting is automatic).  The spread over rlims is the systematic.
    """
    resp = np.full((len(rlims), colors.shape[1]), np.nan)
    dm = 2.5 * np.log10(1.0 + dmu)
    for i, rl in enumerate(rlims):
        for c in range(colors.shape[1]):
            v0 = valid[:, c] & (rmag < rl)
            v1 = valid[:, c] & (rmag < rl + dm)
            resp[i, c] = (colors[v1, c].mean() - colors[v0, c].mean()) / dmu
    return resp.mean(0), resp.std(0), resp


# -------------------------------------------------------------------- fit ---
def joint_fit(rp, delta, sigma, mc, cnames):
    """Offsets from the anchor bins, then per-bin 2-parameter linear fits.

    Returns (O, Oerr, rows); each row is a dict for one bin with r_p below
    RP_ANCHOR: E, dmu, their covariance, chi2 for joint / dust-only /
    mag-only models, and the per-color pieces for contour plotting.
    """
    iuse = np.array([cnames.index(c) for c in COLS_USE])
    outer = rp >= RP_ANCHOR
    w = 1.0 / sigma ** 2
    O = (w[:, outer] * delta[:, outer]).sum(1) / w[:, outer].sum(1)
    Oerr = 1.0 / np.sqrt(w[:, outer].sum(1))

    rows = []
    for j in np.where(~outer)[0]:
        cs = np.array([c for c in iuse
                       if not (cnames[c] == "z-W1" and rp[j] < WISE_RP_MIN)])
        d = delta[cs, j] - O[cs]
        s = np.sqrt(sigma[cs, j] ** 2 + Oerr[cs] ** 2)
        A = np.column_stack([DR_COEF[cs], mc[cs]])
        W = 1.0 / s ** 2
        cov = np.linalg.inv((A.T * W) @ A)
        th = cov @ (A.T * W) @ d
        chi2 = (W * (A @ th - d) ** 2).sum()
        one = {}
        for lab, v in (("dust", DR_COEF[cs]), ("mag", mc[cs])):
            amp = (W * v * d).sum() / (W * v ** 2).sum()
            one[lab] = (amp, (W * (v * amp - d) ** 2).sum())
        rows.append(dict(j=j, th=th, cov=cov, chi2=chi2, dof=len(cs) - 2,
                         chi2_dust=one["dust"][1], chi2_mag=one["mag"][1],
                         cs=cs, d=d, s=s, A=A))
    return O, Oerr, rows


def contour_fig(fg, rp, rows, png, picks=(0.047, 0.082, 0.254, 0.447)):
    """Delta-chi2 = 2.30 / 6.17 contours in the E(B-V)-dmu plane."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping contour figure")
        return
    sel = [min(rows, key=lambda r: abs(np.log(rp[r["j"]] / p)))
           for p in picks]
    fig, axes = plt.subplots(1, len(sel), figsize=(3.3 * len(sel), 3.4))
    for ax, row in zip(np.atleast_1d(axes), sel):
        th, se = row["th"], np.sqrt(np.diag(row["cov"]))
        Eg = np.linspace(th[0] - 5 * se[0], th[0] + 5 * se[0], 160)
        Mg = np.linspace(th[1] - 5 * se[1], th[1] + 5 * se[1], 160)
        EE, MM = np.meshgrid(Eg, Mg)
        c2 = np.zeros_like(EE)
        for k in range(len(row["d"])):
            c2 += ((row["A"][k, 0] * EE + row["A"][k, 1] * MM
                    - row["d"][k]) / row["s"][k]) ** 2
        c2 -= c2.min()
        ax.contour(Eg * 1e3, Mg, c2, levels=[2.30, 6.17],
                   colors="C0", linestyles=["-", "--"])
        ax.plot(th[0] * 1e3, th[1], "x", color="C3", ms=8)
        ax.axhline(0, lw=0.5, color="0.6")
        ax.axvline(0, lw=0.5, color="0.6")
        ax.set_title(f"{fg.upper()}  $r_p$ = {rp[row['j']]*1e3:.0f} kpc"
                     f"  ({len(row['d'])} colors)", fontsize=10)
        ax.set_xlabel(r"E(B$-$V)  [$10^{-3}$ mag]")
    np.atleast_1d(axes)[0].set_ylabel(r"$\delta\mu$")
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    print(f"contours -> {png}")


# ----------------------------------------------------------------- driver ---
def run(fg, base, out):
    f = f"{out}/paircount_{fg}_qso.npz"
    d = np.load(f)
    rp, delta, sigma = d["rp_mid"], d["delta"], d["sigma"]
    cnames = [str(c) for c in d["colors"]]
    print(f"== {fg}: {f}")

    print("== magnification response from the QSO catalog")
    colors, valid, rmag = load_qso_mags(base)
    mc, mc_sys, resp = mag_response(colors, valid, rmag)
    for c, m, s in zip(cnames, mc, mc_sys):
        print(f"  M({c:>4s}) = {m:+.4f} +- {s:.4f} (sys over rlims "
              f"{RLIMS})   [dust ratio dR = "
              f"{DR_COEF[cnames.index(c)]:+.3f}]")

    O, Oerr, rows = joint_fit(rp, delta, sigma, mc, cnames)
    print(f"\n== offsets O_c from rp > {RP_ANCHOR} Mpc:  "
          + "  ".join(f"{c} {o:+.2e}+-{e:.1e}"
                      for c, o, e in zip(cnames, O, Oerr)))

    print(f"\n== joint fit: {fg}  (colors {COLS_USE}; z-W1 dropped inside "
          f"{WISE_RP_MIN*1e3:.0f} kpc)")
    print("   rp[kpc] nc   E(B-V)  +-err      dmu   +-err    corr"
          "   chi2/dof   dust-only  mag-only   E_naive")
    naive = d["ebv"]
    for r in rows:
        se = np.sqrt(np.diag(r["cov"]))
        corr = r["cov"][0, 1] / (se[0] * se[1])
        j = r["j"]
        print(f"  {rp[j]*1e3:8.0f} {len(r['d']):2d}  {r['th'][0]:+9.2e}"
              f" {se[0]:8.1e}  {r['th'][1]:+9.2e} {se[1]:8.1e}  {corr:+5.2f}"
              f"   {r['chi2']:5.1f}/{r['dof']}   {r['chi2_dust']:7.1f}"
              f"   {r['chi2_mag']:7.1f}   {naive[j]:+9.2e}")

    os.makedirs(out, exist_ok=True)
    jj = np.array([r["j"] for r in rows])
    fo = f"{out}/dustmag_{fg}.npz"
    np.savez(
        fo, rp=rp[jj],
        E=np.array([r["th"][0] for r in rows]),
        Eerr=np.array([np.sqrt(r["cov"][0, 0]) for r in rows]),
        dmu=np.array([r["th"][1] for r in rows]),
        dmuerr=np.array([np.sqrt(r["cov"][1, 1]) for r in rows]),
        cov=np.array([r["cov"] for r in rows]),
        chi2=np.array([r["chi2"] for r in rows]),
        dof=np.array([r["dof"] for r in rows]),
        chi2_dust=np.array([r["chi2_dust"] for r in rows]),
        chi2_mag=np.array([r["chi2_mag"] for r in rows]),
        O=O, Oerr=Oerr, Mc=mc, Mc_sys=mc_sys, Mc_resp=resp,
        meta=json.dumps(dict(fg=fg, cols=COLS_USE, wise_rp_min=WISE_RP_MIN,
                             rp_anchor=RP_ANCHOR, rlims=RLIMS, dmu_step=DMU,
                             date=time.strftime("%F %T"))))
    print(f"saved -> {fo}")
    contour_fig(fg, rp, rows, f"{out}/dustmag_{fg}.png")
    return rows


# --------------------------------------------------------------- selftest ---
def selftest():
    """Inject known E(B-V), dmu and offsets into synthetic colors; recover.

    Pass criterion: every fitted parameter within 4.5 sigma of the truth
    (the per-bin fits are linear, so pulls should be ~N(0,1))."""
    rng = np.random.default_rng(7)
    _, rp = pc.rp_bins()
    cnames = [f"{a}-{b}" for a, b in pc.COLORS]
    mc = np.array([0.12, -0.04, 0.08, 0.30])
    E = 0.02 * (rp / 0.1) ** -0.8
    mu = 0.05 * (rp / 0.1) ** -1.0
    E[rp >= RP_ANCHOR] = mu[rp >= RP_ANCHOR] = 0.0
    O = np.array([3e-4, -2e-4, 1e-4, -4e-4])
    sigma = np.full((4, len(rp)), 5e-4)
    delta = (DR_COEF[:, None] * E + mc[:, None] * mu + O[:, None]
             + rng.normal(0, 1, sigma.shape) * sigma)
    # g-z is the sum of g-r and r-z by construction in real data; the fit
    # ignores it, so leaving it independent noise here is harmless.
    Ofit, Oerr, rows = joint_fit(rp, delta, sigma, mc, cnames)
    pulls = []
    for r in rows:
        se = np.sqrt(np.diag(r["cov"]))
        pulls += [(r["th"][0] - E[r["j"]]) / se[0],
                  (r["th"][1] - mu[r["j"]]) / se[1]]
        print(f"  rp={rp[r['j']]*1e3:6.0f} kpc  "
              f"E {r['th'][0]:+.4f} (true {E[r['j']]:+.4f})  "
              f"dmu {r['th'][1]:+.4f} (true {mu[r['j']]:+.4f})")
    pulls += list((Ofit - O) / Oerr)
    worst = np.abs(pulls).max()
    print(f"  worst pull: {worst:.2f} sigma over {len(pulls)} parameters")
    if worst < 4.5:
        print("SELFTEST PASS")
        return 0
    print("SELFTEST FAIL")
    return 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fg", help="profile name: bgs, lrg, bgs_mq, bgsdes, "
                                 "... (reads paircount_<fg>_qso.npz)")
    ap.add_argument("--base", default=BASE_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.fg:
        ap.error("--fg required (or --selftest)")
    run(a.fg, a.base, a.out)
