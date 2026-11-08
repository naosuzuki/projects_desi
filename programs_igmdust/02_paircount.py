#!/usr/bin/env python3
"""
02_paircount.py — igmdust: the galaxy-dust reddening estimator (Menard+10).

Measures the mean background-QSO color excess <delta_c>(r_p) around DESI DR1
foreground galaxies (BGS_ANY or LRG LSS clustering catalogs), relative to the
same measurement around the LSS randoms, in bins of projected PHYSICAL radius
r_p at the foreground redshift.  Combines the colors into E(B-V)(r_p) with
MW-law extinction coefficients, and estimates errors by spatial jackknife.

    delta_c(r_p) = <c>_pairs(fg) - <c>_pairs(random),   c in {g-r, r-z, g-z, z-W1}
    E(B-V)(r_p)  = sum_c w_c (R_a-R_b) delta_c / sum_c w_c (R_a-R_b)^2,
                   w_c = 1/sigma_c^2 (jackknife)

WHERE TO RUN: astro4 (catalogs on /data/astrofs2_1; never run jobs on astrofs4).
    nohup nice -n 19 python3 02_paircount.py --fg bgs > pc_bgs.log 2>&1 &

Usage:
    python3 02_paircount.py --selftest          # no data needed: inject + recover
    python3 02_paircount.py --fg bgs            # QSO colors around BGS (pilot)
    python3 02_paircount.py --fg lrg
    python3 02_paircount.py --fg bgs --max-fg 200000   # quick subsampled pass

Inputs (from 00_download_dr1_catalogs.sh):
    $BASE/{BGS_ANY,LRG}_{NGC,SGC}_clustering.dat.fits        foregrounds
    $BASE/{BGS_ANY,LRG}_{NGC,SGC}_{0,1}_clustering.ran.fits  randoms
    $BASE/QSO_cat_iron_cumulative_v0.fits                    backgrounds
Output:
    $OUT/paircount_<fg>_qso.npz  + table to stdout
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.cosmology import Planck18 as cosmo
from scipy.spatial import cKDTree

# ----------------------------------------------------------------- config ---
BASE_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/dr1"
OUT_DEFAULT = os.path.expanduser("~/data/igmdust/paircount")

RP_MIN, RP_MAX, NBIN_PER_DEX = 0.02, 10.0, 4          # Mpc (physical)
Z_BUFFER = 0.2                                        # z_bg > z_fg + buffer
FG_ZRANGE = {"bgs": (0.05, 0.45), "lrg": (0.40, 1.10)}
ZSLICE = 0.025                                        # fg z-slice for tree queries
CHUNK = 50_000                                        # fg per tree query chunk

# extinction coefficients (R_V=3.1 MW law; as in 01_feasibility.py)
R_COEF = {"g": 3.214, "r": 2.165, "z": 1.211, "W1": 0.184}
COLORS = [("g", "r"), ("r", "z"), ("g", "z"), ("z", "W1")]

# jackknife regions: healpy nside if available, else RA/sin(dec) grid
JK_NSIDE = 8
JK_GRID_DEG = 7.5


# -------------------------------------------------------------- utilities ---
def unitvec(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    c = np.cos(dec)
    return np.column_stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)])


def da_interp():
    """Fast z -> angular-diameter distance [Mpc] interpolator."""
    zg = np.linspace(1e-3, 4.0, 800)
    dg = cosmo.angular_diameter_distance(zg).value
    return lambda z: np.interp(z, zg, dg)


def assign_regions(ra, dec):
    """Jackknife region label per object (healpy pixel or coarse grid)."""
    try:
        import healpy as hp
        return hp.ang2pix(JK_NSIDE, ra, dec, lonlat=True)
    except ImportError:
        i = np.floor(ra / JK_GRID_DEG).astype(np.int64)
        j = np.floor((np.sin(np.radians(dec)) + 1.0)
                     / (JK_GRID_DEG / 57.3)).astype(np.int64)
        return i * 1000 + j


def getcol(hdu, *names):
    cols = {c.upper() for c in hdu.columns.names}
    for n in names:
        if n.upper() in cols:
            return np.asarray(hdu.data[n])
    raise KeyError(f"none of {names} in {hdu.columns.names}")


# ------------------------------------------------------------ data loading ---
def load_lss(base, sample, kind="dat", nran=2):
    """Concatenate NGC+SGC clustering catalogs (data or randoms)."""
    stem = {"bgs": "BGS_ANY", "lrg": "LRG"}[sample]
    files = []
    for cap in ("NGC", "SGC"):
        if kind == "dat":
            files.append(f"{base}/{stem}_{cap}_clustering.dat.fits")
        else:
            files += [f"{base}/{stem}_{cap}_{i}_clustering.ran.fits"
                      for i in range(nran)]
    ra, dec, z, w = [], [], [], []
    for f in files:
        with fits.open(f) as h:
            d = h[1]
            ra.append(getcol(d, "RA"))
            dec.append(getcol(d, "DEC"))
            z.append(getcol(d, "Z"))
            try:
                w.append(getcol(d, "WEIGHT"))
            except KeyError:
                w.append(np.ones(len(ra[-1])))
        print(f"  {os.path.basename(f)}: {len(ra[-1]):,} rows")
    return (np.concatenate(ra), np.concatenate(dec),
            np.concatenate(z), np.concatenate(w))


def load_qso_merged(path):
    """Merged background (05_qso_merge + 06_sweepmatch): DESI + SDSS +
    Quaia QSOs with uniform LS DR11 photometry.  PSF morphology only;
    a color is valid when both bands have positive flux AND ivar.
    NOTE for 04: M_c there is computed from the DESI selection alone —
    acceptable around BGS where dmu is consistent with zero, wrong for
    precision magnification work with this mixed-selection sample."""
    with fits.open(path) as h:
        d = h[1].data
        ra = np.asarray(d["RA"], float)
        dec = np.asarray(d["DEC"], float)
        z = np.asarray(d["Z"], float)
        src = np.asarray(d["SRC"], int)
        sel = (np.asarray(d["MATCHED"], bool)
               & (np.char.strip(d["LS_TYPE"].astype(str)) == "PSF")
               & np.isfinite(z) & (z > 0.3) & (z < 4.0))
        mags = {}
        for b in ("G", "R", "Z", "W1"):
            flux = np.asarray(d[f"FLUX_{b}"], float)
            ivar = np.asarray(d[f"IVAR_{b}"], float)
            mwt = np.asarray(d[f"MWT_{b}"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = flux / mwt
            ok = (np.isfinite(flux) & (flux > 0)
                  & np.isfinite(ivar) & (ivar > 0))
            m = np.full(len(flux), np.nan)
            m[ok] = 22.5 - 2.5 * np.log10(flux[ok])
            mags[b if b == "W1" else b.lower()] = m
    colors = np.column_stack([mags[a] - mags[b] for a, b in COLORS])
    valid = np.isfinite(colors)
    colors = np.where(valid, colors, 0.0)
    cnt = np.bincount(src[sel], minlength=3)
    print(f"  merged QSO: {sel.sum():,} / {len(ra):,} after cuts "
          f"(DESI {cnt[0]:,}, SDSS {cnt[1]:,}, Quaia {cnt[2]:,}); "
          f"validity {valid[sel].mean(0).round(3)}")
    return ra[sel], dec[sel], z[sel], colors[sel], valid[sel]


def load_qso(base):
    """QSO backgrounds: positions, z, MW-dereddened colors + validity."""
    f = f"{base}/QSO_cat_iron_cumulative_v0.fits"
    with fits.open(f) as h:
        d = h[1]
        ra = getcol(d, "TARGET_RA", "RA")
        dec = getcol(d, "TARGET_DEC", "DEC")
        z = getcol(d, "Z")
        names = {c.upper() for c in d.columns.names}

        sel = np.isfinite(z) & (z > 0.3) & (z < 4.0)
        if "SURVEY" in names:
            sel &= np.char.strip(d.data["SURVEY"].astype(str)) == "main"
        if "PROGRAM" in names:
            sel &= np.char.strip(d.data["PROGRAM"].astype(str)) == "dark"
        for balcol in ("BI_CIV",):                      # BAL cut if available
            if balcol in names:
                bi = np.asarray(d.data[balcol], dtype=float)
                sel &= ~(np.isfinite(bi) & (bi > 0))
                print(f"  BAL cut on {balcol}: keeping {sel.sum():,}")

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

    colors = np.column_stack([mags[a] - mags[b] for a, b in COLORS])
    valid = np.isfinite(colors)
    colors = np.where(valid, colors, 0.0)
    print(f"  QSO: {sel.sum():,} / {len(ra):,} after cuts; "
          f"color validity {valid[sel].mean(0).round(3)}")
    return ra[sel], dec[sel], z[sel], colors[sel], valid[sel]


# -------------------------------------------------------------- estimator ---
def rp_bins():
    nb = int(np.ceil(NBIN_PER_DEX * np.log10(RP_MAX / RP_MIN)))
    edges = np.logspace(np.log10(RP_MIN), np.log10(RP_MAX), nb + 1)
    return edges, np.sqrt(edges[:-1] * edges[1:])


def accumulate(fg_ra, fg_dec, fg_z, fg_w, fg_reg, bg_tree, bg_vec, bg_z,
               colors, valid, edges, nreg, reg_index, label=""):
    """Stack sum(w) and sum(w*c) per (jackknife region, r_p bin, color)."""
    ncol, nbin = colors.shape[1], len(edges) - 1
    sumw = np.zeros((ncol, nreg * nbin))
    sumwc = np.zeros((ncol, nreg * nbin))
    npair = np.zeros(nreg * nbin)
    da = da_interp()
    logedges = np.log(edges)

    order = np.argsort(fg_z)
    zlo, zhi = fg_z.min(), fg_z.max()
    t0, done = time.time(), 0
    for za in np.arange(zlo, zhi, ZSLICE):
        i0, i1 = np.searchsorted(fg_z[order], [za, za + ZSLICE])
        if i1 == i0:
            continue
        idx = order[i0:i1]
        d_slice = da(np.maximum(fg_z[idx], 1e-3))
        theta_max = RP_MAX / d_slice.min()
        chord = 2 * np.sin(min(theta_max, np.pi / 2) / 2)
        vec = unitvec(fg_ra[idx], fg_dec[idx])
        for c0 in range(0, len(idx), CHUNK):
            sl = slice(c0, min(c0 + CHUNK, len(idx)))
            try:
                nbrs = bg_tree.query_ball_point(vec[sl], chord, workers=-1)
            except TypeError:                            # older scipy
                nbrs = bg_tree.query_ball_point(vec[sl], chord)
            cnt = np.fromiter((len(n) for n in nbrs), int, len(nbrs))
            if cnt.sum() == 0:
                continue
            fi = np.repeat(np.arange(len(nbrs)), cnt)
            bi = np.concatenate([np.asarray(n, int) for n in nbrs if n])
            gi = idx[sl][fi]                             # global fg index
            cosang = np.clip(np.einsum("ij,ij->i", vec[sl][fi], bg_vec[bi]),
                             -1, 1)
            rp = np.arccos(cosang) * d_slice[sl][fi]
            keep = ((rp >= RP_MIN) & (rp < RP_MAX)
                    & (bg_z[bi] > fg_z[gi] + Z_BUFFER))
            fi, bi, gi, rp = fi[keep], bi[keep], gi[keep], rp[keep]
            if len(rp) == 0:
                continue
            b = np.clip(np.searchsorted(logedges, np.log(rp), "right") - 1,
                        0, nbin - 1)
            key = reg_index[fg_reg[gi]] * nbin + b
            w = fg_w[gi]
            np.add.at(npair, key, 1.0)
            for c in range(ncol):
                wv = w * valid[bi, c]
                np.add.at(sumw[c], key, wv)
                np.add.at(sumwc[c], key, wv * colors[bi, c])
        done += i1 - i0
        print(f"  [{label}] z<{za+ZSLICE:.3f}  fg done {done:,}/{len(fg_z):,}"
              f"  pairs {npair.sum():,.0f}  {time.time()-t0:.0f}s",
              flush=True)
    shape = (ncol, nreg, nbin)
    return sumw.reshape(shape), sumwc.reshape(shape), npair.reshape(nreg, nbin)


def jackknife(dat, ran):
    """delta_c per bin + jackknife sigma from leave-one-region-out."""
    sw_d, swc_d, _ = dat
    sw_r, swc_r, _ = ran
    with np.errstate(invalid="ignore", divide="ignore"):
        delta = (swc_d.sum(1) / sw_d.sum(1)) - (swc_r.sum(1) / sw_r.sum(1))
        nreg = sw_d.shape[1]
        dk = np.array([
            (swc_d.sum(1) - swc_d[:, k]) / (sw_d.sum(1) - sw_d[:, k])
            - (swc_r.sum(1) - swc_r[:, k]) / (sw_r.sum(1) - sw_r[:, k])
            for k in range(nreg)])                       # (K, ncol, nbin)
        good = np.isfinite(dk)
        n_eff = good.sum(0)
        mean_k = np.nanmean(np.where(good, dk, np.nan), 0)
        var = (n_eff - 1) / np.maximum(n_eff, 2) * np.nansum(
            np.where(good, (dk - mean_k) ** 2, 0), 0)
    return delta, np.sqrt(var)


def solve_ebv(delta, sig):
    """Per-bin E(B-V) from the color vector (inverse-variance weights)."""
    dR = np.array([R_COEF[a] - R_COEF[b] for a, b in COLORS])[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = 1.0 / sig ** 2
        num = np.nansum(w * dR * delta, 0)
        den = np.nansum(w * dR ** 2, 0)
        return num / den, 1.0 / np.sqrt(den)


# ------------------------------------------------------------------ driver ---
def run(fg_sample, base, out, max_fg, nran, tag="", qso_cat=None):
    edges, rp_mid = rp_bins()
    print(f"== foregrounds: {fg_sample}")
    fra, fdec, fz, fw = load_lss(base, fg_sample, "dat")
    print("== randoms")
    rra, rdec, rz, rw = load_lss(base, fg_sample, "ran", nran=nran)
    print("== backgrounds: QSO")
    if qso_cat:
        qra, qdec, qz, qcol, qval = load_qso_merged(qso_cat)
    else:
        qra, qdec, qz, qcol, qval = load_qso(base)
    return run_arrays(fg_sample, (fra, fdec, fz, fw), (rra, rdec, rz, rw),
                      (qra, qdec, qz, qcol, qval), edges, rp_mid,
                      out, max_fg, tag)


def run_arrays(name, fg, ran, bg, edges, rp_mid, out, max_fg=0, tag=""):
    fra, fdec, fz, fw = fg
    rra, rdec, rz, rw = ran
    qra, qdec, qz, qcol, qval = bg

    zlo, zhi = FG_ZRANGE.get(name, (fz.min(), fz.max()))
    fcut = (fz >= zlo) & (fz <= zhi)
    rcut = (rz >= zlo) & (rz <= zhi)
    fra, fdec, fz, fw = fra[fcut], fdec[fcut], fz[fcut], fw[fcut]
    rra, rdec, rz, rw = rra[rcut], rdec[rcut], rz[rcut], rw[rcut]
    if max_fg and len(fz) > max_fg:
        rng = np.random.default_rng(1)
        pick = rng.choice(len(fz), max_fg, replace=False)
        fra, fdec, fz, fw = fra[pick], fdec[pick], fz[pick], fw[pick]
        rmax = min(len(rz), 5 * max_fg)                  # keep randoms ~5x fg
        rpick = rng.choice(len(rz), rmax, replace=False)
        rra, rdec, rz, rw = rra[rpick], rdec[rpick], rz[rpick], rw[rpick]
        print(f"  subsampled fg to {max_fg:,}, randoms to {rmax:,}")
    print(f"  fg {len(fz):,}  ran {len(rz):,}  bg {len(qz):,}"
          f"  z_fg [{zlo},{zhi}]")

    freg = assign_regions(fra, fdec)
    rreg = assign_regions(rra, rdec)
    labels = np.unique(np.concatenate([freg, rreg]))
    nreg = len(labels)
    lut = np.zeros(labels.max() + 1, dtype=int)
    lut[labels] = np.arange(nreg)
    print(f"  jackknife regions: {nreg}")

    bg_vec = unitvec(qra, qdec)
    tree = cKDTree(bg_vec)
    dat = accumulate(fra, fdec, fz, fw, freg, tree, bg_vec, qz, qcol, qval,
                     edges, nreg, lut, label="data")
    rnd = accumulate(rra, rdec, rz, rw, rreg, tree, bg_vec, qz, qcol, qval,
                     edges, nreg, lut, label="random")

    delta, sig = jackknife(dat, rnd)
    ebv, ebv_err = solve_ebv(delta, sig)
    npair_d, npair_r = dat[2].sum(0), rnd[2].sum(0)

    cnames = [f"{a}-{b}" for a, b in COLORS]
    print(f"\n== result: {name}{tag}  (delta color = fg - random, mag)")
    hdr = "   rp[kpc]   npairs_d " + "".join(
        f"  d({c:>5s}) +-sig    " for c in cnames) + "   E(B-V)  +-err"
    print(hdr)
    for i, rp in enumerate(rp_mid):
        row = f"  {rp*1e3:8.0f} {npair_d[i]:11.3e}"
        for c in range(len(COLORS)):
            row += f"  {delta[c, i]:+9.2e} {sig[c, i]:8.1e}"
        row += f"  {ebv[i]:+9.2e} {ebv_err[i]:8.1e}"
        print(row)

    os.makedirs(out, exist_ok=True)
    fo = f"{out}/paircount_{name}{tag}_qso.npz"
    np.savez(fo, rp_edges=edges, rp_mid=rp_mid, delta=delta, sigma=sig,
             ebv=ebv, ebv_err=ebv_err, npair_data=npair_d, npair_ran=npair_r,
             colors=np.array(cnames),
             # per-region accumulations: enable covariance-weighted
             # re-estimation (color-color covariance from jackknife)
             jk_sw_d=dat[0], jk_swc_d=dat[1],
             jk_sw_r=rnd[0], jk_swc_r=rnd[1],
             meta=json.dumps(dict(fg=name, zrange=[float(zlo), float(zhi)],
                                  zbuffer=Z_BUFFER, nreg=int(nreg),
                                  date=time.strftime("%F %T"))))
    print(f"\nsaved -> {fo}")
    return rp_mid, delta, sig, ebv, ebv_err


# ------------------------------------------------------------- self test ----
def selftest():
    """Synthetic end-to-end check: inject a Menard-like profile, recover it.

    10x10 deg patch; amplitude boosted x50 over reality so the recovery is
    high-S/N; rp limited to < 3 Mpc for speed. Passes if the stacked E(B-V)
    matches the injected profile to 15% where S/N > 5.
    """
    global RP_MAX
    rp_max_save, RP_MAX = RP_MAX, 3.0
    rng = np.random.default_rng(42)
    amp, slope = 0.06, -0.86                     # injected E(B-V) at 100 kpc

    nfg, nbg, nrn = 20_000, 100_000, 200_000
    def scatter(n):
        ra = rng.uniform(150, 160, n)
        dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(-5)),
                                               np.sin(np.radians(5)), n)))
        return ra, dec
    fra, fdec = scatter(nfg)
    fz = rng.uniform(0.15, 0.30, nfg)
    qra, qdec = scatter(nbg)
    qz = rng.uniform(0.8, 2.5, nbg)
    rra, rdec = scatter(nrn)
    rz = rng.permutation(np.repeat(fz, 10))

    # inject: per QSO, sum E(B-V) from all fg within RP_MAX, then add to colors
    da = da_interp()
    fvec, qvec = unitvec(fra, fdec), unitvec(qra, qdec)
    qtree = cKDTree(qvec)
    ebv_q = np.zeros(nbg)
    dfg = da(fz)
    chord = 2 * np.sin(RP_MAX / dfg.min() / 2)
    nbrs = qtree.query_ball_point(fvec, chord)
    for i, nn in enumerate(nbrs):
        nn = np.asarray(nn, int)
        if len(nn) == 0:
            continue
        ang = np.arccos(np.clip(qvec[nn] @ fvec[i], -1, 1))
        rp = np.maximum(ang * dfg[i], 0.02)
        m = (rp < RP_MAX) & (qz[nn] > fz[i] + Z_BUFFER)
        ebv_q[nn[m]] += amp * (rp[m] / 0.1) ** slope
    dR = np.array([R_COEF[a] - R_COEF[b] for a, b in COLORS])
    qcol = rng.normal(0, 0.15, (nbg, len(COLORS))) + ebv_q[:, None] * dR
    qval = np.ones_like(qcol, dtype=bool)

    edges, rp_mid = rp_bins()
    out = os.path.expanduser("~/data/igmdust/paircount")
    rp, delta, sig, ebv, err = run_arrays(
        "selftest", (fra, fdec, fz, np.ones(nfg)),
        (rra, rdec, rz, np.ones(nrn)),
        (qra, qdec, qz, qcol, qval), edges, rp_mid, out, tag="")

    truth = amp * (rp / 0.1) ** slope
    snr = truth / err
    test = snr > 5
    frac = np.abs(ebv[test] - truth[test]) / truth[test]
    print("\n== selftest comparison (bins with expected S/N>5):")
    for r, t, e, s, f in zip(rp[test] * 1e3, truth[test], ebv[test],
                             err[test], frac):
        print(f"  rp={r:7.0f} kpc  injected={t:.3e}  recovered={e:.3e}"
              f" +-{s:.1e}  frac.dev={f:.2%}")
    RP_MAX = rp_max_save
    if np.all(frac < 0.15):
        print("SELFTEST PASS")
        return 0
    print("SELFTEST FAIL (deviation > 15%)")
    return 1


# --------------------------------------------------------------------- cli ---
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fg", choices=["bgs", "lrg"], help="foreground sample")
    ap.add_argument("--base", default=BASE_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--max-fg", type=int, default=0,
                    help="subsample foregrounds (quick look)")
    ap.add_argument("--nran", type=int, default=2, help="random files per cap")
    ap.add_argument("--qso-cat", default=None,
                    help="merged background (06 output) instead of DESI-only;"
                         " output tag gains _mq")
    ap.add_argument("--chunk", type=int, default=CHUNK,
                    help="fg per tree-query chunk (lower on 20G SLURM "
                         "with the denser merged background)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    CHUNK = a.chunk
    if a.selftest:
        sys.exit(selftest())
    if not a.fg:
        ap.error("--fg required (or --selftest)")
    run(a.fg, a.base, a.out, a.max_fg, a.nran,
        tag=("_mq" if a.qso_cat else "")
            + (f"_sub{a.max_fg}" if a.max_fg else ""),
        qso_cat=a.qso_cat)
