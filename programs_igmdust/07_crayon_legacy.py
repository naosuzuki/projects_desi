#!/usr/bin/env python3
"""
07_crayon_legacy.py — igmdust/Legacy: the standard-crayon recipe with
DESI spectroscopic LRGs as crayons on Legacy Surveys DR11 photometry.

The McCleary-style measurement, upgraded by spectroscopy: the crayons
are DR1 LRG LSS galaxies (0.4 < z < 1.1) whose red-sequence locus is
anchored on their SPECTROSCOPIC redshifts, so the photo-z reddening
absorption of the DES version is exactly zero (f_abs = 1) and no
correction factor is needed.  Photometry: DR11 sweep FLUX_{G,R,Z,W1}
attached by 06_sweepmatch (run with --in on the LRG positions file that
--make-input writes).  Colors {g-r, r-z, g-z, z-W1}; foregrounds are
the DR1 BGS galaxies and their LSS randoms; radii 10 kpc - 45 Mpc.

WHERE TO RUN: astrofs4/SLURM.
    python3 07_crayon_legacy.py --make-input <out.fits>   # positions
    python3 07_crayon_legacy.py                           # measurement
    python3 07_crayon_legacy.py --selftest
Output: $OUT/paircount_bgslrgcrayon_qso.npz
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "paircount02", os.path.join(_here, "02_paircount.py"))
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)

EXT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
CRAYON_ZRANGE = (0.40, 1.10)
LOCUS_DZ = 0.01
LOCUS_NMIN = 200
MAGERR_IVAR_MIN = 0.0                # ivar>0 is the validity requirement


def make_input(base, out):
    ra, dec, z, w = pc.load_lss(base, "lrg", "dat")
    Table(dict(RA=ra, DEC=dec, Z=z, WEIGHT=w)).write(out, overwrite=True)
    print(f"  LRG positions -> {out} ({len(ra):,})")


def locus_residuals(z, colors, valid):
    """Color minus the spec-z red-sequence locus (median per dz bin)."""
    zb = np.arange(CRAYON_ZRANGE[0], CRAYON_ZRANGE[1] + LOCUS_DZ, LOCUS_DZ)
    zc = 0.5 * (zb[:-1] + zb[1:])
    res = np.full_like(colors, np.nan)
    for c in range(colors.shape[1]):
        med = np.full(len(zc), np.nan)
        for i in range(len(zc)):
            m = valid[:, c] & (z >= zb[i]) & (z < zb[i + 1])
            if m.sum() >= LOCUS_NMIN:
                med[i, ] = np.median(colors[m, c])
        ok = np.isfinite(med)
        res[:, c] = colors[:, c] - np.interp(z, zc[ok], med[ok])
        r = res[valid[:, c], c]
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        a, b = pc.COLORS[c]
        print(f"  locus {a}-{b}: {ok.sum()} z-bins, scatter {sig:.3f} mag")
    return res


def load_crayons(path):
    """LRG crayons with DR11 photometry (06_sweepmatch output)."""
    with fits.open(path) as h:
        d = h[1].data
        ra = np.asarray(d["RA"], float)
        dec = np.asarray(d["DEC"], float)
        z = np.asarray(d["Z"], float)
        sel = (np.asarray(d["MATCHED"], bool) & np.isfinite(z)
               & (z >= CRAYON_ZRANGE[0]) & (z <= CRAYON_ZRANGE[1]))
        mag, ok = {}, {}
        for b in ("G", "R", "Z", "W1"):
            flux = np.asarray(d[f"FLUX_{b}"], float)
            ivar = np.asarray(d[f"IVAR_{b}"], float)
            mwt = np.asarray(d[f"MWT_{b}"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = flux / mwt
            k = b if b == "W1" else b.lower()
            ok[k] = (np.isfinite(flux) & (flux > 0)
                     & np.isfinite(ivar) & (ivar > MAGERR_IVAR_MIN))
            m = np.full(len(flux), np.nan)
            m[ok[k]] = 22.5 - 2.5 * np.log10(flux[ok[k]])
            mag[k] = m
    colors = np.column_stack([mag[a] - mag[b] for a, b in pc.COLORS])
    valid = np.column_stack([ok[a] & ok[b] for a, b in pc.COLORS])
    print(f"  LRG crayons: {sel.sum():,} matched in z range; "
          f"validity {valid[sel].mean(0).round(3)}")
    res = locus_residuals(z[sel], colors[sel], valid[sel])
    v = valid[sel] & np.isfinite(res)
    return ra[sel], dec[sel], z[sel], np.where(v, res, 0.0), v


def main(a):
    name = "bgslrgcrayon"
    pc.FG_ZRANGE[name] = pc.FG_ZRANGE["bgs"]
    pc.RP_MIN, pc.RP_MAX = a.rp_min, a.rp_max
    pc.CHUNK = 3_000
    edges, rp_mid = pc.rp_bins()

    print("== crayons: DESI DR1 LRG spec-z x DR11 photometry")
    bg = load_crayons(a.crayons)
    print("== foregrounds: DR1 BGS + LSS randoms")
    fra, fdec, fz, fw = pc.load_lss(a.base, "bgs", "dat")
    rra, rdec, rz, rw = pc.load_lss(a.base, "bgs", "ran", nran=a.nran)

    tag = f"_sub{a.max_fg}" if a.max_fg else ""
    if a.nshard > 1:
        tag += f"_s{a.nshard}.{a.shard}"
    pc.run_arrays(name, (fra, fdec, fz, fw), (rra, rdec, rz, rw), bg,
                  edges, rp_mid, a.out, a.max_fg, tag=tag,
                  nshard=a.nshard, shard=a.shard)


def selftest():
    """Inject E(B-V)=0.05 on half of synthetic spec-z crayons; recover."""
    rng = np.random.default_rng(13)
    n = 200_000
    z = rng.uniform(*CRAYON_ZRANGE, n)
    dR = np.array([pc.R_COEF[a] - pc.R_COEF[b] for a, b in pc.COLORS])
    locus = np.column_stack([0.9 + 1.1 * z, 0.8 + 0.9 * z,
                             1.7 + 2.0 * z, 1.1 + 0.6 * z])
    ebv = np.where(rng.random(n) < 0.5, 0.05, 0.0)
    colors = locus + rng.normal(0, 0.06, (n, 4)) + ebv[:, None] * dR
    valid = np.ones((n, 4), bool)
    res = locus_residuals(z, colors, valid)
    d = res[ebv > 0].mean(0) - res[ebv == 0].mean(0)
    frac = np.abs(d - dR * 0.05) / (dR * 0.05)
    for i, (a, b) in enumerate(pc.COLORS):
        print(f"  {a}-{b}: recovered {d[i]:+.5f}  "
              f"injected {0.05 * dR[i]:+.5f}  frac.dev {frac[i]:.2%}")
    print("SELFTEST PASS" if np.all(frac < 0.05) else "SELFTEST FAIL")
    return 0 if np.all(frac < 0.05) else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--base", default=pc.BASE_DEFAULT)
    ap.add_argument("--out", default=pc.OUT_DEFAULT)
    ap.add_argument("--crayons", default=f"{EXT}/lrg_dr11.fits")
    ap.add_argument("--make-input", default=None,
                    help="write LRG positions fits for 06_sweepmatch")
    ap.add_argument("--max-fg", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1,
                    help="split fg+randoms into N interleaved shards "
                         "(merge with 02b_merge_shards --name bgslrgcrayon)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nran", type=int, default=2)
    ap.add_argument("--rp-min", type=float, default=0.01)
    ap.add_argument("--rp-max", type=float, default=45.0)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if a.make_input:
        make_input(a.base, a.make_input)
        sys.exit(0)
    main(a)
