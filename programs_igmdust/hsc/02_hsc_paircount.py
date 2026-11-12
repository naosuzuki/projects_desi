#!/usr/bin/env python3
"""
02_hsc_paircount.py — igmdust/HSC: the Legacy recipe on HSC DR4 Wide.

DESI BGS spectroscopic foregrounds (+ LSS randoms, both cut to the HSC
footprint) x HSC quiescent "crayon" background galaxies (CAS jobs
669956/669957: i<26, extended, Mizuki photoz 0.4-1.7), or the HSC
point-source star null (job 669937).  Inherits the three validated DES
ingredients: CSFD re-correction, wide radii (10 kpc - 45 Mpc), and the
saved jackknife accumulations for the covariance-ML re-solve.

Reuses projects_desi 02_paircount machinery with HSC extinction
coefficients (A_b/E(B-V), SSP standard a_* convention) and grizy
colors {g-r, r-i, i-z, z-y}; finer jackknife (nside 32) for the
~1,100 deg^2 footprint.

WHERE TO RUN: astro4 (catalogs on /home/suzuki/data/HSC, shared NFS):
  ( nohup python3 02_hsc_paircount.py --shard 0 --nshard 4 ... & )
NOTE: HSC DR4 is collaboration-internal; confirm publication rights
before HSC numbers go in a paper.
"""
import argparse
import importlib.util
import os
import re

import numpy as np
from astropy.io import fits

_here = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pc = _load("paircount02", os.path.expanduser(
    "~/github/projects_desi/programs_igmdust/02_paircount.py"))

# default: the download location (on the astrofs2_1-served home volume).
# On astrofs4 pass --hsc-dir pointing at a stage on a DIFFERENT disk
# server (04_slurm_hsc_astrofs4.sbatch rsyncs it there first), so HSC
# reads do not contend with the sweep scan hammering astrofs2_1.
HSC_DEFAULT = "/home/suzuki/data/HSC/catalog"
GAL_NAMES = ["hsc_dr4_igmdust_669956", "hsc_dr4_igmdust_669957"]
STAR_NAME = "hsc_dr4_igmdust_669937"

# HSC SSP extinction coefficients A_b/E(B-V) (the a_* column convention)
R_HSC = {"g": 3.240, "r": 2.276, "i": 1.633, "z": 1.263, "y": 1.075}
COLORS_HSC = [("g", "r"), ("r", "i"), ("i", "z"), ("z", "y")]
BANDS = ["g", "r", "i", "z", "y"]


def _rawtable(f):
    """Memory-mapped raw view of a BINTABLE's rows (bypasses FITS_rec).

    astropy materializes ~2.4x the whole table in RAM on ANY slice of
    h[1].data (measured: 43 GB for the 18 GB crayon file -- the
    2026-11-08 astrofs4 OOM), so map the packed big-endian records
    directly and let the chunk loop read ~0.5 GB at a time.
    """
    tmap = {"L": "u1", "B": "u1", "I": ">i2", "J": ">i4", "K": ">i8",
            "E": ">f4", "D": ">f8"}
    with fits.open(f, memmap=True) as h:
        hd = h[1].header.copy()
        off = h[1].fileinfo()["datLoc"]
    if hd.get("PCOUNT", 0):
        raise ValueError(f"{f}: heap/VLA tables unsupported")
    names, fmts, scaled = [], [], set()
    for i in range(1, hd["TFIELDS"] + 1):
        tform = str(hd[f"TFORM{i}"]).strip()
        m = re.fullmatch(r"(\d*)([LABIJKED])", tform)
        if not m:
            raise ValueError(f"{f}: TFORM{i}={tform} unsupported")
        rep, code = int(m.group(1) or "1"), m.group(2)
        name = str(hd.get(f"TTYPE{i}", f"col{i}")).strip()
        if code == "A":
            fmt = f"S{rep}"
        elif rep == 1:
            fmt = tmap[code]
        else:
            fmt = (tmap[code], (rep,))
        if hd.get(f"TSCAL{i}", 1) != 1 or hd.get(f"TZERO{i}", 0) != 0:
            scaled.add(name)
        names.append(name)
        fmts.append(fmt)
    dt = np.dtype({"names": names, "formats": fmts})
    if dt.itemsize != hd["NAXIS1"]:
        raise ValueError(f"{f}: parsed row {dt.itemsize} B != NAXIS1 "
                         f"{hd['NAXIS1']}")
    mm = np.memmap(f, dtype=dt, mode="r", offset=off,
                   shape=(hd["NAXIS2"],))
    return mm, scaled


def load_hsc(files, imax, zrange=None, risk_max=None, quiescent=None,
             star_z=None, chunk=4_000_000):
    """RA/DEC/z/colors/validity from the CAS fits, chunked."""
    out = {k: [] for k in ("ra", "dec", "z", "col", "val")}
    kept = tot = 0
    for f in files:
        mm, scaled = _rawtable(f)
        names = set(mm.dtype.names)
        used = ({"ra", "dec", "i_cmodel_mag"}
                | {f"{b}_cmodel_mag" for b in BANDS}
                | {f"{b}_cmodel_magerr" for b in BANDS}
                | {f"a_{b}" for b in BANDS})
        if star_z is None:
            used.add("photoz_best")
        bad = used & scaled
        if bad:
            raise ValueError(f"{f}: TSCAL/TZERO columns {sorted(bad)}")
        for j0 in range(0, len(mm), chunk):
            d = np.array(mm[j0:j0 + chunk])      # one sequential read
            tot += len(d)
            imag = np.asarray(d["i_cmodel_mag"], float)
            sel = np.isfinite(imag) & (imag < imax) & (imag > 15)
            if zrange is not None:
                z = np.asarray(d["photoz_best"], float)
                sel &= (z >= zrange[0]) & (z <= zrange[1])
            if risk_max is not None and "photoz_risk_best" in names:
                r = np.asarray(d["photoz_risk_best"], float)
                sel &= np.isfinite(r) & (r < risk_max)
            if quiescent is not None and {"sfr", "stellar_mass"} <= names:
                sfr = np.asarray(d["sfr"], float)
                sm = np.asarray(d["stellar_mass"], float)
                ok = np.isfinite(sfr) & np.isfinite(sm) & (sm > 0) \
                    & (sfr >= 0)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ssfr = np.where(ok, sfr / np.maximum(sm, 1e-30),
                                    np.inf)
                sel &= ok & (ssfr < 10.0 ** quiescent)
            mags, errs = {}, {}
            for b in BANDS:
                m = np.asarray(d[f"{b}_cmodel_mag"], float)
                a = np.asarray(d[f"a_{b}"], float)
                mags[b] = m - a                     # SFD-dereddened
                errs[b] = np.asarray(d[f"{b}_cmodel_magerr"], float)
            col = np.column_stack([mags[x] - mags[y]
                                   for x, y in COLORS_HSC])
            err_ok = np.column_stack(
                [np.isfinite(errs[x]) & (errs[x] < 0.5)
                 & np.isfinite(errs[y]) & (errs[y] < 0.5)
                 for x, y in COLORS_HSC])
            val = np.isfinite(col) & err_ok
            out["ra"].append(np.asarray(d["ra"], float)[sel])
            out["dec"].append(np.asarray(d["dec"], float)[sel])
            zc = (np.full(int(sel.sum()), star_z) if star_z is not None
                  else np.asarray(d["photoz_best"], float)[sel])
            out["z"].append(zc)
            out["col"].append(np.where(val, col, 0.0)[sel])
            out["val"].append(val[sel])
            kept += int(sel.sum())
            print(f"  {os.path.basename(f)}: {kept:,}/{tot:,} kept",
                  flush=True)
        del mm
    ra = np.concatenate(out["ra"])
    dec = np.concatenate(out["dec"])
    z = np.concatenate(out["z"])
    col = np.vstack(out["col"])
    val = np.vstack(out["val"])
    print(f"  HSC backgrounds: {len(ra):,} "
          f"(validity {val.mean(0).round(3)})")
    return ra, dec, z, col, val


def _pix(ra, dec, nside):
    """Ring-less nested-free pixelization: RA/sin(dec) grid at ~nside
    resolution (no healpy on astro4; pure numpy, same spirit as 02)."""
    step = 58.6 / nside                                # ~pixel scale deg
    i = np.floor(ra / step).astype(np.int64)
    j = np.floor((np.sin(np.radians(dec)) + 1.0) / (step / 57.3)
                 ).astype(np.int64)
    return i * 100_000 + j


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fg", choices=["bgs", "lrg"], default="bgs")
    ap.add_argument("--mode", choices=["crayons", "stars"],
                    default="crayons")
    ap.add_argument("--imax", type=float, default=24.0)
    ap.add_argument("--zrange", type=float, nargs=2, default=[0.7, 1.7])
    ap.add_argument("--risk-max", type=float, default=0.5)
    ap.add_argument("--quiescent", type=float, default=-10.7,
                    help="log10 sSFR/yr cut (Mizuki); 99 disables")
    ap.add_argument("--rp-min", type=float, default=0.01)
    ap.add_argument("--rp-max", type=float, default=45.0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=64,
                    help="fg per tree query; HSC bg density (~7k-20k "
                         "per deg2) x 45 Mpc aperture means ~1e5 pairs "
                         "per fg -- 64 keeps query buffers ~2 GB")
    ap.add_argument("--hsc-dir", default=HSC_DEFAULT,
                    help="directory holding the three CAS fits (point at "
                         "a stage on another disk server to split I/O)")
    ap.add_argument("--base",
                    default="/data/astrofs2_1/suzuki/data/DESI/catalog/dr1")
    ap.add_argument("--out", default=os.path.expanduser(
        "~/data/igmdust/paircount"))
    a = ap.parse_args()

    # HSC configuration of the shared machinery
    pc.R_COEF = R_HSC
    pc.COLORS = COLORS_HSC
    pc.JK_NSIDE = 32
    pc.JK_GRID_DEG = 1.8
    for k in ("hscbgs", "hscbgsstars", "hsclrg", "hsclrgstars"):
        pc.FG_ZRANGE[k] = pc.FG_ZRANGE[a.fg]
    pc.RP_MIN, pc.RP_MAX = a.rp_min, a.rp_max
    pc.CHUNK = a.chunk
    edges, rp_mid = pc.rp_bins()

    print(f"== HSC backgrounds: {a.mode} (from {a.hsc_dir})")
    if a.mode == "crayons":
        q = None if a.quiescent > 90 else a.quiescent
        bra, bdec, bz, bcol, bval = load_hsc(
            [f"{a.hsc_dir}/{n}" for n in GAL_NAMES], a.imax,
            zrange=a.zrange, risk_max=a.risk_max, quiescent=q)
        name = f"hsc{a.fg}"
    else:
        bra, bdec, bz, bcol, bval = load_hsc(
            [f"{a.hsc_dir}/{STAR_NAME}"], 22.5, star_z=99.0)
        name = f"hsc{a.fg}stars"

    # CSFD re-correction with HSC coefficients (DES ingredient 1)
    c6 = _load("csfd06", os.path.expanduser(
        "~/github/projects_desi/programs_igmdust/des/06_csfd.py"))
    dE = c6.delta_ebv(bra, bdec)
    dRc = np.array([R_HSC[x] - R_HSC[y] for x, y in COLORS_HSC])
    bcol = np.where(bval, bcol + dE[:, None] * dRc[None, :], 0.0)
    print(f"  CSFD recorrection (HSC R): mean {dE.mean():+.5f}, "
          f"rms {dE.std():.5f} mag")

    print(f"== foregrounds: {a.fg}")
    fra, fdec, fz, fw = pc.load_lss(a.base, a.fg, "dat")
    rra, rdec, rz, rw = pc.load_lss(a.base, a.fg, "ran", nran=2)
    fin = np.isin(_pix(fra, fdec, 64), np.unique(_pix(bra, bdec, 64)))
    rin = np.isin(_pix(rra, rdec, 64), np.unique(_pix(bra, bdec, 64)))
    print(f"  HSC footprint cut: fg {fin.sum():,}/{len(fra):,}  "
          f"ran {rin.sum():,}/{len(rra):,}")
    fg = (fra[fin], fdec[fin], fz[fin], fw[fin])
    ran = (rra[rin], rdec[rin], rz[rin], rw[rin])

    tag = "csfd" + (f"_s{a.nshard}.{a.shard}" if a.nshard > 1 else "")
    pc.run_arrays(name, fg, ran, (bra, bdec, bz, bcol, bval),
                  edges, rp_mid, a.out, tag=tag,
                  nshard=a.nshard, shard=a.shard)


if __name__ == "__main__":
    main()
