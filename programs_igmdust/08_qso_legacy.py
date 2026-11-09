#!/usr/bin/env python3
"""
08_qso_legacy.py — igmdust/Legacy: the quasar-recipe star null test.

Stars carry no cosmological reddening, so running the Menard estimator
with DR11 PSF stars in place of quasars yields a pure imaging-systematics
template: any nonzero "profile" around BGS foregrounds measures
data-vs-randoms photometric systematics (the designed arbiter of the
constant per-color offsets seen in the first BGS x QSO measurement).

Two modes:
  --make-stars   scan the DR11 sweeps (parallel), select TYPE=PSF stars
                 with 17 < r < 21.5, thin to ~STARS_PER_SWEEP per sweep,
                 write $EXT/stars_dr11.fits (grzW1 fluxes + MW transm.)
  (default)      run the estimator: BGS foregrounds + LSS randoms vs the
                 star 'background' (z set far beyond every foreground so
                 the redshift buffer never rejects a pair)

WHERE TO RUN: astrofs4/SLURM (make-stars is I/O bound like 06).
Output: $EXT/stars_dr11.fits, $OUT/paircount_bgsstars_qso.npz
"""
import argparse
import glob
import importlib.util
import os
import sys
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

_here = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_here, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pc = _load("paircount02", "02_paircount.py")
sw = _load("sweepmatch06", "06_sweepmatch.py")

EXT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
RMAG_RANGE = (17.0, 21.5)
STARS_PER_SWEEP = 2000               # ~4M stars over 2,126 sweeps
BANDS = ("G", "R", "Z", "W1")


def stars_from_sweep(fname):
    b = sw.parse_bounds(fname)
    if b is None:
        return None
    try:
        with fits.open(fname, memmap=True) as h:
            d = h[1].data
            psf = np.char.strip(d["TYPE"].astype(str)) == "PSF"
            flux_r = np.asarray(d["FLUX_R"], float)
            mwt_r = np.asarray(d["MW_TRANSMISSION_R"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                rmag = 22.5 - 2.5 * np.log10(flux_r / mwt_r)
            keep = np.where(psf & np.isfinite(rmag)
                            & (rmag > RMAG_RANGE[0])
                            & (rmag < RMAG_RANGE[1]))[0]
            if len(keep) == 0:
                return None
            rng = np.random.default_rng(abs(hash(fname)) % 2**32)
            if len(keep) > STARS_PER_SWEEP:
                keep = rng.choice(keep, STARS_PER_SWEEP, replace=False)
            cols = dict(RA=np.asarray(d["RA"][keep], float),
                        DEC=np.asarray(d["DEC"][keep], float))
            for band in BANDS:
                cols[f"FLUX_{band}"] = np.asarray(
                    d[f"FLUX_{band}"][keep], float)
                cols[f"IVAR_{band}"] = np.asarray(
                    d[f"FLUX_IVAR_{band}"][keep], float)
                cols[f"MWT_{band}"] = np.asarray(
                    d[f"MW_TRANSMISSION_{band}"][keep], float)
            return Table(cols)
    except OSError:
        return None


def make_stars(ls, out, nproc):
    sweeps = sorted(glob.glob(f"{ls}/south/sweep-*.fits")
                    + glob.glob(f"{ls}/north/sweep-*.fits"))
    sweeps = [s for s in sweeps if sw.parse_bounds(s)]
    print(f"  scanning {len(sweeps)} sweeps for PSF stars "
          f"r in {RMAG_RANGE} ...")
    parts = []
    with Pool(nproc) as pool:
        for i, t in enumerate(pool.imap_unordered(stars_from_sweep,
                                                  sweeps, chunksize=4)):
            if t is not None:
                parts.append(t)
            if (i + 1) % 200 == 0:
                n = sum(len(p) for p in parts)
                print(f"  [{i+1}/{len(sweeps)}] stars so far {n:,}",
                      flush=True)
    t = vstack(parts)
    t.write(out, overwrite=True)
    print(f"  stars -> {out} ({len(t):,})")


def load_stars(path):
    with fits.open(path) as h:
        d = h[1].data
        ra = np.asarray(d["RA"], float)
        dec = np.asarray(d["DEC"], float)
        mag, ok = {}, {}
        for band in BANDS:
            flux = np.asarray(d[f"FLUX_{band}"], float)
            ivar = np.asarray(d[f"IVAR_{band}"], float)
            mwt = np.asarray(d[f"MWT_{band}"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = flux / mwt
            k = band if band == "W1" else band.lower()
            ok[k] = (np.isfinite(flux) & (flux > 0)
                     & np.isfinite(ivar) & (ivar > 0))
            m = np.full(len(flux), np.nan)
            m[ok[k]] = 22.5 - 2.5 * np.log10(flux[ok[k]])
            mag[k] = m
    colors = np.column_stack([mag[a] - mag[b] for a, b in pc.COLORS])
    valid = np.column_stack([ok[a] & ok[b] for a, b in pc.COLORS])
    z = np.full(len(ra), 9.9)        # never rejected by the z buffer
    print(f"  stars: {len(ra):,}; validity {valid.mean(0).round(3)}")
    return ra, dec, z, np.where(valid, np.nan_to_num(colors), 0.0), valid


def run_null(a):
    name = "bgsstars"
    pc.FG_ZRANGE[name] = pc.FG_ZRANGE["bgs"]
    pc.RP_MIN, pc.RP_MAX = a.rp_min, a.rp_max
    pc.CHUNK = 3_000
    edges, rp_mid = pc.rp_bins()
    print("== backgrounds: DR11 PSF stars (null template)")
    bg = load_stars(a.stars)
    print("== foregrounds: DR1 BGS + LSS randoms")
    fg = pc.load_lss(a.base, "bgs", "dat")
    ran = pc.load_lss(a.base, "bgs", "ran", nran=a.nran)
    pc.run_arrays(name, fg, ran, bg, edges, rp_mid, a.out,
                  a.max_fg, tag=f"_sub{a.max_fg}" if a.max_fg else "")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--make-stars", action="store_true")
    ap.add_argument("--ls", default=sw.LS_DEFAULT)
    ap.add_argument("--stars", default=f"{EXT}/stars_dr11.fits")
    ap.add_argument("--base", default=pc.BASE_DEFAULT)
    ap.add_argument("--out", default=pc.OUT_DEFAULT)
    ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--max-fg", type=int, default=0)
    ap.add_argument("--nran", type=int, default=2)
    ap.add_argument("--rp-min", type=float, default=0.01)
    ap.add_argument("--rp-max", type=float, default=45.0)
    a = ap.parse_args()
    if a.make_stars:
        make_stars(a.ls, a.stars, a.nproc)
        sys.exit(0)
    run_null(a)
