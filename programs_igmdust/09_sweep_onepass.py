#!/usr/bin/env python3
"""
09_sweep_onepass.py — igmdust/Legacy: ONE pass over the DR11 sweeps
serving every remaining consumer, so the 1.9 TB is read once, not once
per task (the scans are NFS-I/O limited; compute is negligible).

Per sweep (read once into memory): (a) 1" match of the LRG crayon
positions inside the sweep's block (06 logic), (b) PSF-star selection
for the null test (08 logic).  Outputs lrg_dr11.fits + stars_dr11.fits,
identical in format to what 06/08 would have produced separately.

WHERE TO RUN: astrofs4/SLURM:
    sbatch -J igmdust_1ps 03_slurm_paircount.sbatch onepass
"""
import argparse
import glob
import importlib.util
import os
import time
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
st = _load("qsolegacy08", "08_qso_legacy.py")

EXT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
_ctx = {}


def _init(ra, dec):
    _ctx["ra"], _ctx["dec"] = ra, dec
    _ctx["vec"] = sw.unitvec(ra, dec)


def one_sweep(fname):
    """Returns (match-dict or None, star Table or None)."""
    b = sw.parse_bounds(fname)
    if b is None:
        return None, None
    # (a) LRG match, reusing 06's worker against our shared positions
    sw._qso["ra"], sw._qso["dec"] = _ctx["ra"], _ctx["dec"]
    sw._qso["vec"] = _ctx["vec"]
    try:
        out = sw.match_one(fname)
    except OSError:
        out = None
    # (b) stars from the same (OS-cached) file
    stars = st.stars_from_sweep(fname)
    return out, stars


def run(a):
    q = Table.read(a.inp)
    ra = np.asarray(q["RA"], float)
    dec = np.asarray(q["DEC"], float)
    n = len(q)
    sweeps = sorted(glob.glob(f"{a.ls}/south/sweep-*.fits")
                    + glob.glob(f"{a.ls}/north/sweep-*.fits"))
    sweeps = [s for s in sweeps if sw.parse_bounds(s)]
    print(f"  one pass over {len(sweeps)} sweeps: "
          f"matching {n:,} LRGs + extracting stars")

    res = {"dist": np.full(n, np.inf), "type": np.zeros(n, "U4")}
    for band in sw.BANDS:
        for k in (f"flux_{band}", f"ivar_{band}", f"mwt_{band}"):
            res[k] = np.full(n, np.nan)
    star_parts, t0 = [], time.time()
    with Pool(a.nproc, initializer=_init, initargs=(ra, dec)) as pool:
        for i, (out, stars) in enumerate(
                pool.imap_unordered(one_sweep, sweeps, chunksize=2)):
            if out is not None:
                take = out["dist"] < res["dist"][out["idx"]]
                ix = out["idx"][take]
                for k in res:
                    if k != "dist":
                        res[k][ix] = out[k][take]
                res["dist"][ix] = out["dist"][take]
            if stars is not None:
                star_parts.append(stars)
            if (i + 1) % 100 == 0:
                print(f"  [{i+1}/{len(sweeps)}] matched so far "
                      f"{np.isfinite(res['dist']).sum():,}  "
                      f"{time.time()-t0:.0f}s", flush=True)

    matched = np.isfinite(res["dist"])
    print(f"  matched {matched.sum():,} / {n:,} ({matched.mean():.1%})")
    for k, v in res.items():
        q[k.upper() if k != "type" else "LS_TYPE"] = v
    q["MATCHED"] = matched
    q.write(a.out, overwrite=True)
    print(f"saved -> {a.out}")
    tstars = vstack(star_parts)
    tstars.write(a.stars, overwrite=True)
    print(f"  stars -> {a.stars} ({len(tstars):,})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--in", dest="inp", default=f"{EXT}/lrg_positions.fits")
    ap.add_argument("--out", default=f"{EXT}/lrg_dr11.fits")
    ap.add_argument("--stars", default=f"{EXT}/stars_dr11.fits")
    ap.add_argument("--ls", default=sw.LS_DEFAULT)
    ap.add_argument("--nproc", type=int, default=6)
    run(ap.parse_args())
