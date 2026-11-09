#!/usr/bin/env python3
"""
11_sweepgal.py — igmdust: photometric foreground galaxies from DR11 sweeps.

WHY (power budget, 2026-11-08): the per-pair QSO color scatter measured in
paircount_bgs_mq is 0.33-0.93 mag (intrinsic quasar SED diversity), ~8x the
feasibility assumption, so imaging depth buys nothing per pair -- power is
pair counts.  Inside 500 kpc we are shot-limited and BGS (~370/deg^2) is
~17x sparser than Menard+10's 24M photometric galaxies (~6,300/deg^2).
This script builds the Menard-replica foreground: every r<21 galaxy in the
DR11 sweeps (Menard used 17<i<21), ~17x the BGS pair counts.

    --mode galaxies   scan the 2,126 sweeps -> photgal_dr11.fits (RA,DEC,RMAG)
    --mode randoms    cut downloaded DR11 imaging randoms with the SAME
                      footprint cuts -> photran_dr11.fits (RA,DEC)

Selection (both modes where applicable): NOBS_G/R/Z >= 1 and
MASKBITS clean (BRIGHT|ALLMASK_GRZ|MEDIUM|GALAXY|CLUSTER); galaxies
additionally TYPE not in {PSF, DUP} and 16 < r < 21 (MW-dereddened).

WHERE TO RUN: astrofs4 SLURM (see 13_slurm_menard.sbatch) for galaxies;
randoms mode anywhere the download sits (cheap).
"""
import argparse
import glob
import os
import time
from multiprocessing import Pool

import numpy as np
from astropy.io import fits

LS_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11"
EXT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
MASK = (1 << 1) | (1 << 5) | (1 << 6) | (1 << 7) \
    | (1 << 11) | (1 << 12) | (1 << 13)
RMAG_RANGE = (16.0, 21.0)


def one_sweep(fname):
    """(ra, dec, rmag) of selected galaxies in one sweep, or None."""
    try:
        with fits.open(fname, memmap=True) as h:
            d = h[1].data
            typ = np.char.strip(d["TYPE"].astype(str))
            sel = (typ != "PSF") & (typ != "DUP")
            for b in ("G", "R", "Z"):
                sel &= np.asarray(d[f"NOBS_{b}"]) >= 1
            sel &= (np.asarray(d["MASKBITS"]) & MASK) == 0
            flux = np.asarray(d["FLUX_R"], float)
            mwt = np.asarray(d["MW_TRANSMISSION_R"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = flux / mwt
            ok = np.isfinite(flux) & (flux > 0)
            rmag = np.full(len(flux), np.nan)
            rmag[ok] = 22.5 - 2.5 * np.log10(flux[ok])
            sel &= ok & (rmag > RMAG_RANGE[0]) & (rmag < RMAG_RANGE[1])
            return (np.asarray(d["RA"], float)[sel],
                    np.asarray(d["DEC"], float)[sel],
                    rmag[sel].astype(np.float32))
    except OSError as e:
        print(f"  WARN: {os.path.basename(fname)}: {e}", flush=True)
        return None


def run_galaxies(a):
    sweeps = sorted(glob.glob(f"{a.ls}/south/sweep/11.0/sweep-*.fits")
                    + glob.glob(f"{a.ls}/north/sweep/11.0/sweep-*.fits")
                    + glob.glob(f"{a.ls}/south/sweep-*.fits")
                    + glob.glob(f"{a.ls}/north/sweep-*.fits"))
    print(f"  scanning {len(sweeps)} sweeps for "
          f"{RMAG_RANGE[0]} < r < {RMAG_RANGE[1]} galaxies")
    parts, tot, t0 = [], 0, time.time()
    with Pool(a.nproc) as pool:
        for i, out in enumerate(pool.imap_unordered(one_sweep, sweeps,
                                                    chunksize=2)):
            if out is not None and len(out[0]):
                parts.append(out)
                tot += len(out[0])
            if (i + 1) % 50 == 0:
                print(f"  [{i+1}/{len(sweeps)}] matched so far {tot:,}  "
                      f"{time.time()-t0:.0f}s", flush=True)
    ra = np.concatenate([p[0] for p in parts])
    dec = np.concatenate([p[1] for p in parts])
    rmag = np.concatenate([p[2] for p in parts])
    cols = [fits.Column(name="RA", format="D", array=ra),
            fits.Column(name="DEC", format="D", array=dec),
            fits.Column(name="RMAG", format="E", array=rmag)]
    fits.BinTableHDU.from_columns(cols).writeto(a.out, overwrite=True)
    print(f"  matched {len(ra):,} galaxies")
    print(f"saved -> {a.out}")


def run_randoms(a):
    files = sorted(glob.glob(a.randoms))
    if not files:
        raise SystemExit(f"no randoms match {a.randoms}")
    ras, decs, tot, t0 = [], [], 0, time.time()
    for fi, f in enumerate(files):
        with fits.open(f, memmap=True) as h:
            n = h[1].header["NAXIS2"]
            step = 20_000_000
            for j0 in range(0, n, step):
                d = h[1].data[j0:j0 + step]
                sel = (np.asarray(d["MASKBITS"]) & MASK) == 0
                for b in ("G", "R", "Z"):
                    sel &= np.asarray(d[f"NOBS_{b}"]) >= 1
                ras.append(np.asarray(d["RA"], float)[sel])
                decs.append(np.asarray(d["DEC"], float)[sel])
                tot += int(sel.sum())
                print(f"  [{fi+1}/{len(files)}] matched so far {tot:,}  "
                      f"{time.time()-t0:.0f}s", flush=True)
    ra, dec = np.concatenate(ras), np.concatenate(decs)
    if a.max_rows and len(ra) > a.max_rows:
        rng = np.random.default_rng(2)
        pick = rng.choice(len(ra), a.max_rows, replace=False)
        ra, dec = ra[pick], dec[pick]
        print(f"  thinned uniformly to {a.max_rows:,}")
    cols = [fits.Column(name="RA", format="D", array=ra),
            fits.Column(name="DEC", format="D", array=dec)]
    fits.BinTableHDU.from_columns(cols).writeto(a.out, overwrite=True)
    print(f"  kept {len(ra):,} randoms")
    print(f"saved -> {a.out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--mode", choices=["galaxies", "randoms"],
                    default="galaxies")
    ap.add_argument("--ls", default=LS_DEFAULT)
    ap.add_argument("--randoms", default=f"{LS_DEFAULT}/randoms/"
                                         "randoms-1-*.fits",
                    help="glob of downloaded DR11 imaging randoms")
    ap.add_argument("--out", default=None)
    ap.add_argument("--nproc", type=int, default=10)
    ap.add_argument("--max-rows", type=int, default=150_000_000,
                    help="uniform thinning cap for randoms (0 = keep all)")
    a = ap.parse_args()
    if a.out is None:
        a.out = f"{EXT}/photgal_dr11.fits" if a.mode == "galaxies" \
            else f"{EXT}/photran_dr11.fits"
    (run_galaxies if a.mode == "galaxies" else run_randoms)(a)
