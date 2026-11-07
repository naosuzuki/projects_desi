#!/usr/bin/env python3
"""
05_qso_merge.py — igmdust: the merged QSO background sample.

Combines three quasar catalogs into one background list (positions and
redshifts only; uniform Legacy Surveys DR11 photometry is attached
afterwards by 06_sweepmatch.py, so survey-specific photometry never
enters the estimator):

  SRC 0  DESI DR1 iron     spec-z, main-dark, BAL-cleaned (as 02)
  SRC 1  SDSS DR16Q v4     spec-z, IS_QSO_FINAL / ZWARNING / BAL_PROB cuts
  SRC 2  Quaia G<20.5      Gaia-unWISE spectrophotometric z, zerr cut

Duplicates are removed by a 1" positional match with priority
DESI > SDSS > Quaia (spec-z wins; among spec-z, DESI wins).
Common cut: 0.3 < z < 4.0 (matching 02_paircount).

WHERE TO RUN: astro4 (reads /data/astrofs2_1; ~2 min).
    python3 05_qso_merge.py
    python3 05_qso_merge.py --selftest
Output: $EXT/qso_merged.fits  (RA, DEC, Z, ZERR, SRC)
"""
import argparse
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.spatial import cKDTree

BASE_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/dr1"
EXT_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
MATCH_ARCSEC = 1.0
ZRANGE = (0.3, 4.0)
QUAIA_ZERR_MAX = 0.15
SRC_NAMES = ["DESI", "SDSS", "Quaia"]


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    c = np.cos(dec)
    return np.column_stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)])


def getcol(data, names, *cands):
    for c in cands:
        if c.upper() in names:
            return np.asarray(data[c])
    return None


def load_desi(base):
    """Same selection as 02_paircount.load_qso, positions+z only."""
    with fits.open(f"{base}/QSO_cat_iron_cumulative_v0.fits") as h:
        d = h[1].data
        names = {c.upper() for c in h[1].columns.names}
        ra = getcol(d, names, "TARGET_RA", "RA")
        dec = getcol(d, names, "TARGET_DEC", "DEC")
        z = np.asarray(d["Z"], float)
        sel = np.isfinite(z) & (z > ZRANGE[0]) & (z < ZRANGE[1])
        if "SURVEY" in names:
            sel &= np.char.strip(d["SURVEY"].astype(str)) == "main"
        if "PROGRAM" in names:
            sel &= np.char.strip(d["PROGRAM"].astype(str)) == "dark"
        if "BI_CIV" in names:
            bi = np.asarray(d["BI_CIV"], float)
            sel &= ~(np.isfinite(bi) & (bi > 0))
    print(f"  DESI: {sel.sum():,} / {len(z):,}")
    return ra[sel], dec[sel], z[sel], np.zeros(sel.sum())


def load_dr16q(ext):
    with fits.open(f"{ext}/DR16Q_v4.fits") as h:
        d = h[1].data
        names = {c.upper() for c in h[1].columns.names}
        ra = getcol(d, names, "RA")
        dec = getcol(d, names, "DEC")
        z = np.asarray(getcol(d, names, "Z"), float)
        sel = np.isfinite(z) & (z > ZRANGE[0]) & (z < ZRANGE[1])
        for col, cond in (("IS_QSO_FINAL", lambda v: v > 0),
                          ("ZWARNING", lambda v: v == 0)):
            v = getcol(d, names, col)
            if v is not None:
                sel &= cond(np.asarray(v))
        bp = getcol(d, names, "BAL_PROB")
        if bp is not None:
            sel &= ~(np.isfinite(np.asarray(bp, float))
                     & (np.asarray(bp, float) >= 0.5))
    print(f"  SDSS DR16Q: {sel.sum():,} / {len(z):,}")
    return ra[sel], dec[sel], z[sel], np.zeros(sel.sum())


def load_quaia(ext):
    with fits.open(f"{ext}/quaia_G20.5.fits") as h:
        d = h[1].data
        names = {c.upper() for c in h[1].columns.names}
        ra = getcol(d, names, "RA")
        dec = getcol(d, names, "DEC")
        z = np.asarray(getcol(d, names, "REDSHIFT_QUAIA", "Z"), float)
        ze = getcol(d, names, "REDSHIFT_QUAIA_ERR", "ZERR")
        ze = (np.asarray(ze, float) if ze is not None
              else np.zeros(len(z)))
        sel = (np.isfinite(z) & (z > ZRANGE[0]) & (z < ZRANGE[1])
               & np.isfinite(ze) & (ze < QUAIA_ZERR_MAX))
    print(f"  Quaia (zerr<{QUAIA_ZERR_MAX}): {sel.sum():,} / {len(z):,}")
    return ra[sel], dec[sel], z[sel], ze[sel]


def dedup(samples):
    """Keep higher-priority sources; drop later ones within MATCH_ARCSEC."""
    chord = 2 * np.sin(np.radians(MATCH_ARCSEC / 3600) / 2)
    ra, dec, z, ze, src = (np.empty(0),) * 4 + (np.empty(0, int),)
    for isrc, (r, d, zz, zerr) in enumerate(samples):
        if len(ra):
            tree = cKDTree(unitvec(ra, dec))
            dist, _ = tree.query(unitvec(r, d), k=1,
                                 distance_upper_bound=chord)
            new = ~np.isfinite(dist)
        else:
            new = np.ones(len(r), bool)
        print(f"  {SRC_NAMES[isrc]}: +{new.sum():,} new "
              f"({(~new).sum():,} duplicates dropped)")
        ra = np.concatenate([ra, r[new]])
        dec = np.concatenate([dec, d[new]])
        z = np.concatenate([z, zz[new]])
        ze = np.concatenate([ze, zerr[new]])
        src = np.concatenate([src, np.full(new.sum(), isrc)])
    return ra, dec, z, ze, src


def run(a):
    print("== loading")
    samples = [load_desi(a.base), load_dr16q(a.ext), load_quaia(a.ext)]
    print("== dedup (1\", DESI > SDSS > Quaia)")
    ra, dec, z, ze, src = dedup(samples)
    t = Table(dict(RA=ra, DEC=dec, Z=z, ZERR=ze, SRC=src.astype(np.int16)))
    t.meta["HISTORY"] = f"05_qso_merge.py {time.strftime('%F %T')}"
    fo = f"{a.ext}/qso_merged.fits"
    t.write(fo, overwrite=True)
    print(f"  total {len(t):,} -> {fo}")


def selftest():
    """Priority dedup on a synthetic overlap: SDSS/Quaia clones of DESI
    positions must be dropped; offset ones kept."""
    rng = np.random.default_rng(5)
    n = 2000
    ra = rng.uniform(100, 110, n)
    dec = rng.uniform(-5, 5, n)
    z = rng.uniform(0.5, 3.0, n)
    zero = np.zeros(n)
    # SDSS: first half = clones shifted 0.3" (dup), second half shifted 10"
    off_dup = 0.3 / 3600 / np.cos(np.radians(dec[:1000]))
    off_new = 10.0 / 3600 / np.cos(np.radians(dec[1000:]))
    sdss = (np.concatenate([ra[:1000] + off_dup, ra[1000:] + off_new]),
            dec.copy(), z.copy(), zero.copy())
    ra2, dec2, z2, ze2, src = dedup([(ra, dec, z, zero), sdss])
    n_sdss = (src == 1).sum()
    print(f"  kept from SDSS: {n_sdss} (expect 1000)")
    if n_sdss == 1000 and len(ra2) == 3000:
        print("SELFTEST PASS")
        return 0
    print("SELFTEST FAIL")
    return 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--base", default=BASE_DEFAULT)
    ap.add_argument("--ext", default=EXT_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    sys.exit(selftest() if a.selftest else run(a))
