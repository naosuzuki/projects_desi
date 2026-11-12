#!/usr/bin/env python3
"""
05_qso_menard_des.py — igmdust/DES: the Menard+2010 recipe on the DES
footprint with quasars.

Menard, Scranton, Fukugita & Richards (2010) measured <delta color> of
SDSS quasars around SDSS galaxies.  This is the same recipe with DES-area
data: foregrounds are the redMaGiC high-density galaxies (0.15<z<0.45,
weights + their own randoms, as in the crayon runs), backgrounds are the
merged quasar catalog (DESI + SDSS + Quaia) with DECam g,r,i,z photometry
attached from the LS DR11 sweeps by 06_sweepmatch (same camera and
filters as DES; swapping in DES-pipeline magnitudes is a later
refinement requiring an authenticated crossmatch).  Colors
{g-r, r-i, i-z, g-i}; PSF morphology; validity needs flux>0 and ivar>0
in both bands.  Radii 10 kpc - 45 Mpc as for the crayon runs.

WHERE TO RUN: astrofs4/SLURM, AFTER 06_sweepmatch has written
qso_merged_dr11.fits:
    sbatch -J igmdust_mqd -d afterok:<job60>,singleton \
        03_slurm_paircount.sbatch menarddes
Output: $OUT/paircount_menarddes_qso.npz
"""
import argparse
import importlib.util
import os

import numpy as np
from astropy.io import fits

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "crayon01", os.path.join(_here, "01_crayon_paircount.py"))
m01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m01)

QCAT_DEFAULT = ("/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext/"
                "qso_merged_dr11.fits")
FG_ZRANGE_MD = (0.15, 0.45)


def load_qso_des(path):
    """Merged QSOs with DECam griz from the DR11 sweep match."""
    with fits.open(path) as h:
        d = h[1].data
        ra = np.asarray(d["RA"], float)
        dec = np.asarray(d["DEC"], float)
        z = np.asarray(d["Z"], float)
        src = np.asarray(d["SRC"], int)
        sel = (np.asarray(d["MATCHED"], bool)
               & (np.char.strip(d["LS_TYPE"].astype(str)) == "PSF")
               & np.isfinite(z) & (z > 0.3) & (z < 4.0))
        mag, ok = {}, {}
        for b in "GRIZ":
            flux = np.asarray(d[f"FLUX_{b}"], float)
            ivar = np.asarray(d[f"IVAR_{b}"], float)
            mwt = np.asarray(d[f"MWT_{b}"], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = flux / mwt
            ok[b.lower()] = (np.isfinite(flux) & (flux > 0)
                             & np.isfinite(ivar) & (ivar > 0))
            m = np.full(len(flux), np.nan)
            m[ok[b.lower()]] = 22.5 - 2.5 * np.log10(flux[ok[b.lower()]])
            mag[b.lower()] = m
    colors = np.column_stack([mag[a] - mag[b]
                              for a, b in m01.COLORS_DES])
    valid = np.column_stack([ok[a] & ok[b] for a, b in m01.COLORS_DES])
    colors = np.where(valid, np.nan_to_num(colors), 0.0)
    cnt = np.bincount(src[sel & valid.any(1)], minlength=3)
    print(f"  QSO bg: {(sel & valid.any(1)).sum():,} with DECam griz "
          f"(DESI {cnt[0]:,}, SDSS {cnt[1]:,}, Quaia {cnt[2]:,})")
    return ra[sel], dec[sel], z[sel], colors[sel], valid[sel]


def main(a):
    pc = m01.load_pc(a.pc02)
    pc.R_COEF["i"] = m01.R_I
    pc.COLORS = m01.COLORS_DES
    pc.JK_NSIDE = m01.JK_NSIDE
    pc.FG_ZRANGE["menarddes"] = FG_ZRANGE_MD
    pc.RP_MIN, pc.RP_MAX = a.rp_min, a.rp_max
    pc.CHUNK = 3_000
    edges, rp_mid = pc.rp_bins()

    print("== backgrounds: merged QSOs x DECam griz (DR11 sweeps)")
    bg = load_qso_des(a.qcat)
    print("== foregrounds: DES redMaGiC highdens + randoms")
    fg, ran = m01.load_redmagic_fg(a.des_dir)

    pc.run_arrays("menarddes", fg, ran, bg, edges, rp_mid, a.out,
                  a.max_fg, tag=f"_sub{a.max_fg}" if a.max_fg else "")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--qcat", default=QCAT_DEFAULT)
    ap.add_argument("--des-dir", default=m01.DES_DIR_DEFAULT)
    ap.add_argument("--pc02", default=m01.PC02_DEFAULT)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-fg", type=int, default=0)
    ap.add_argument("--rp-min", type=float, default=0.01)
    ap.add_argument("--rp-max", type=float, default=45.0)
    a = ap.parse_args()
    a.out = a.out or m01.load_pc(a.pc02).OUT_DEFAULT
    main(a)
