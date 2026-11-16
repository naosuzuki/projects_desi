#!/usr/bin/env python3
"""32_rmg_dr11match.py — Path A step 1: DR11 colors for her crayons.

Reads the DES Y3 redMaGiC highlum_highz crayon catalog, writes an
uppercase RA/DEC table, and hands it to 06_sweepmatch.py, which walks
the DR11 sweeps and attaches FLUX_{G,R,I,Z,W1} + ivars + MW
transmissions within 1".  The output feeds des/01_crayon_paircount.py
--dr11-colors: the SAME crayons and foregrounds as the completed DES
reproduction, but with their colors measured by the Legacy (Tractor)
reduction instead of DES MOF — the third photometry pipeline on the
same galaxies (SOF vs MOF changed the small-scale amplitude x2).

WHERE TO RUN: astro4 (sweeps on /data/astrofs2_1; I/O bound, ~2-5 h):
    nohup python3 32_rmg_dr11match.py > ~/slurmlogs/igmdust_rdm_match.log 2>&1 &
Output: $DES_DIR/rmg_crayons_dr11.fits
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

HERE = os.path.dirname(os.path.abspath(__file__))
DES_DIR = "/home/suzuki/data/DES/catalog/redmagic"
CRAYON_FILE = "y3a2_gold2.2.1_redmagic_highlum_highz.fits"

ap = argparse.ArgumentParser()
ap.add_argument("--des-dir", default=DES_DIR)
ap.add_argument("--nproc", type=int, default=8)
a = ap.parse_args()

prep = f"{a.des_dir}/rmg_crayons_radec.fits"
out = f"{a.des_dir}/rmg_crayons_dr11.fits"

with fits.open(f"{a.des_dir}/{CRAYON_FILE}") as h:
    d = h[1].data
    t = Table(dict(RA=np.asarray(d["ra"], float),
                   DEC=np.asarray(d["dec"], float),
                   COADD_OBJECT_ID=np.asarray(d["coadd_object_id"],
                                              np.int64),
                   ZREDMAGIC=np.asarray(d["zredmagic"], float)))
t.write(prep, overwrite=True)
print(f"prep: {len(t):,} crayons -> {prep}", flush=True)

rc = subprocess.call(
    [sys.executable, f"{HERE}/06_sweepmatch.py",
     "--in", prep, "--out", out, "--nproc", str(a.nproc)])
sys.exit(rc)
