#!/usr/bin/env python3
"""
00b_fetch_y3gold_colors.py — igmdust/DES: fetch griz photometry for the Y3A2
redMaGiC highlum_highz "crayon" galaxies from Y3 GOLD via the NOIRLab Astro
Data Lab TAP service (anonymous, sync).

WHY: the redMaGiC value-added catalog carries only mag_z + zredmagic; the
reddening estimator needs griz colors. redMaGiC coadd_object_ids are Y3
Gold ids (DR1 coadd namespace — they do NOT exist in des_dr2), so we pull
sof_cm_mag_corrected_[griz] (dereddened SOF composite-model photometry,
the standard Y3 choice) from des_dr1.y3_gold (850,028 rows -> ~100 MB).

WHERE TO RUN: astro4 (direct internet; this is a network job, NOT a cluster
job). From ~/github/projects_desi/programs_igmdust/des:
    nohup nice -n 19 python3 00b_fetch_y3gold_colors.py > ~/des_y3gold_colors.log 2>&1 &

Chunked (5,000 ids per query, 4 parallel workers, 3 retries each); CSV part
files land in $DEST/parts_dr2colors/ (the Data Lab sync TAP does not serve
FITS) and are stitched to $DEST/redmagic_highlum_highz_y3gold_colors.fits at
the end. Idempotent: existing part files are skipped on rerun.
"""
import os
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

DEST = "/home/suzuki/data/DES/catalog/redmagic"
SRC = f"{DEST}/y3a2_gold2.2.1_redmagic_highlum_highz.fits"
PARTS = f"{DEST}/parts_dr2colors"
OUT = f"{DEST}/redmagic_highlum_highz_y3gold_colors.fits"
TAP = "https://datalab.noirlab.edu/tap/sync"
CHUNK = 5_000
WORKERS = 4

TABLE = "des_dr1.y3_gold"
COLS_DERED = ("coadd_object_id, "
              "sof_cm_mag_corrected_g, sof_cm_mag_corrected_r, "
              "sof_cm_mag_corrected_i, sof_cm_mag_corrected_z, "
              "sof_cm_mag_err_g, sof_cm_mag_err_r, "
              "sof_cm_mag_err_i, sof_cm_mag_err_z, "
              "flags_gold, extended_class_sof")
COLS_PLAIN = ("coadd_object_id, mag_auto_g, mag_auto_r, mag_auto_i, "
              "mag_auto_z, magerr_auto_g, magerr_auto_r, magerr_auto_i, "
              "magerr_auto_z, flags_gold, extended_class_sof")
COLS_MOF = ("coadd_object_id, "
            "mof_cm_mag_corrected_g, mof_cm_mag_corrected_r, "
            "mof_cm_mag_corrected_i, mof_cm_mag_corrected_z, "
            "mof_cm_mag_err_g, mof_cm_mag_err_r, "
            "mof_cm_mag_err_i, mof_cm_mag_err_z, "
            "flags_gold, extended_class_sof")

# --mof: fetch MOF photometry (McCleary's choice) to a separate file;
# selection columns (flags_gold, extended_class_sof) stay identical so
# only the photometry changes between iterations.
MOF = "--mof" in sys.argv
if MOF:
    PARTS = f"{DEST}/parts_mofcolors"
    OUT = f"{DEST}/redmagic_highlum_highz_y3gold_mofcolors.fits"


def tap_query(adql, timeout=300):
    """Run a sync TAP query, return the decoded CSV text (validated)."""
    data = urllib.parse.urlencode(dict(
        REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql)).encode()
    req = urllib.request.Request(TAP, data=data)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        text = r.read().decode("utf-8", "replace")
    if text.lstrip().startswith("<?xml") or "QUERY_STATUS" in text[:2000]:
        raise RuntimeError(f"TAP error: {text[:300]}")
    return text


def pick_columns():
    """Probe whether the _dered columns exist; fall back to plain."""
    for cols in ((COLS_MOF,) if MOF else (COLS_DERED, COLS_PLAIN)):
        try:
            t = Table.read(tap_query(f"SELECT TOP 1 {cols} FROM {TABLE}",
                                     timeout=120), format="ascii.csv")
            assert len(t.colnames) >= 10
            lab = ("MOF corrected" if cols is COLS_MOF else
                   "SOF corrected" if cols is COLS_DERED else "mag_auto")
            print(f"columns OK: {lab}")
            return cols
        except Exception as e:                          # noqa: BLE001
            print(f"probe failed ({e}); trying fallback")
    sys.exit(f"no usable column set on {TABLE}")


def fetch_chunk(args):
    i, ids, cols = args
    part = f"{PARTS}/part_{i:04d}.csv"
    if os.path.exists(part) and os.path.getsize(part) > 0:
        return f"{i}: cached"
    idlist = ",".join(map(str, ids))
    adql = (f"SELECT {cols} FROM {TABLE} "
            f"WHERE coadd_object_id IN ({idlist})")
    for attempt in range(3):
        try:
            text = tap_query(adql)
            n = len(Table.read(text, format="ascii.csv"))        # validate
            with open(part + ".tmp", "w") as f:
                f.write(text)
            os.replace(part + ".tmp", part)
            return f"{i}: {n} rows"
        except Exception as e:                          # noqa: BLE001
            time.sleep(10 * (attempt + 1))
            err = e
    return f"{i}: FAILED after 3 tries ({err})"


def main():
    os.makedirs(PARTS, exist_ok=True)
    ids = np.sort(np.asarray(fits.getdata(SRC)["coadd_object_id"], np.int64))
    print(f"{len(ids):,} redMaGiC ids -> {int(np.ceil(len(ids)/CHUNK))} chunks")
    cols = pick_columns()
    jobs = [(i, ids[i * CHUNK:(i + 1) * CHUNK], cols)
            for i in range(int(np.ceil(len(ids) / CHUNK)))]
    t0, nfail = time.time(), 0
    with ThreadPoolExecutor(WORKERS) as ex:
        for k, msg in enumerate(ex.map(fetch_chunk, jobs)):
            nfail += "FAILED" in msg
            if k % 10 == 0 or "FAILED" in msg:
                print(f"[{k+1}/{len(jobs)}] {msg}  {time.time()-t0:.0f}s",
                      flush=True)
    if nfail:
        sys.exit(f"{nfail} chunks failed -- rerun to retry (parts cached)")

    print("stitching...")
    tabs = [Table.read(f"{PARTS}/part_{i:04d}.csv", format="ascii.csv")
            for i in range(len(jobs))]
    full = vstack(tabs, metadata_conflicts="silent")
    full.write(OUT, overwrite=True)
    print(f"{len(full):,} / {len(ids):,} matched -> {OUT}")
    print(f"DES_DR2COLORS_DONE {time.strftime('%F %T')}")


if __name__ == "__main__":
    main()
