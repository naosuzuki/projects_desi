#!/usr/bin/env python3
"""
01_build_manifests.py
----------------------
Build Globus transfer manifests (and per-type target lists) for collecting
DESI DR2 ("loa") spectra:

  * STAR + GALAXY spectra from the HEALPIX coadds   (zall-pix-loa.fits)
  * QSO spectra from the TILE cumulative coadds      (zall-tilecumulative-loa.fits)

Selection (per the agreed scope):
  * all surveys / programs            (no SURVEY/PROGRAM filter)
  * ZCAT_PRIMARY == True              (one best row per target, no duplicates)
  * ZWARN == 0                        (good redshift / classification only)

Outputs (written to OUTDIR):
  healpix_coadd_files.batch    SRC DST lines for `globus transfer --batch`
  tiles_coadd_files.batch      SRC DST lines for `globus transfer --batch`
  healpix_star_targets.csv     TARGETID,SURVEY,PROGRAM,HEALPIX
  healpix_galaxy_targets.csv   TARGETID,SURVEY,PROGRAM,HEALPIX
  tiles_qso_targets.csv        TARGETID,TILEID,PETAL_LOC,LASTNIGHT
  manifest_summary.txt         counts + estimated transfer size

Nothing is downloaded here. Review manifest_summary.txt, then transfer with
02_transfer_spectra.sh.
"""
import os
import sys
import numpy as np
import fitsio

# --- Paths -------------------------------------------------------------------
HEALPIX_CAT = os.environ.get("ZPIX_CAT",  "/Users/suzuki/desiredux/zall-pix-loa.fits")
TILE_CAT    = os.environ.get("ZTILE_CAT", "/Volumes/exdisk2/desiredux/DR2/zall-tilecumulative-loa.fits")
OUTDIR      = os.environ.get("OUTDIR",    os.path.join(os.path.dirname(os.path.abspath(__file__)), "manifests"))

# Remote (NERSC) and local (exdisk2) roots. The .batch files hold paths
# RELATIVE to these roots, which are passed to `globus transfer` on the CLI.
REDUX_REL   = "spectro/redux/loa"                       # under each endpoint root
NERSC_BASE  = "/global/cfs/cdirs/desi/public/dr2/" + REDUX_REL
LOCAL_BASE  = os.environ.get("LOCAL_BASE", "/Volumes/exdisk2/desiredux/DR2/" + REDUX_REL)

# Read the catalog in row-chunks to bound memory on the 64M-row table.
CHUNK = 5_000_000

os.makedirs(OUTDIR, exist_ok=True)


def _decode(arr):
    """bytes -> str for FITS string columns; leave others alone."""
    if arr.dtype.kind == "S":
        return np.char.strip(arr.astype("U"))
    return arr


def read_columns(path, hdu, cols):
    """Read selected columns from a (possibly huge) FITS binary table in chunks."""
    f = fitsio.FITS(path)
    h = f[hdu]
    n = h.get_nrows()
    out = {c: [] for c in cols}
    for lo in range(0, n, CHUNK):
        hi = min(lo + CHUNK, n)
        rows = h[cols][lo:hi]
        for c in cols:
            out[c].append(_decode(rows[c]))
        print(f"    read rows {lo:,}–{hi:,} / {n:,}", file=sys.stderr)
    f.close()
    return {c: np.concatenate(v) for c, v in out.items()}


def write_target_csv(path, header, columns):
    """Write a small CSV of selected targets (columns: list of equal-length arrays)."""
    arrs = list(zip(*columns))  # not used; kept explicit below
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for row in zip(*columns):
            fh.write(",".join(str(x) for x in row) + "\n")


# =============================================================================
# HEALPIX: STAR + GALAXY
# =============================================================================
def build_healpix():
    print(f"[healpix] reading {HEALPIX_CAT}", file=sys.stderr)
    cols = ["TARGETID", "SURVEY", "PROGRAM", "HEALPIX",
            "SPECTYPE", "ZWARN", "ZCAT_PRIMARY"]
    d = read_columns(HEALPIX_CAT, "ZCATALOG", cols)

    good = (d["ZCAT_PRIMARY"]) & (d["ZWARN"] == 0)
    is_star = good & (d["SPECTYPE"] == "STAR")
    is_gal  = good & (d["SPECTYPE"] == "GALAXY")
    print(f"[healpix] primary+ZWARN0: {good.sum():,}  "
          f"STAR={is_star.sum():,}  GALAXY={is_gal.sum():,}", file=sys.stderr)

    # per-type target lists (for the later extraction step)
    for name, sel in (("star", is_star), ("galaxy", is_gal)):
        write_target_csv(
            os.path.join(OUTDIR, f"healpix_{name}_targets.csv"),
            ["TARGETID", "SURVEY", "PROGRAM", "HEALPIX"],
            [d["TARGETID"][sel], d["SURVEY"][sel], d["PROGRAM"][sel], d["HEALPIX"][sel]],
        )

    # unique coadd files = union of star+galaxy locations (coadds mix types)
    sel = is_star | is_gal
    keys = np.stack([d["SURVEY"][sel], d["PROGRAM"][sel],
                     d["HEALPIX"][sel].astype("U12")], axis=1)
    uniq = np.unique(keys, axis=0)
    print(f"[healpix] unique coadd files: {len(uniq):,}", file=sys.stderr)

    batch = os.path.join(OUTDIR, "healpix_coadd_files.batch")
    with open(batch, "w") as fh:
        for survey, program, hpx in uniq:
            hpx_i = int(hpx)
            rel = (f"healpix/{survey}/{program}/{hpx_i // 100}/{hpx_i}/"
                   f"coadd-{survey}-{program}-{hpx_i}.fits")
            fh.write(f"{NERSC_BASE}/{rel} {LOCAL_BASE}/{rel}\n")
    return {"n_files": len(uniq), "n_star": int(is_star.sum()),
            "n_galaxy": int(is_gal.sum()), "batch": batch}


# =============================================================================
# TILES (cumulative): QSO
# =============================================================================
def build_tiles():
    if not os.path.exists(TILE_CAT):
        print(f"[tiles] catalog not present yet: {TILE_CAT} — skipping.", file=sys.stderr)
        return None
    print(f"[tiles] reading {TILE_CAT}", file=sys.stderr)
    # NOTE: verify these column names against the tile catalog header.
    cols = ["TARGETID", "TILEID", "PETAL_LOC", "LASTNIGHT",
            "SPECTYPE", "ZWARN", "ZCAT_PRIMARY"]
    d = read_columns(TILE_CAT, "ZCATALOG", cols)

    is_qso = (d["ZCAT_PRIMARY"]) & (d["ZWARN"] == 0) & (d["SPECTYPE"] == "QSO")
    print(f"[tiles] QSO (primary+ZWARN0): {is_qso.sum():,}", file=sys.stderr)

    write_target_csv(
        os.path.join(OUTDIR, "tiles_qso_targets.csv"),
        ["TARGETID", "TILEID", "PETAL_LOC", "LASTNIGHT"],
        [d["TARGETID"][is_qso], d["TILEID"][is_qso],
         d["PETAL_LOC"][is_qso], d["LASTNIGHT"][is_qso]],
    )

    keys = np.stack([d["TILEID"][is_qso].astype("U12"),
                     d["PETAL_LOC"][is_qso].astype("U2"),
                     d["LASTNIGHT"][is_qso].astype("U8")], axis=1)
    uniq = np.unique(keys, axis=0)
    print(f"[tiles] unique coadd (petal) files: {len(uniq):,}", file=sys.stderr)

    batch = os.path.join(OUTDIR, "tiles_coadd_files.batch")
    with open(batch, "w") as fh:
        for tileid, petal, night in uniq:
            rel = (f"tiles/cumulative/{int(tileid)}/{night}/"
                   f"coadd-{int(petal)}-{int(tileid)}-thru{night}.fits")
            fh.write(f"{NERSC_BASE}/{rel} {LOCAL_BASE}/{rel}\n")
    return {"n_files": len(uniq), "n_qso": int(is_qso.sum()), "batch": batch}


def main():
    hp = build_healpix()
    tl = build_tiles()

    lines = ["DESI DR2 (loa) spectra collection — manifest summary", "=" * 52, ""]
    lines += [f"HEALPIX (STAR+GALAXY):",
              f"  STAR targets   : {hp['n_star']:,}",
              f"  GALAXY targets : {hp['n_galaxy']:,}",
              f"  coadd files    : {hp['n_files']:,}   -> {hp['batch']}", ""]
    if tl:
        lines += [f"TILES cumulative (QSO):",
                  f"  QSO targets    : {tl['n_qso']:,}",
                  f"  coadd files    : {tl['n_files']:,}   -> {tl['batch']}", ""]
    else:
        lines += ["TILES cumulative (QSO): SKIPPED (catalog not downloaded yet)", ""]
    lines += ["Next: run 02_estimate_size.sh to sample real file sizes,",
              "then 02_transfer_spectra.sh to start the Globus transfer."]

    summary = os.path.join(OUTDIR, "manifest_summary.txt")
    with open(summary, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote manifests to {OUTDIR}")


if __name__ == "__main__":
    main()
