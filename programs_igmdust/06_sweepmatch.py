#!/usr/bin/env python3
"""
06_sweepmatch.py — igmdust: attach Legacy Surveys DR11 photometry to the
merged QSO background (05_qso_merge.py).

Walks the DR11 sweep files (south + north), box-selects the merged QSOs
inside each sweep's RA/Dec block (parsed from the filename), matches them
to sweep sources within 1", and records FLUX_{G,R,I,Z,W1}, FLUX_IVAR_*,
MW_TRANSMISSION_* and TYPE.  Every QSO thus gets photometry from ONE
uniform imaging reduction, whichever catalog it came from.

Startup gate: waits (up to WAIT_MAX_S) for the DR11_SWEEPS_DONE markers
of both download shards (00b on astro4 + 03b SLURM job), then proceeds —
with a warning if a marker never appeared (QSOs on missing sweeps simply
stay unmatched and are flagged).

WHERE TO RUN: astrofs4/SLURM (10 cores; I/O bound, ~2-5 h over NFS):
    sbatch -d afterok:<dl job>,singleton 03_slurm_paircount.sbatch sweepmatch
On astro4 ALWAYS pass --mbps: unpaced readers (8, and 3 as well) saturate
its link to astrofs2 and stall every login on the host (2026-11-16).
A killed run resumes from <out>.ckpt.npz (written every 15 min).
Output: $EXT/qso_merged_dr11.fits
"""
import argparse
import glob
import os
import re
import sys
import time
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
from astropy.table import Table

EXT_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"
LS_DEFAULT = "/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11"
MARKERS = [os.path.expanduser("~/dr11_dl_shard1.log"),
           *sorted(glob.glob(os.path.expanduser(
               "~/slurmlogs/igmdust_dl_*.log")))]
WAIT_MAX_S = 4 * 3600
MATCH_ARCSEC = 1.0
BANDS = ["G", "R", "I", "Z", "W1"]
NPROC = 6
CKPT_S = 900
NAME_RE = re.compile(
    r"sweep-(\d{3})([mp])(\d{3})-(\d{3})([mp])(\d{3})\.fits$")


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    c = np.cos(dec)
    return np.column_stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)])


def parse_bounds(fname):
    """sweep-290m035-295m030.fits -> (290, -35, 295, -30)."""
    m = NAME_RE.search(os.path.basename(fname))
    if not m:
        return None
    r1, s1, d1, r2, s2, d2 = m.groups()
    sgn = {"m": -1, "p": 1}
    return (int(r1), sgn[s1] * int(d1), int(r2), sgn[s2] * int(d2))


def wait_for_downloads():
    t0 = time.time()
    while time.time() - t0 < WAIT_MAX_S:
        done = [f for f in MARKERS if os.path.exists(f)
                and "DR11_SWEEPS_DONE" in open(f, errors="ignore").read()]
        if len(done) >= 2:
            print(f"  both download shards done ({len(done)} markers)")
            return
        print(f"  waiting for download markers ({len(done)}/2) "
              f"[{(time.time()-t0)/60:.0f} min]", flush=True)
        time.sleep(300)
    print("  WARNING: proceeding without all DONE markers; "
          "unmatched QSOs will be flagged")


_qso = {}


def _init(ra, dec, mbps=0.0):
    _qso["ra"], _qso["dec"] = ra, dec
    _qso["vec"] = unitvec(ra, dec)
    _qso["mbps"] = mbps


def paced_read(fname, mbps):
    """Pull fname into the page cache at <= mbps MB/s.

    The memmap pass in match_one faults in the whole file as fast as NFS
    delivers it; a few such readers saturate the link to the file server
    and every other NFS client operation on the host (home directories
    included) queues behind them.  Reading in short paced bursts first
    leaves gaps for everyone else; the memmap pass then hits the cache.
    """
    blk = 8 << 20
    buf = bytearray(blk)
    t0, n = time.time(), 0
    with open(fname, "rb", buffering=0) as f:
        while True:
            k = f.readinto(buf)
            if not k:
                break
            n += k
            lag = n / (mbps * 1e6) - (time.time() - t0)
            if lag > 0:
                time.sleep(lag)


def _job(fname):
    return os.path.basename(fname), match_one(fname)


def match_one(fname):
    """Match QSOs inside this sweep's block; return (idx, dist, cols)."""
    b = parse_bounds(fname)
    if b is None:
        return None
    ra, dec = _qso["ra"], _qso["dec"]
    sel = np.where((ra >= b[0]) & (ra < b[2])
                   & (dec >= b[1]) & (dec < b[3]))[0]
    if len(sel) == 0:
        return None
    from scipy.spatial import cKDTree
    chord = 2 * np.sin(np.radians(MATCH_ARCSEC / 3600) / 2)
    if _qso.get("mbps"):
        paced_read(fname, _qso["mbps"])
    with fits.open(fname, memmap=True) as h:
        d = h[1].data
        names = {c.upper() for c in h[1].columns.names}
        svec = unitvec(np.asarray(d["RA"], float),
                       np.asarray(d["DEC"], float))
        dist, qi = cKDTree(_qso["vec"][sel]).query(
            svec, k=1, distance_upper_bound=chord)
        hit = np.isfinite(dist)
        if not hit.any():
            return None
        rows = np.where(hit)[0]
        # best (nearest) sweep row per QSO
        order = np.argsort(dist[rows])
        rows, qrows = rows[order], qi[rows][order]
        first = np.unique(qrows, return_index=True)[1]
        rows, qrows = rows[first], qrows[first]
        out = {"idx": sel[qrows],
               "dist": np.degrees(2 * np.arcsin(dist[rows] / 2)) * 3600,
               "type": np.char.strip(d["TYPE"][rows].astype(str))}
        for b_ in BANDS:
            for stem, key in ((f"FLUX_{b_}", f"flux_{b_}"),
                              (f"FLUX_IVAR_{b_}", f"ivar_{b_}"),
                              (f"MW_TRANSMISSION_{b_}", f"mwt_{b_}")):
                out[key] = (np.asarray(d[stem][rows], float)
                            if stem in names
                            else np.full(len(rows), np.nan))
    return out


def run(a):
    wait_for_downloads()
    inp = a.inp or f"{a.ext}/qso_merged.fits"
    t0 = time.time()
    while not os.path.exists(inp) and time.time() - t0 < 7200:
        print(f"  waiting for {inp}", flush=True)
        time.sleep(300)
    q = Table.read(inp)
    n = len(q)
    print(f"  merged QSOs: {n:,}")
    sweeps = sorted(glob.glob(f"{a.ls}/south/sweep-*.fits")
                    + glob.glob(f"{a.ls}/north/sweep-*.fits"))
    sweeps = [s for s in sweeps if parse_bounds(s)]
    print(f"  sweeps: {len(sweeps)}")

    ra = np.asarray(q["RA"], float)
    dec = np.asarray(q["DEC"], float)
    res = {"dist": np.full(n, np.inf), "type": np.zeros(n, "U4")}
    for b in BANDS:
        for k in (f"flux_{b}", f"ivar_{b}", f"mwt_{b}"):
            res[k] = np.full(n, np.nan)

    fo = a.out or f"{a.ext}/qso_merged_dr11.fits"
    ck = fo + ".ckpt.npz"
    done = []
    if os.path.exists(ck):
        with np.load(ck) as z:
            if int(z["n"]) == n:
                for k in res:
                    res[k] = z[k]
                done = z["done"].tolist()
                print(f"  resume: {len(done)} sweeps already done ({ck})",
                      flush=True)
            else:
                print(f"  ignoring {ck}: made for a different input")
    seen = set(done)
    todo = [s for s in sweeps if os.path.basename(s) not in seen]

    def checkpoint():
        tmp = fo + ".ckpt.tmp.npz"
        np.savez(tmp, n=n, done=np.array(done), **res)
        os.replace(tmp, ck)

    t0, tck = time.time(), time.time()
    per = a.mbps / a.nproc if a.mbps else 0.0
    with Pool(a.nproc, initializer=_init, initargs=(ra, dec, per)) as pool:
        for name, out in pool.imap_unordered(
                _job, todo, chunksize=1 if a.mbps else 2):
            if out is not None:
                take = out["dist"] < res["dist"][out["idx"]]
                ix = out["idx"][take]
                for k in res:
                    if k != "dist":
                        res[k][ix] = out[k][take]
                res["dist"][ix] = out["dist"][take]
            done.append(name)
            if len(done) % 100 == 0:
                print(f"  [{len(done)}/{len(sweeps)}] matched so far "
                      f"{(np.isfinite(res['dist'])).sum():,} "
                      f" {time.time()-t0:.0f}s", flush=True)
            if time.time() - tck > CKPT_S:
                checkpoint()
                tck = time.time()

    matched = np.isfinite(res["dist"])
    print(f"  matched {matched.sum():,} / {n:,} "
          f"({matched.mean():.1%}); median sep "
          f"{np.median(res['dist'][matched]):.3f}\"")
    for k, v in res.items():
        q[k.upper() if k != "type" else "LS_TYPE"] = v
    q["MATCHED"] = matched
    q.write(fo, overwrite=True)
    print(f"saved -> {fo}")
    if os.path.exists(ck):
        os.remove(ck)


def selftest():
    """Name parsing + box-select + nearest-match on a synthetic sweep."""
    assert parse_bounds("sweep-290m035-295m030.fits") == (290, -35, 295, -30)
    assert parse_bounds("sweep-000p005-005p010.fits") == (0, 5, 5, 10)
    assert parse_bounds("sweep-foo.fits") is None
    rng = np.random.default_rng(2)
    nq = 500
    qra = rng.uniform(290.1, 294.9, nq)
    qdec = rng.uniform(-34.9, -30.1, nq)
    _init(qra, qdec)
    # synthetic sweep: every even QSO present (offset 0.4"), plus field junk
    off = 0.4 / 3600
    sra = np.concatenate([qra[::2] + off * 0.7, rng.uniform(290, 295, 5000)])
    sdec = np.concatenate([qdec[::2] + off * 0.7, rng.uniform(-35, -30, 5000)])
    cols = dict(RA=sra, DEC=sdec, TYPE=np.full(len(sra), "PSF"))
    for b in BANDS:
        cols[f"FLUX_{b}"] = np.arange(len(sra), dtype=float)
        cols[f"FLUX_IVAR_{b}"] = np.ones(len(sra))
        cols[f"MW_TRANSMISSION_{b}"] = np.ones(len(sra))
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        f = f"{td}/sweep-290m035-295m030.fits"
        Table(cols).write(f)
        out = match_one(f)
    ok = (out is not None and len(out["idx"]) == 250
          and np.array_equal(np.sort(out["idx"]), np.arange(0, nq, 2))
          and np.array_equal(out["flux_G"][np.argsort(out["idx"])],
                             np.arange(250.0)))
    print(f"  matched {0 if out is None else len(out['idx'])}/250 expected")
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ext", default=EXT_DEFAULT)
    ap.add_argument("--in", dest="inp", default=None,
                    help="input catalog with RA/DEC (default: the merged "
                         "QSO list); any table works, e.g. LRG positions")
    ap.add_argument("--out", default=None)
    ap.add_argument("--ls", default=LS_DEFAULT)
    ap.add_argument("--nproc", type=int, default=NPROC)
    ap.add_argument("--mbps", type=float, default=0.0,
                    help="total read-rate cap in MB/s (0 = unpaced); set "
                         "it on a shared host whose NFS link the readers "
                         "would otherwise saturate")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    sys.exit(selftest() if a.selftest else run(a))
