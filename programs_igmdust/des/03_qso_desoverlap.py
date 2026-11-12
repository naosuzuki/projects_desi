#!/usr/bin/env python3
"""
03_qso_desoverlap.py — igmdust/DES: the apples-to-apples footprint test.

The DES crayon profile (01) sees no reddening beyond 82 kpc, while the
full-footprint QSO profile (projects_desi 02+04) tracks the Menard
extrapolation.  Two explanations remain: QSO-side systematics, or sky
variance (the DES overlap holds only 6% of the BGS sample).  This script
decides between them by re-running the QSO color estimator with the SAME
foregrounds as the crayon run: DR1 BGS galaxies and randoms cut to the
DES Y3 Gold footprint.  QSO backgrounds, colors and cuts are identical
to the main run, and the joint dust+magnification fit (04) is applied to
the result.  If the QSO signal survives here while the crayons see
nothing, it is a QSO systematic; if it vanishes, it was sky variance.

WHERE TO RUN: astrofs4/SLURM (10 cores) or astro4.
    python3 03_qso_desoverlap.py
Output: $OUT/paircount_bgsdes_qso.npz + dustmag_bgsdes.{npz,png}
"""
import argparse
import importlib.util
import os

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "crayon01", os.path.join(_here, "01_crayon_paircount.py"))
m01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m01)

PC04_DEFAULT = os.path.join(os.path.dirname(m01.PC02_DEFAULT),
                            "04_dustmag.py")


def main(a):
    pc = m01.load_pc(a.pc02)
    pc.FG_ZRANGE["bgsdes"] = pc.FG_ZRANGE["bgs"]
    pc.CHUNK = 10_000              # 20G SLURM limit (see 01, same geometry)
    edges, rp_mid = pc.rp_bins()

    good = m01.load_footprint(a.des_dir)
    print("== foregrounds: DESI DR1 BGS, DES overlap")
    fra, fdec, fz, fw = pc.load_lss(a.base, "bgs", "dat")
    m = m01.in_footprint(good, fra, fdec, "BGS")
    fra, fdec, fz, fw = fra[m], fdec[m], fz[m], fw[m]
    print("== randoms: DESI DR1 BGS LSS randoms, DES overlap")
    rra, rdec, rz, rw = pc.load_lss(a.base, "bgs", "ran", nran=a.nran)
    m = m01.in_footprint(good, rra, rdec, "randoms")
    rra, rdec, rz, rw = rra[m], rdec[m], rz[m], rw[m]

    print("== backgrounds: QSO (identical to the main run)")
    qso = pc.load_qso(a.base)
    pc.run_arrays("bgsdes", (fra, fdec, fz, fw), (rra, rdec, rz, rw),
                  qso, edges, rp_mid, a.out,
                  a.max_fg, tag=f"_sub{a.max_fg}" if a.max_fg else "")

    print("\n== joint dust+magnification fit on the overlap profile")
    spec4 = importlib.util.spec_from_file_location("dustmag04", a.pc04)
    m04 = importlib.util.module_from_spec(spec4)
    spec4.loader.exec_module(m04)
    m04.run("bgsdes", a.base, a.out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--pc02", default=m01.PC02_DEFAULT)
    ap.add_argument("--pc04", default=PC04_DEFAULT)
    ap.add_argument("--des-dir", default=m01.DES_DIR_DEFAULT)
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-fg", type=int, default=0)
    ap.add_argument("--nran", type=int, default=2)
    a = ap.parse_args()
    pc_probe = m01.load_pc(a.pc02)
    a.base = a.base or pc_probe.BASE_DEFAULT
    a.out = a.out or pc_probe.OUT_DEFAULT
    main(a)
