#!/usr/bin/env python3
"""
04_redmapper_crayons.py — igmdust/DES: dust around redMaPPer clusters,
measured with redMaGiC standard crayons (the McCleary massive-halo case).

Foregrounds: DES Y3 redMaPPer v6.4.22+2 clusters, the official lgt20
(lambda > 20) cosmology selection from the release file's index group
(19,388 clusters), with the release's own weighted randoms (ztrue).
Backgrounds: the redMaGiC highlum_highz standard crayons of 01 (locus
residuals about zredmagic; multiply profiles by 1/f_abs = 1/0.806 from
02_speczcalib when quoting de-absorbed reddening).

Estimator and bins identical to the other crayon runs (projects_desi 02
machinery; physical r_p at the cluster z_lambda; z buffer 0.2).

WHERE TO RUN: astrofs4/SLURM (10 cores, ~15 min) or astro4:
    sbatch -J igmdust_rm -d singleton 03_slurm_paircount.sbatch redmapper
Output: $OUT/paircount_redmapper_qso.npz
"""
import argparse
import importlib.util
import os

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "crayon01", os.path.join(_here, "01_crayon_paircount.py"))
m01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m01)

H5_DEFAULT = ("/home/suzuki/data/DES/catalog/redmapper/"
              "y3_redmapper_v6.4.22+2_release.h5")
INDEX = "index/redmapper/lgt20"      # official lambda>20 cosmology sample
FG_ZRANGE_RM = (0.10, 0.70)          # crayons are 0.4-0.9 with 0.2 buffer
MAX_RAN = 1_500_000                  # plenty: ~75x the cluster count


def load_clusters(h5path, lam_min=0.0, lam_max=np.inf):
    """lgt20 selection, optionally further cut in richness; the release
    randoms carry the richness they were drawn for (lambda_in), so the
    same cut applies to them."""
    import h5py
    with h5py.File(h5path, "r") as f:
        sel = f[f"{INDEX}/select"][:]
        g = f["catalog/cluster"]
        fra = g["ra"][:][sel]
        fdec = g["dec"][:][sel]
        fz = np.asarray(g["z_lambda"][:][sel], float)
        lam = np.asarray(g["lambda_chisq"][:][sel], float)
        keep = (lam >= lam_min) & (lam < lam_max)
        fra, fdec, fz, lam = fra[keep], fdec[keep], fz[keep], lam[keep]
        print(f"  clusters (lgt20, {lam_min}<=lambda<{lam_max}): "
              f"{keep.sum():,}  z_lambda [{fz.min():.2f},{fz.max():.2f}]  "
              f"lambda [{lam.min():.0f},{lam.max():.0f}]")
        rsel = f[f"{INDEX}/random_select"][:]
        r = f["randoms"]
        rlam = np.asarray(r["lambda_in"][:][rsel], float)
        rkeep = (rlam >= lam_min) & (rlam < lam_max)
        rra = r["ra"][:][rsel][rkeep]
        rdec = r["dec"][:][rsel][rkeep]
        rz = np.asarray(r["ztrue"][:][rsel], float)[rkeep]
        rw = np.asarray(r["weight"][:][rsel], float)[rkeep]
    if len(rra) > MAX_RAN:
        rng = np.random.default_rng(3)
        pick = rng.choice(len(rra), MAX_RAN, replace=False)
        rra, rdec, rz, rw = rra[pick], rdec[pick], rz[pick], rw[pick]
    print(f"  randoms (release, weighted, same lambda cut): {len(rra):,}")
    return ((fra, fdec, fz, np.ones(len(fra))), (rra, rdec, rz, rw))


def main(a):
    pc = m01.load_pc(a.pc02)
    pc.R_COEF["i"] = m01.R_I
    pc.COLORS = m01.COLORS_DES
    pc.JK_NSIDE = m01.JK_NSIDE
    name = "redmapper" + a.tag
    pc.FG_ZRANGE[name] = FG_ZRANGE_RM
    pc.RP_MIN = a.rp_min
    pc.RP_MAX = a.rp_max
    # 20G SLURM limit (see 01); beyond 10 Mpc the low-z query radii grow
    # to ~5 deg, so halve the chunk again.
    pc.CHUNK = 5_000 if a.rp_max > 10 else 10_000
    edges, rp_mid = pc.rp_bins()

    print("== crayons: DES redMaGiC highlum_highz + Y3 Gold colors")
    cra, cdec, cz, cres, cval = m01.load_crayons(a.des_dir)
    print("== foregrounds: DES Y3 redMaPPer clusters")
    fg, ran = load_clusters(a.h5, a.lambda_min, a.lambda_max)

    pc.run_arrays(name, fg, ran, (cra, cdec, cz, cres, cval),
                  edges, rp_mid, a.out,
                  a.max_fg, tag=f"_sub{a.max_fg}" if a.max_fg else "")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--h5", default=H5_DEFAULT)
    ap.add_argument("--des-dir", default=m01.DES_DIR_DEFAULT)
    ap.add_argument("--pc02", default=m01.PC02_DEFAULT)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-fg", type=int, default=0)
    ap.add_argument("--lambda-min", type=float, default=0.0)
    ap.add_argument("--lambda-max", type=float, default=np.inf)
    ap.add_argument("--rp-min", type=float, default=0.02)
    ap.add_argument("--rp-max", type=float, default=10.0,
                    help="outer radius [Mpc]; 45 covers 30 h^-1 Mpc com.")
    ap.add_argument("--tag", default="",
                    help="suffix for the output profile name")
    a = ap.parse_args()
    a.out = a.out or m01.load_pc(a.pc02).OUT_DEFAULT
    main(a)
