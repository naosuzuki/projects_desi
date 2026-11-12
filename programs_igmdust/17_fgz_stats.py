"""17_fgz_stats.py -- foreground redshift distributions for the profile
figures: weighted n(z), mean, median, 16/84 quantiles per foreground
sample (DESI BGS/LRG, DES redMaGiC highdens, DES redMaPPer lgt20),
using exactly the z ranges and weights of the paircount runs.

Prints one JSON blob; the figure script consumes it for insets and
annotations (DES vs DESI Legacy foreground-redshift comparison).

Run on astro4 (catalogs on /data + ~/data):
    python3 17_fgz_stats.py
"""
import json

import numpy as np
from astropy.io import fits

OUT = {}


def stats(z, w, zrange, nbin=60):
    m = (z >= zrange[0]) & (z <= zrange[1]) & np.isfinite(z)
    z = z[m]
    w = np.ones(z.size) if w is None else np.asarray(w, float)[m]
    order = np.argsort(z)
    zs, ws = z[order], w[order]
    cw = np.cumsum(ws)
    cw /= cw[-1]
    h, e = np.histogram(z, bins=nbin, range=zrange, weights=w, density=True)
    return dict(n=int(z.size),
                mean=float(np.average(z, weights=w)),
                median=float(np.interp(0.5, cw, zs)),
                q16=float(np.interp(0.16, cw, zs)),
                q84=float(np.interp(0.84, cw, zs)),
                zrange=list(zrange),
                hist=[round(float(x), 4) for x in h],
                edges=[round(float(x), 4) for x in e])


# DESI DR1 LSS foregrounds (clustering weights), paircount z cuts
BASE = "/data/astrofs2_1/suzuki/data/DESI/catalog/dr1"
for name, stem, zr in (("bgs", "BGS_ANY", (0.05, 0.45)),
                       ("lrg", "LRG", (0.40, 1.10))):
    zz, ww = [], []
    for cap in ("NGC", "SGC"):
        with fits.open(f"{BASE}/{stem}_{cap}_clustering.dat.fits") as h:
            d = h[1].data
            zz.append(np.asarray(d["Z"], float))
            try:
                ww.append(np.asarray(d["WEIGHT"], float))
            except KeyError:
                ww.append(np.ones(len(zz[-1])))
    OUT[name] = stats(np.concatenate(zz), np.concatenate(ww), zr)

# DES redMaGiC high-density foregrounds (systematics weights), 0.15-0.45
with fits.open("/home/suzuki/data/DES/catalog/redmagic/"
               "y3a2_gold2.2.1_redmagic_highdens.fits") as h:
    d = h[1].data
    OUT["redmagic"] = stats(np.asarray(d["zredmagic"], float),
                            np.asarray(d["weight"], float), (0.15, 0.45))

# DES redMaPPer lgt20 cosmology clusters, run z range 0.10-0.70
import h5py  # noqa: E402  (heavy import last)
with h5py.File("/home/suzuki/data/DES/catalog/redmapper/"
               "y3_redmapper_v6.4.22+2_release.h5", "r") as f:
    sel = f["index/redmapper/lgt20/select"][:]
    g = f["catalog/cluster"]
    OUT["redmapper"] = stats(np.asarray(g["z_lambda"][:][sel], float),
                             None, (0.10, 0.70))

# Menard replica: photometric foregrounds carry a single effective z
OUT["photgal"] = dict(n=135952763, zeff=0.30,
                      note="r<21 sweep galaxies; single z_eff by design")

print(json.dumps(OUT))
