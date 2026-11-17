"""17b_fgz_mcc.py -- weighted n(z) of the Path B foreground (02_paircount
--fg-mcc: BGS passive half & 0.15<z<0.45 & M_r<-20.2), in the JSON
schema of 17_fgz_stats.py, for the figure inset of 34_mccleary_legacy.py.
The selection below mirrors run() in 02_paircount.py line for line.

Run on astro4 (catalogs on /data):
    python3 17b_fgz_mcc.py > fgz_stats_bgsmcc.json
"""
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, f"{HERE}/{fn}")
    m = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [fn]
    try:
        spec.loader.exec_module(m)
    finally:
        sys.argv = argv
    return m


pc = load("paircount02", "02_paircount.py")

stdout, sys.stdout = sys.stdout, sys.stderr      # keep stdout pure JSON
fra, fdec, fz, fw, fm = pc.load_lss(pc.BASE_DEFAULT, "bgs", "dat",
                                    want_mag=True, mag_band="g-r,r")
gr, rmag = fm[:, 0], fm[:, 1]
zlo, zhi = 0.15, 0.45
Mr = (rmag - 25.0 - 5.0 * np.log10(np.clip(pc.dl_interp()(fz), 1e-3, None))
      + 2.5 * np.log10(1.0 + fz))
inz = (fz >= zlo) & (fz <= zhi) & np.isfinite(gr) & (Mr < -20.2)
nzsl = 25
zed = pc.weighted_quantiles(fz[inz], fw[inz], np.linspace(0, 1, nzsl + 1))
isl = np.clip(np.searchsorted(zed, fz) - 1, 0, nzsl - 1)
red = np.zeros(len(fz), bool)
for i in range(nzsl):
    m = inz & (isl == i)
    if m.any():
        red[m] = gr[m] >= pc.weighted_quantiles(gr[m], fw[m], [0.5])[0]
keep = inz & red
z, w = fz[keep], fw[keep]

order = np.argsort(z)
cw = np.cumsum(w[order])
cw /= cw[-1]
h, e = np.histogram(z, bins=60, range=(zlo, zhi), weights=w, density=True)
out = dict(bgsmcc=dict(
    n=int(z.size), mean=float(np.average(z, weights=w)),
    median=float(np.interp(0.5, cw, z[order])),
    q16=float(np.interp(0.16, cw, z[order])),
    q84=float(np.interp(0.84, cw, z[order])),
    zrange=[zlo, zhi],
    Mr_median=float(np.median(Mr[keep])),
    hist=[round(float(x), 4) for x in h],
    edges=[round(float(x), 4) for x in e]))
sys.stdout = stdout
print(json.dumps(out))
