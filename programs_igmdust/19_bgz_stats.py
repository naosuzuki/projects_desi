"""19_bgz_stats.py -- BACKGROUND redshift distributions for the
figure insets, companion to 17_fgz_stats.py (which does foregrounds):
weighted n(z), mean, median, 16/84 quantiles for each background
sample, using exactly the catalogs and cuts of the paircount runs.

  descrayon : DES redMaGiC highlum_highz crayons, 0.40 < zredmagic
              < 0.90 (01_crayon_paircount.py CRAYON_ZRANGE) -- the
              backgrounds of the McCleary reproduction (fig:mccrepro).
  lrgcrayon : DESI DR1 LRG spec-z crayons on DR11 photometry
              (07_crayon_legacy.py: lrg_dr11.fits, MATCHED, 0.40-1.10)
              -- left panel of fig:crayondeficit.
  qso       : merged DESI+SDSS+Quaia QSOs with LS DR11 photometry
              (02_paircount.load_qso_merged cuts: MATCHED, PSF,
              0.3 < z < 4.0) -- right panel backgrounds everywhere.

Prints one JSON blob (same schema as 17); the figure scripts
(20_mccrepro_fig.py, 21_crayon_deficit_fig.py) consume it.

Run on astro4 (catalogs on /data + ~/data):
    python3 19_bgz_stats.py
"""
import json

import numpy as np
from astropy.io import fits

OUT = {}


def stats(z, w, zrange, nbin=60):       # identical to 17_fgz_stats.py
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


# DES crayons: redMaGiC highlum_highz, unweighted (01_crayon_paircount
# uses no per-crayon weights; the Y3-Gold color-validity join is color-,
# not redshift-, selective, so zredmagic alone carries the n(z)).
with fits.open("/home/suzuki/data/DES/catalog/redmagic/"
               "y3a2_gold2.2.1_redmagic_highlum_highz.fits") as h:
    OUT["descrayon"] = stats(np.asarray(h[1].data["zredmagic"], float),
                             None, (0.40, 0.90))

EXT = "/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext"

# Legacy LRG spec-z crayons (07_crayon_legacy.load_crayons cuts)
with fits.open(f"{EXT}/lrg_dr11.fits") as h:
    d = h[1].data
    zz = np.asarray(d["Z"], float)
    mm = np.asarray(d["MATCHED"], bool)
    OUT["lrgcrayon"] = stats(zz[mm], None, (0.40, 1.10))

# Merged QSO backgrounds (02_paircount.load_qso_merged cuts)
with fits.open(f"{EXT}/qso_merged_dr11.fits") as h:
    d = h[1].data
    zz = np.asarray(d["Z"], float)
    mm = (np.asarray(d["MATCHED"], bool)
          & (np.char.strip(d["LS_TYPE"].astype(str)) == "PSF"))
    OUT["qso"] = stats(zz[mm], None, (0.30, 4.00), nbin=74)

print(json.dumps(OUT))
