#!/usr/bin/env python3
"""
06_csfd.py — igmdust/DES: CSFD-vs-SFD Milky Way re-correction (pure numpy).

The Y3 Gold SOF *corrected* magnitudes subtract SFD98 extinction.  SFD's
unresolved extragalactic contamination (CIB + galaxy light) means that
subtraction removes real galaxy-correlated dust within the SFD beam
(6.1 arcmin ~ 1.6 Mpc at z~0.3) of every foreground --- the prime
suspect for the 0.15-1.4 Mpc deficit against McCleary's profile, who
used the decontaminated CSFD map (Chiang 2023).

delta_ebv(ra, dec) returns E_SFD - E_CSFD at the given equatorial
positions; adding dR_color * delta_ebv to a color re-expresses the
SFD-corrected photometry in the CSFD convention.

Maps (astro4: /home/suzuki/data/dustmaps):
  csfd_ebv.fits        nside 2048 RING, galactic
  lambda_sfd_ebv.fits  nside  512 NESTED, galactic ('TEMPERATURE'=E(B-V))
No healpy: RING ang2pix from 01 (validated bit-identical) plus the
NESTED version below (validated against healpy on 10^6 points).
"""
import importlib.util
import os

import numpy as np
from astropy.io import fits

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "crayon01", os.path.join(_here, "01_crayon_paircount.py"))
_m01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m01)
ang2pix_ring = _m01.ang2pix_ring

MAPDIR = "/home/suzuki/data/dustmaps"
# DES Y3 fiducial A_b / E(B-V) used by the SOF 'corrected' dereddening
R_DES = dict(g=3.186, r=2.140, i=1.569, z=1.196)

# J2000 equatorial -> galactic rotation
_EQ2GAL = np.array([
    [-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
    [+0.4941094278755837, -0.4448296299600112, +0.7469822444972189],
    [-0.8676661490190047, -0.1980763734312015, +0.4559837761750669]])


def eq2gal(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    c = np.cos(dec)
    v = np.stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)])
    g = _EQ2GAL @ v
    el = np.degrees(np.arctan2(g[1], g[0])) % 360.0
    b = np.degrees(np.arcsin(np.clip(g[2], -1, 1)))
    return el, b


def _spread_bits(x):
    """Interleave zeros between the bits of x (x < 2^16)."""
    x = np.asarray(x, np.int64)
    x = (x | (x << 8)) & 0x00FF00FF
    x = (x | (x << 4)) & 0x0F0F0F0F
    x = (x | (x << 2)) & 0x33333333
    x = (x | (x << 1)) & 0x55555555
    return x


def ang2pix_nest(nside, ra_deg, dec_deg):
    """HEALPix NESTED ang2pix in pure numpy (healpy algorithm)."""
    z = np.sin(np.radians(dec_deg))
    tt = np.mod(np.radians(ra_deg), 2 * np.pi) / (np.pi / 2)
    ix = np.empty(len(z), np.int64)
    iy = np.empty(len(z), np.int64)
    face = np.empty(len(z), np.int64)

    eq = np.abs(z) <= 2.0 / 3.0
    t1 = nside * (0.5 + tt[eq])
    t2 = nside * (z[eq] * 0.75)
    jp = np.floor(t1 - t2).astype(np.int64)
    jm = np.floor(t1 + t2).astype(np.int64)
    shift = int(np.log2(nside))
    ifp, ifm = jp >> shift, jm >> shift
    face[eq] = np.where(ifp == ifm, (ifp & 3) + 4,
                        np.where(ifp < ifm, ifp & 3, (ifm & 3) + 8))
    ix[eq] = jm & (nside - 1)
    iy[eq] = nside - (jp & (nside - 1)) - 1

    po = ~eq
    ntt = np.floor(tt[po]).astype(np.int64).clip(0, 3)
    tp = tt[po] - ntt
    tmp = nside * np.sqrt(3.0 * (1.0 - np.abs(z[po])))
    jp = np.floor(tp * tmp).astype(np.int64).clip(0, nside - 1)
    jm = np.floor((1.0 - tp) * tmp).astype(np.int64).clip(0, nside - 1)
    north = z[po] > 0
    face[po] = np.where(north, ntt, ntt + 8)
    ix[po] = np.where(north, nside - jm - 1, jp)
    iy[po] = np.where(north, nside - jp - 1, jm)

    return (face * nside * nside + _spread_bits(ix)
            + (_spread_bits(iy) << 1))


def _read_healpix(path, field):
    with fits.open(path) as h:
        return np.asarray(h[1].data[field], np.float32).ravel()


_cache = {}


def delta_ebv(ra_deg, dec_deg, mapdir=MAPDIR):
    """E_SFD - E_CSFD [mag] at equatorial positions (galactic lookup)."""
    if "sfd" not in _cache:
        _cache["sfd"] = _read_healpix(
            f"{mapdir}/lambda_sfd_ebv.fits", "TEMPERATURE")
        _cache["csfd"] = _read_healpix(f"{mapdir}/csfd_ebv.fits", "T")
    el, b = eq2gal(ra_deg, dec_deg)
    sfd = _cache["sfd"][ang2pix_nest(512, el, b)]
    csfd = _cache["csfd"][ang2pix_ring(2048, el, b)]
    return sfd - csfd
