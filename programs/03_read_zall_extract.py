import numpy as np
import pandas as pd
import os
import sys

desilibdir = os.environ['DESILIB']
sys.path.append(desilibdir)
import desi_db

# -----------------------------------------------------------------------------
# Field footprints (RA/Dec in degrees)
# -----------------------------------------------------------------------------
COSMOS_RA_MIN, COSMOS_RA_MAX = 148.0, 152.0
COSMOS_DEC_MIN, COSMOS_DEC_MAX = 0.0, 4.0

SXDS_RA_MIN, SXDS_RA_MAX = 33.5, 38.0
SXDS_DEC_MIN, SXDS_DEC_MAX = -6.2, -3.0


def to_native(arr):
    """Return a numpy array in native byte order.

    FITS files (and therefore the zall arrays) are big-endian.  Pandas's
    Cython routines refuse to take rows from big-endian columns on a
    little-endian machine, raising:
        ValueError: Big-endian buffer not supported on little-endian compiler
    Converting to native byte order at construction time avoids that.
    """
    a = np.asarray(arr)
    if a.dtype.kind in 'biufc' and a.dtype.byteorder not in ('=', '|'):
        return a.astype(a.dtype.newbyteorder('='))
    return a


def extract_field(df, ra_min, ra_max, dec_min, dec_max, output_csv):
    """Filter df by an RA/Dec box and write the subset to output_csv."""
    sub = df[(df['ra']  >= ra_min)  & (df['ra']  <= ra_max) &
             (df['dec'] >= dec_min) & (df['dec'] <= dec_max)]
    print(f'{output_csv}: {len(sub)} rows '
          f'(RA {ra_min}-{ra_max}, Dec {dec_min}-{dec_max})')
    sub.to_csv(output_csv, index=False)
    return sub


# -----------------------------------------------------------------------------
# Read the full zall-pix catalog
# -----------------------------------------------------------------------------
zall = desi_db.zall()
zall.read_pix()

df = pd.DataFrame({'targetid':       to_native(zall.targetid),
                   'ra':             to_native(zall.target_ra),
                   'dec':            to_native(zall.target_dec),
                   'spectype':       to_native(zall.spectype),
                   'healpix':        to_native(zall.healpix),
                   'survey':         to_native(zall.survey),
                   'program':        to_native(zall.program),
                   'spgrpval':       to_native(zall.spgrpval),
                   'z':              to_native(zall.z),
                   'zerr':           to_native(zall.zerr),
                   'zwarn':          to_native(zall.zwarn),
                   'chi2':           to_native(zall.chi2),
                   'subtype':        to_native(zall.subtype),
                   'objtype':        to_native(zall.objtype),
                   'gaia_g':         to_native(zall.gaia_g),
                   'coadd_numexp':   to_native(zall.coadd_numexp),
                   'coadd_numnight': to_native(zall.coadd_numnight),
                   'coadd_numtile':  to_native(zall.coadd_numtile),
                   'coadd_exptime':  to_native(zall.coadd_exptime)})

df.to_csv('zall-pix-loa.csv', index=False)

# -----------------------------------------------------------------------------
# Extract COSMOS and SXDS subsets
# -----------------------------------------------------------------------------
extract_field(df, COSMOS_RA_MIN, COSMOS_RA_MAX,
                  COSMOS_DEC_MIN, COSMOS_DEC_MAX,
                  'zall-pix-loa-cosmos.csv')

extract_field(df, SXDS_RA_MIN, SXDS_RA_MAX,
                  SXDS_DEC_MIN, SXDS_DEC_MAX,
                  'zall-pix-loa-sxds.csv')

del df
