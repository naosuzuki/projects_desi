"""Galaxy stacking grid: bin rules shared by 07_assign_cells.py and 08_passC_stack.py (ms.tex Table galgrid).

Axes: z (0.05, 0-1.6) x M (1 mag, band-switched absolute magnitude from the whole-galaxy flux:
M_R z<0.35, M_G 0.35-0.72, M_U 0.72-1.45, M_NUV 1.45-1.6) x color (0.25 mag, dereddened LS total fluxes,
g-r for z<0.5, r-z above) x line EW (log2 bins; Halpha z<0.46, [OII] above) x Dn4000 (0.25 bins).
"""
import numpy as np

DZ = 0.05
ZMAX = 1.6
M_MIN, M_MAX = -24, -16            # MB index range (floor of M)
C_MIN, C_MAX = 0, 7                # CB index range (floor of color/0.25)
EW_EDGES = np.array([1, 2, 4, 8, 16, 32, 64])      # 8 bins: <1 ... >64 A (bin 0 also holds EW<=0 / no line)
DN_EDGES = np.array([1.15, 1.40, 1.65, 1.90])      # 5 bins
LS_EXT = {'G': 3.214, 'R': 2.165, 'Z': 1.211}
Z_LINE_SWITCH = 0.46               # Halpha below, [OII] above
KEY = ['ZB', 'MB', 'CB', 'EB', 'DB']


def mband(z):
    z = np.asarray(z)
    return np.where(z < 0.35, 'R', np.where(z < 0.72, 'G', np.where(z < 1.45, 'U', 'NUV')))


def pick_mabs(z, mabs_r, mabs_g, mabs_u, mabs_nuv):
    z = np.asarray(z)
    return np.where(z < 0.35, mabs_r, np.where(z < 0.72, mabs_g, np.where(z < 1.45, mabs_u, mabs_nuv)))


def color(z, flux_g, flux_r, flux_z, ebv):
    ebv = np.nan_to_num(np.asarray(ebv, float), nan=0.0)
    with np.errstate(all='ignore'):
        mg = 22.5 - 2.5 * np.log10(flux_g) - LS_EXT['G'] * ebv
        mr = 22.5 - 2.5 * np.log10(flux_r) - LS_EXT['R'] * ebv
        mz = 22.5 - 2.5 * np.log10(flux_z) - LS_EXT['Z'] * ebv
    c = np.where(np.asarray(z) < 0.5, mg - mr, mr - mz)
    return np.where(np.isfinite(c), c, np.nan)


def line_ew(z, ew_halpha, ew_oii):
    return np.where(np.asarray(z) < Z_LINE_SWITCH, ew_halpha, ew_oii)


def indices(z, mabs, col, ew, dn):
    """Return dict of int index arrays (ZB, MB, CB, EB, DB) and a boolean 'valid' (inside the grid, finite)."""
    z = np.asarray(z, float); mabs = np.asarray(mabs, float); col = np.asarray(col, float)
    ok = np.isfinite(z) & np.isfinite(mabs) & np.isfinite(col) & (z > 0) & (z < ZMAX)
    zb = np.floor(np.where(ok, z, 0) / DZ).astype(int)
    mb = np.floor(np.where(ok, mabs, 0)).astype(int)
    cb = np.floor(np.where(ok, col, 0) / 0.25).astype(int)
    eb = np.searchsorted(EW_EDGES, np.nan_to_num(np.asarray(ew, float), nan=-99.0))
    db = np.searchsorted(DN_EDGES, np.nan_to_num(np.asarray(dn, float), nan=-9.0))
    ok &= (mb >= M_MIN) & (mb <= M_MAX) & (cb >= C_MIN) & (cb <= C_MAX)
    return dict(ZB=zb, MB=mb, CB=cb, EB=eb, DB=db, valid=ok)


def cell_code(ix):
    """Pack the five indices into one int64 key (unique for the grid ranges used)."""
    return (((ix['ZB'].astype(np.int64) * 100 + (ix['MB'] + 50)) * 100 + (ix['CB'] + 50)) * 10 + ix['EB']) * 10 + ix['DB']
