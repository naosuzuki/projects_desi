"""Native-pixel store for v2 (Nao, 2026-10-29): during the shared DESI coadd read (05_stack.py), keep the ORIGINAL spectra
of every quasar pass-Q row and every DESI main star (SPECTYPE STAR, a superset of the stellar candidates), so v2 rest-frame products (exact sub-pixel shift, one
flux-conserving resampling from the native pixels to a uniform log grid) are built from exdisk3 without reading exdisk2 again.

Stored as read: combine_cameras output (B/R/Z share the 0.8 A grid, ivar-weighted in the overlaps, no resampling),
no flux calibration, no extinction (CAL_A/CAL_B, MEAN_PSF_TO_FIBER_SPECFLUX, EBV come from the tables). A masked
pixel has IVAR = 0. One file per coadd: STORE/<coadd>.npz with
    TARGETID int64 [n], KIND int8 [n] (bits: 1 quasar pass-Q row, 2 stellar step-6 candidate, 4 SPECTYPE STAR), W0, DW (wavelength = W0 + DW * arange(npix)),
    FLUX float32 [n, npix], IVAR float32 [n, npix]   (DESI units 1e-17 erg/s/cm2/A)
A coadd without such a file has not been stored yet (native_todo() lists them for the re-read).
"""
import os
import numpy as np
import pandas as pd

STORE = '/Volumes/exdisk3/sdssdesi_native/desi'
QSO = os.path.expanduser('~/data/sdssdesi_products/quasar/qso_master_desi.parquet')
STARS = os.path.expanduser('~/data/sdssdesi_products/stars/step6_selection_desi_0.5x0.1_gspphot.parquet')


TARGETS = os.path.expanduser('~/data/sdssdesi_products/native_targets_desi.parquet')   # FILE, TARGETID, KIND (10.2 M rows)


def targets():
    """{coadd file: (TARGETID, KIND)} from TARGETS (built 2026-10-29 from qso_master_desi, the stellar selection and
    zall-pix-loa main dark+bright SPECTYPE STAR)."""
    t = pd.read_parquet(TARGETS)
    return {f: (g.TARGETID.to_numpy(np.int64), g.KIND.to_numpy(np.int8)) for f, g in t.groupby('FILE')}


def out(fname, store=STORE):
    return os.path.join(store, os.path.basename(fname).replace('.fits', '') + '.npz')


class NativeHook:
    def __init__(self, store=STORE):
        self.store = store; os.makedirs(store, exist_ok=True)
        self.by_file = targets()

    def process(self, fname, tid, w, flux, ivar):
        fname = os.path.basename(fname); p = out(fname, self.store)
        if fname not in self.by_file or os.path.exists(p):
            return 0
        want, kind = self.by_file[fname]
        pos = pd.Series(np.arange(len(tid)), index=np.asarray(tid, np.int64))
        pos = pos[~pos.index.duplicated()]
        keep = np.isin(want, pos.index); want, kind = want[keep], kind[keep]
        rows = pos.loc[want].to_numpy()
        dw = round(float(w[1] - w[0]), 6)
        tmp = p + '.tmp.npz'
        np.savez(tmp, TARGETID=want, KIND=kind, W0=float(w[0]), DW=dw,
                 FLUX=flux[rows].astype(np.float32), IVAR=ivar[rows].astype(np.float32))
        os.replace(tmp, p)
        return len(want)


def native_todo(files, store=STORE):
    """coadd files (basenames) that still need a read for the native store."""
    by = set(targets())
    return [f for f in files if os.path.basename(f) in by and not os.path.exists(out(f, store))]
