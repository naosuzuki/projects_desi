"""Star hook for a DESI coadd read that someone else already does (the quasar stack, programs_quasars/05_stack.py), so step 6
needs no exdisk2 pass of its own (Nao: one read per coadd for all classes).

The selection (which stars of which coadd: every star of the step-5 cells with S/N > 5, before the per-stack cap) comes from
`step6_stack.py --select-only`: ~/data/sdssdesi_products/stars/
step6_selection_<tag>.parquet (FILE, TARGETID, ...). For a coadd with selected stars, process() rebins those rows onto the log
grid exactly as step6_stack.iter_desi does and writes STORE/<FILE>.npz (TARGETID, RID0, FLUX, IVAR; float32), atomically.
`step6_stack.py --from-store` stacks from STORE and reads from exdisk2 only the coadds that are not there.

    SH = star_hook.StarHook()          # None-safe: StarHook.load() returns None when the selection does not exist
    ... after w, flux, ivar, mask = combine_cameras(f):  SH.process(fname, tid, w, flux, ivar, mask)
"""
import os, sys, json, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
TAG = 'desi_0.5x0.1_gspphot'
SEL = os.path.expanduser(f'~/data/sdssdesi_products/stars/step6_selection_{TAG}.parquet')
ARCHIVE = '/Volumes/exdisk3/sdssdesi_star_spectra/desi'      # exdisk3 (Nao, 2026-10-27): every candidate star, ~35 GB (slow to read back)
STORE = os.path.expanduser('~/data/sdssdesi_stars/specstore_desi')   # internal: the stars step 6 stacks (300 per stack), ~7 GB, fast
SEL300 = SEL.replace('.parquet', '_max300.parquet')
GAUGE = os.path.expanduser('~/data/sdssdesi_status/step_star_stack.json')
GAUGE_EVERY = 25          # coadds per hook instance between gauge updates (the count is the store itself, so shards agree)


class StarHook:
    """writes STORE/<coadd>.npz (the step-6 stars, internal disk) and ARCHIVE/<coadd>.npz (every candidate, exdisk3)."""
    def __init__(self, sel=SEL, store=STORE, sel300=SEL300, archive=ARCHIVE):
        import starlib
        self.passB = starlib.passB
        s = pd.read_parquet(sel, columns=['FILE', 'TARGETID'])
        self.by_file_all = {f: g.TARGETID.to_numpy(np.int64) for f, g in s.groupby('FILE')}
        s = pd.read_parquet(sel300, columns=['FILE', 'TARGETID'])
        self.by_file = {f: g.TARGETID.to_numpy(np.int64) for f, g in s.groupby('FILE')}
        self.store = store; os.makedirs(store, exist_ok=True)
        self.archive = archive
        if archive:
            os.makedirs(archive, exist_ok=True)
        self.nfile = 0; self.nstar = 0
        sp = os.path.join(os.path.dirname(store), 'started_desi.txt')
        if not os.path.exists(sp):
            open(sp, 'w').write(str(time.time()))
        self.started = float(open(sp).read().strip() or time.time())
        self.gauge()

    @classmethod
    def load(cls, sel=SEL, store=STORE):
        return cls(sel, store) if os.path.exists(sel) else None

    def out(self, fname, d=None):
        return os.path.join(d or self.store, os.path.basename(fname).replace('.fits', '') + '.npz')

    def needs(self, fname):
        f = os.path.basename(fname)
        return (f in self.by_file and not os.path.exists(self.out(f))) or \
               (bool(self.archive) and f in self.by_file_all and not os.path.exists(self.out(f, self.archive)))

    def process(self, fname, tid, w, flux, ivar, mask):
        """tid: FIBERMAP TARGETID of all rows; w/flux/ivar/mask: combine_cameras output for all rows."""
        fname = os.path.basename(fname)
        if not self.needs(fname):
            return 0
        pos = pd.Series(np.arange(len(tid)), index=np.asarray(tid, np.int64))
        want = self.by_file_all.get(fname, self.by_file.get(fname)); want = want[np.isin(want, pos.index)]
        if len(want):
            rows = pos.loc[want].to_numpy()
            rows = rows if rows.ndim == 1 else rows.ravel()
            rid, rf, ri = self.passB.rebin(w, flux[rows].astype(float), ivar[rows].astype(float), mask[rows])
            rid0 = int(rid[0]); rf = rf.astype(np.float32); ri = ri.astype(np.float32)
        else:
            rf = ri = np.zeros((0, 1), np.float32); rid0 = 0
        def save(d, keep):
            p = self.out(fname, d); tmp = p + '.tmp.npz'
            np.savez(tmp, TARGETID=want[keep], RID0=rid0, FLUX=rf[keep], IVAR=ri[keep]); os.replace(tmp, p)
        if fname in self.by_file and not os.path.exists(self.out(fname)):
            save(self.store, np.isin(want, self.by_file[fname]))
        if self.archive and fname in self.by_file_all and not os.path.exists(self.out(fname, self.archive)):
            save(self.archive, np.ones(len(want), bool))
        want = want[np.isin(want, self.by_file.get(fname, want[:0]))]
        self.nfile += 1; self.nstar += len(want)
        if self.nfile % GAUGE_EVERY == 0:
            self.gauge()
        return len(want)

    def gauge(self):
        """monitor gauge for step 6: coadds whose star spectra are in the store, out of all coadds with selected stars."""
        try:
            n = sum(1 for e in os.scandir(self.store) if e.name.endswith('.npz') and '.tmp' not in e.name); tot = len(self.by_file)
            note = (f'DESI: step-6 star spectra stored for {n:,} / {tot:,} coadds (inside the quasar-stack read + fill_store.py; '
                    f'no extra pass); step 6 starts when all are stored or at 07:15')
            tmp = f'{GAUGE}.{os.getpid()}.tmp'
            json.dump({'done': n, 'total': tot, 'status': 'running', 'note': note, 'started': self.started}, open(tmp, 'w')); os.replace(tmp, GAUGE)
        except Exception:
            pass
