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
# v2 (Nao, 2026-10-29): one resampling from the native 0.8 A pixels straight onto the REST-frame lattice at each star's
# pipeline RV (programs_templates/restframe.desi_rest; sub-pixel), stored as RID0 = SRID_MIN, FRAME = 'rest'; step6 --v2
# reads it without a further shift. v1 (STORE/ARCHIVE above) stays as it is: observed log grid, shifted in step 6.
STORE_V2 = os.path.expanduser('~/data/sdssdesi_stars/specstore_desi_v2')
ARCHIVE_V2 = '/Volumes/exdisk3/sdssdesi_star_spectra/desi_v2'
GAUGE = os.path.expanduser('~/data/sdssdesi_status/step_star_stack.json')
GAUGE_EVERY = 25          # coadds per hook instance between gauge updates (the count is the store itself, so shards agree)


class StarHook:
    """writes STORE/<coadd>.npz (the step-6 stars, internal disk) and ARCHIVE/<coadd>.npz (every candidate, exdisk3)."""
    def __init__(self, sel=SEL, store=None, sel300=SEL300, archive=None, version=2):
        import starlib
        self.passB = starlib.passB; self.version = version
        if store is None:
            store = STORE_V2 if version == 2 else STORE
        if archive is None:
            archive = ARCHIVE_V2 if version == 2 else ARCHIVE
        if version == 2:
            sys.path.insert(0, os.path.join(HERE, '..', 'programs_templates')); import restframe
            self.restframe = restframe; self.srid_min, self.npix = starlib.SRID_MIN, len(starlib.SWAVE)
            assert abs(restframe.DLOG - starlib.DLOG) < 1e-12, 'restframe and starlib lattices differ'
        s = pd.read_parquet(sel, columns=['FILE', 'TARGETID', 'Z'])
        self.by_file_all = {f: g.TARGETID.to_numpy(np.int64) for f, g in s.groupby('FILE')}
        self.z = dict(zip(s.TARGETID.to_numpy(np.int64), s.Z.to_numpy(float)))
        s = pd.read_parquet(sel300, columns=['FILE', 'TARGETID'])
        self.by_file = {f: g.TARGETID.to_numpy(np.int64) for f, g in s.groupby('FILE')}
        self.store = store; os.makedirs(store, exist_ok=True)
        self.archive = archive
        if archive:
            os.makedirs(archive, exist_ok=True)
        self.nfile = 0; self.nstar = 0
        sp = os.path.join(os.path.dirname(store), 'started_desi.txt' if version == 1 else 'started_desi_v2.txt')
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
            if self.version == 2:            # rest frame at the pipeline RV, one resampling (|z| clipped at 0.004 as in step 6)
                zz = np.clip(np.nan_to_num(np.array([self.z.get(int(t), 0.0) for t in want]), nan=0.0), -0.004, 0.004)
                rf, ri = self.restframe.desi_rest(w, flux[rows], ivar[rows], mask[rows], zz, self.srid_min, self.npix)
                rid0 = self.srid_min
            else:
                rid, rf, ri = self.passB.rebin(w, flux[rows].astype(float), ivar[rows].astype(float), mask[rows])
                rid0 = int(rid[0]); zz = np.zeros(len(want))
            rf = rf.astype(np.float32); ri = ri.astype(np.float32)
        else:
            rf = ri = np.zeros((0, 1), np.float32); rid0 = 0; zz = np.zeros(0)
        frame = 'rest' if self.version == 2 else 'observed'
        def save(d, keep):
            p = self.out(fname, d); tmp = p + '.tmp.npz'
            np.savez(tmp, TARGETID=want[keep], RID0=rid0, FLUX=rf[keep], IVAR=ri[keep], Z=zz[keep], FRAME=frame); os.replace(tmp, p)
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
        write_gauge(self.store, len(self.by_file), self.started)


QSO_STATUS = os.path.expanduser('~/data/sdssdesi_status')
PRIORITY = os.path.expanduser('~/data/sdssdesi_products/stars/star_priority_files.txt')


def eta_epoch(store=STORE, nshard=3, tier1=11288):
    """when the last pending step-6 coadd will have been read by the quasar shards: each shard's read order as in
    programs_quasars/05_stack.py (phase 0 = star tier 1, 1 = other dark, 2 = bright; star rank within a phase), the position of
    its last pending star coadd, the shard's rate since its restart. (The recent arrival rate is useless: almost no star coadds
    arrive while the shards read the dark quasar coadds.) None when it cannot be estimated."""
    now = time.time()
    rank = {l.strip(): i for i, l in enumerate(open(PRIORITY)) if l.strip()}
    have = {e.name[:-4] + '.fits' for e in os.scandir(store) if e.name.endswith('.npz') and '.tmp' not in e.name}
    need = set(pd.read_parquet(SEL300, columns=['FILE']).FILE) - have
    phase = lambda f: 0 if rank.get(f, 10**9) < tier1 else (1 if 'dark' in f else 2)
    worst = 0.0
    for i in range(nshard):
        try:
            fs = sorted({os.path.basename(l.strip()) for l in open(f'{QSO_STATUS}/qso_desi_files_shard{i}of{nshard}.txt') if l.strip()},
                        key=lambda f: (phase(f), rank.get(f, 10**9), f))
            gone = set()                # MISSING in the shard's current run only
            for l in open(f'{QSO_STATUS}/qso_stack_desi{i}.log', errors='replace'):
                if 'files already done' in l: gone = set()
                elif l.startswith('MISSING '): gone.add(l.split()[1])
            pos = [k for k, f in enumerate(fs) if f in need and f not in gone]
            if not pos:
                continue
            j = json.load(open(f'{QSO_STATUS}/step_qso_stack_desi_shard{i}of{nshard}.json'))
            if j.get('status') != 'running':
                return None
            d0 = [int(l.split(';')[1].split()[0]) for l in open(f'{QSO_STATUS}/qso_stack_desi{i}.log', errors='replace') if 'files already done' in l][-1]
            rate = (j['done'] - d0) / max(now - j['started'], 1.0)
            if rate <= 0:
                return None
            worst = max(worst, (pos[-1] + 1 - j['done']) / rate)
        except Exception:
            return None
    return now + worst


def write_gauge(store=STORE, tot=None, started=None):
    """monitor gauge for step 6: coadds whose star spectra are in the store, out of all coadds with step-6 stars; ETA from the
    shards' read order (eta_epoch), not from the recent rate."""
    try:
        if tot is None:
            tot = pd.read_parquet(SEL300, columns=['FILE']).FILE.nunique()
        if started is None:
            started = float(open(os.path.join(os.path.dirname(store), 'started_desi_v2.txt' if store == STORE_V2 else 'started_desi.txt')).read().strip())
        n = sum(1 for e in os.scandir(store) if e.name.endswith('.npz') and '.tmp' not in e.name)
        if store == STORE_V2:            # v2: the shared walk reads uniformly, the monitor's rate is fine
            n = sum(1 for e in os.scandir(store) if e.name.endswith('.npz') and '.tmp' not in e.name)
            note = f'v2 (rest frame, one resampling at the pipeline RV): star spectra stored for {n:,} / {tot:,} coadds, inside the shared DESI walk'
            tmp = f'{GAUGE}.{os.getpid()}.tmp'
            json.dump({'done': n, 'total': tot, 'status': 'running', 'note': note, 'started': started}, open(tmp, 'w')); os.replace(tmp, GAUGE)
            return
        eta = eta_epoch(store)
        note = (f'DESI stacks done 10-28 07:22 (898 stacks, median 230 stars); filling the rest of the 300-per-stack set for an '
                f'optional restack: {n:,} / {tot:,} coadds stored, inside the quasar-stack read (bright coadds come last in its order)'
                + (f'; last one ~{time.strftime("%H:%M", time.localtime(eta))}' if eta else ''))
        tmp = f'{GAUGE}.{os.getpid()}.tmp'
        json.dump({'done': n, 'total': tot, 'status': 'running', 'note': note, 'started': started, 'eta_epoch': eta}, open(tmp, 'w')); os.replace(tmp, GAUGE)
    except Exception:
        pass


if __name__ == '__main__':      # python3 star_hook.py --gauge [SECONDS]: refresh the gauge (once, or every SECONDS until the store is full)
    if '--gauge' in sys.argv:
        every = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 0
        while True:
            write_gauge(); j = json.load(open(GAUGE))
            if not every or j['done'] >= j['total']:
                break
            time.sleep(every)
