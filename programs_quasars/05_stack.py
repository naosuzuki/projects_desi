#!/usr/bin/env python3
"""Quasar stacks (handover sec. 8.24/8.26): every training quasar (04_cells.py TRAIN) is shifted to the rest frame with
Z_SYS, normalized by the mean of its rest spectrum in the cell's normalization window (5100-5130 A for z<0.8,
3000-3050 A for 0.8<=z<2.0, 1275-1290 A for z>=2.0, chosen by the cell's z slice, so one window per stack), and added
to its stack cell. One read per file (DESI coadd on exdisk2 / SDSS spPlate on the internal disk).
Per cell and rest pixel (qsolib.QWAVE, 600-9900 A, dlog 1e-4) the accumulators hold
  s1 = sum f, s2 = sum f^2, n = count (unweighted mean and scatter), wf = sum ivar f, w = sum ivar (ivar-weighted mean).
Normalization: NORM = ivar-weighted mean of the rest flux (F_lam (1+z), cgs) in the window; if < 50 % of the window is
good, the power-law continuum value there (CONT_LOGF3000, ALPHA_LAM; NORM_SRC 1); NORM_SRC 2 = neither (skipped).
The normalized spectra (float16, cropped to coverage) go to --store/<survey>/<file>.npz for medians and EM-PCA.
Checkpoints every 100 files (accumulators + done list, resumable).
Usage:
  python3 05_stack.py --survey sdss --shard 0/4          # spPlates, internal disk
  python3 05_stack.py --survey desi --shard 0/2          # DESI coadds, exdisk2 (read-only)
  python3 05_stack.py --merge                            # all shard accumulators -> qso_stacks.npz + qso_stacks_cells.csv
"""
import os, sys, glob, time, json, argparse, signal, gc, warnings
import numpy as np
import pandas as pd
import fitsio
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qsolib as Q
import importlib
PQS = importlib.import_module('01_passQ_sdss')

HOME = os.path.expanduser('~')
PROD = f'{HOME}/data/sdssdesi_products/quasar'
ACCDIR = f'{PROD}/stack_acc'
STATUS_DIR = f'{HOME}/data/sdssdesi_status'
STORE = '/Volumes/exdisk3/sdssdesi_qso_restframe'
# v2 (--v2, Nao 2026-10-29): exact sub-pixel rest frame (qsolib exact=True); DESI from the native-pixel store written by the
# shared walk (native_hook.py), SDSS from the spPlates with the exact COEFF0; outputs PROD/v2, v1 untouched.
V2 = False
OUT = PROD
SNAME = 'step_qso_stack'
CHECKPOINT_EVERY = 100
MIN_WIN_GOOD = 0.5
DESI_COLS = ['TARGETID', 'PROGRAM', 'FILE', 'Z', 'Z_SYS', 'CAL_A', 'CAL_B', 'EBV', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'CONT_LOGF3000', 'ALPHA_LAM']
SDSS_COLS = ['PLATE', 'MJD', 'FIBERID', 'SURVEY', 'Z', 'Z_SYS', 'CAL_A', 'EBV', 'CONT_LOGF3000', 'ALPHA_LAM']
STOP = {'flag': False}
signal.signal(signal.SIGTERM, lambda s, f: STOP.__setitem__('flag', True))


def members():
    """Training quasars with their stack cell (string label -> integer index, sorted) and the cell's window."""
    c = pd.read_parquet(f'{PROD}/qso_cells_all.parquet', columns=['PLATE', 'MJD', 'FIBERID', 'TARGETID', 'PROGRAM', 'SURVEY', 'TRAIN',
                                                                   'STACK_CELL', 'STACK_ZBIN', 'STACK_N'])
    c = c[c.TRAIN & (c.STACK_CELL != '')].copy()
    cells = np.array(sorted(c.STACK_CELL.unique()))
    c['CELL_I'] = np.searchsorted(cells, c.STACK_CELL.to_numpy()).astype(np.int32)
    zc = 0.5 * (Q_ZEDGES[c.STACK_ZBIN.to_numpy()] + Q_ZEDGES[c.STACK_ZBIN.to_numpy() + 1])
    lo = np.full(len(c), np.nan); hi = np.full(len(c), np.nan)
    for zlo, zhi, a, b in Q.NORM_WINDOWS:
        m = (zc >= zlo) & (zc < zhi); lo[m] = a; hi[m] = b
    c['WIN_LO'], c['WIN_HI'] = lo, hi
    return c, cells


ZE = importlib.import_module('04_cells')
Q_ZEDGES = np.asarray(ZE.ZEDGES, float)


def normalize(rflux, rivar, win_lo, win_hi, logf3000, alpha):
    """Window mean (ivar-weighted) of each row; power-law fallback. Returns norm [n], src [n] (0 window, 1 fit, 2 none)."""
    n = len(rflux); norm = np.full(n, np.nan); src = np.full(n, 2, np.int8)
    for lo, hi in set(zip(win_lo, win_hi)):
        rows = np.where((win_lo == lo) & (win_hi == hi))[0]
        sel = (Q.QWAVE >= lo) & (Q.QWAVE <= hi)
        f = rflux[np.ix_(rows, np.where(sel)[0])]; iv = rivar[np.ix_(rows, np.where(sel)[0])]
        good = iv > 0
        with np.errstate(all='ignore'):
            m = (f * iv).sum(axis=1) / np.where(good.any(axis=1), (iv * good).sum(axis=1), np.nan)
        ok = (good.mean(axis=1) >= MIN_WIN_GOOD) & np.isfinite(m) & (m > 0)
        norm[rows[ok]] = m[ok]; src[rows[ok]] = 0
        fit = 10**(logf3000[rows] + alpha[rows] * np.log10(0.5 * (lo + hi) / 3000.0))
        use = ~ok & np.isfinite(fit) & (fit > 0)
        norm[rows[use]] = fit[use]; src[rows[use]] = 1
    return norm, src


class Acc:
    def __init__(self, path, ncell):
        self.path = path
        if os.path.exists(path):
            z = np.load(path, allow_pickle=False)
            self.a = {k: z[k] for k in ('s1', 's2', 'n', 'wf', 'w')}; self.done = set(z['done'].tolist())
            self.nobj = int(z['nobj']); self.nsrc = z['nsrc']
        else:
            shp = (ncell, Q.QNREST)
            self.a = dict(s1=np.zeros(shp, np.float32), s2=np.zeros(shp, np.float32), n=np.zeros(shp, np.int32),
                          wf=np.zeros(shp, np.float32), w=np.zeros(shp, np.float32))
            self.done = set(); self.nobj = 0; self.nsrc = np.zeros(3, np.int64)

    def add(self, cell_i, fn, ivn):
        good = ivn > 0
        for k in np.unique(cell_i):
            r = cell_i == k; f = np.where(good[r], fn[r], 0.0); g = good[r]; iv = np.where(g, ivn[r], 0.0)
            self.a['s1'][k] += f.sum(0); self.a['s2'][k] += (f * f).sum(0); self.a['n'][k] += g.sum(0)
            self.a['wf'][k] += (iv * f).sum(0); self.a['w'][k] += iv.sum(0)

    def save(self):
        tmp = self.path + '.tmp.npz'
        np.savez(tmp, **self.a, done=np.array(sorted(self.done)), nobj=self.nobj, nsrc=self.nsrc)
        os.replace(tmp, self.path)


def store_file(store, survey, base, ids, cell_i, norm, src, fn, ivn):
    """Normalized spectra, float16, cropped to each row's coverage (START = first QWAVE pixel, OFFSET into FLUX/IVAR)."""
    d = os.path.join(store, survey); os.makedirs(d, exist_ok=True)
    starts, flx, ivs = [], [], []
    for i in range(len(fn)):
        g = np.where(ivn[i] > 0)[0]
        a0, a1 = (g[0], g[-1] + 1) if len(g) else (0, 0)
        starts.append(a0); flx.append(fn[i, a0:a1].astype(np.float16)); ivs.append(ivn[i, a0:a1].astype(np.float16))
    lens = np.array([len(x) for x in flx]); off = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
    out = os.path.join(d, base + '.npz'); tmp = out + '.tmp.npz'
    np.savez(tmp, ID=ids, CELL_I=cell_i, NORM=norm, NORM_SRC=src, START=np.array(starts, np.int32), OFFSET=off,
             FLUX=np.concatenate(flx) if flx else np.zeros(0, np.float16), IVAR=np.concatenate(ivs) if ivs else np.zeros(0, np.float16),
             QRID_MIN=Q.QRID_MIN)
    os.replace(tmp, out)


def finish_rows(acc, store, survey, base, ids, cell_i, rflux, rivar, sub):
    norm, src = normalize(rflux, rivar, sub.WIN_LO.to_numpy(), sub.WIN_HI.to_numpy(),
                          sub.CONT_LOGF3000.to_numpy(dtype=float), sub.ALPHA_LAM.to_numpy(dtype=float))
    ok = src < 2
    fn = np.where(ok[:, None], rflux / np.where(ok, norm, 1.0)[:, None], 0.0).astype(np.float32)
    ivn = np.where(ok[:, None], rivar * np.where(ok, norm, 1.0)[:, None]**2, 0.0).astype(np.float32)
    acc.add(cell_i[ok], fn[ok], ivn[ok])
    acc.nobj += int(ok.sum()); acc.nsrc += np.bincount(src, minlength=3)
    if store:
        store_file(store, survey, base, ids, cell_i, norm, src, fn, ivn)


NoMaskFITS = Q.NoMaskFITS      # moved to qsolib (shared with the walk 02_passQ_desi.py)


def star_hook():
    """Stellar thread's hook (programs_stars/star_hook.py): the same coadd read also feeds its step-6 star spectra
    (one read per coadd for all tracks). None when the module or its selection file is missing."""
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'programs_stars'))
        import star_hook as SHM
        return SHM.StarHook()
    except Exception as e:
        print('star hook not available:', repr(e), flush=True)
        return None


def run_desi(files, mem, acc, store, t0, status):
    global native_hook
    import native_hook
    SH = NH = None
    if not V2:                                                    # hooks only on an exdisk2 read
        SH = star_hook()
        try:
            NH = native_hook.NativeHook()                         # v2: native-pixel spectra of quasars + stars on exdisk3
        except Exception as e:
            print('native hook not available:', repr(e), flush=True)
    print('star hook:', 'on' if SH else 'off', '| native hook:', 'on' if NH else 'off', flush=True)
    for k, (fname, sub) in enumerate(files):
        if STOP['flag']:
            break
        if fname in acc.done:
            continue
        if PHASE.get(fname, 0) == 2 and not os.path.exists(acc.path + '.dark_done'):
            acc.save(); open(acc.path + '.dark_done', 'w').write(time.strftime('%F %T'))   # provisional merge point
        if V2:                                                    # native-pixel store on exdisk3 (the shared v2 walk)
            npz = native_hook.out(fname)
            if not os.path.exists(npz):
                print('NOT STORED', fname, flush=True); continue
            z = np.load(npz)
            tid = z['TARGETID']; flux = z['FLUX']; ivar = z['IVAR']; mask = (ivar <= 0).astype(np.int8)
            w = float(z['W0']) + float(z['DW']) * np.arange(flux.shape[1]); del z
        else:
            path = Q.passB.local_path(PATHS[fname]) if fname in PATHS else None   # 969 overflow coadds live on exdisk1
            if path is None:
                print('MISSING', fname, flush=True); continue
            f = fitsio.FITS(path)
            tid = f['FIBERMAP'].read(columns=['TARGETID'])['TARGETID']
            w, flux, ivar, mask, _ = Q.passB.combine_cameras(NoMaskFITS(f)); f.close()
        if SH:
            try:
                SH.process(fname, tid, w, flux, ivar, mask)
            except Exception as e:                                # a star-hook failure must not stop the quasar stacks
                print('star hook error', fname, repr(e), flush=True)
        if NH:
            try:
                NH.process(fname, tid, w, flux, ivar)
            except Exception as e:                                # nor a native-store failure
                print('native hook error', fname, repr(e), flush=True)
        pos = pd.Series(np.arange(len(tid)), index=tid.astype(np.int64))
        sub = sub[sub.TARGETID.isin(pos.index)]
        rows = pos.loc[sub.TARGETID.to_numpy()].to_numpy()
        df = sub.copy(); df['Z'] = df['Z_SYS']                  # rest frame with the systemic redshift
        for i0 in range(0, len(df), 2000):
            s = slice(i0, i0 + 2000); r = rows[s]
            rflux, rivar = Q.quasar_rest(df.iloc[s], w, flux[r], ivar[r], mask[r], exact=V2)
            finish_rows(acc, store, 'desi', fname.replace('.fits', '') + (f'_{i0}' if i0 else ''), df.TARGETID.to_numpy()[s],
                        df.CELL_I.to_numpy()[s], rflux, rivar, df.iloc[s])
        del w, flux, ivar, mask; gc.collect()
        acc.done.add(fname)
        status(fname)


def run_sdss(files, mem, acc, store, t0, status):
    for k, (fname, sub) in enumerate(files):
        if STOP['flag']:
            break
        if fname in acc.done:
            continue
        path = PATHS.get(fname)
        if path is None:
            print('MISSING', fname, flush=True); continue
        F = fitsio.FITS(path)
        h = F[0].read_header(); rid0 = int(round(float(h['COEFF0']) / Q.DLOG))
        fib = sub.FIBERID.to_numpy().astype(int) - 1
        flux = F[0].read()[fib].astype(float); ivar = F[1].read()[fib].astype(float); mask = F[2].read()[fib]
        F.close()
        rflux, rivar = Q.rest_from_loggrid(rid0, flux, ivar, mask, sub.Z_SYS.to_numpy(dtype=float),
                                           np.nan_to_num(sub.CAL_A.to_numpy(dtype=float), nan=1.0), sub.EBV.to_numpy(dtype=float),
                                           coeff0=float(h['COEFF0']) if V2 else None)
        ids = (sub.PLATE.to_numpy().astype(np.int64) * 100000 + sub.MJD.to_numpy().astype(np.int64)) * 10000 + sub.FIBERID.to_numpy().astype(np.int64)
        finish_rows(acc, store, sub.SURVEY.iloc[0], fname, ids, sub.CELL_I.to_numpy(), rflux, rivar, sub)
        del flux, ivar, mask, rflux, rivar
        acc.done.add(fname)
        status(fname)


PATHS = {}
PHASE = {}
STAR_PRIORITY = os.path.expanduser('~/data/sdssdesi_products/stars/star_priority_files.txt')
STAR_TIER1 = 11288


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', choices=['sdss', 'desi'])
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--store', default=None, help="normalized rest spectra; '' to skip (default v1 STORE, v2 STORE + '_v2')")
    ap.add_argument('--nfiles', type=int, default=0, help='test: only the first N files of this shard')
    ap.add_argument('--accdir', default=None, help='default v1 PROD/stack_acc, v2 PROD/v2/stack_acc')
    ap.add_argument('--v2', action='store_true', help='v2: exact sub-pixel rest frame, DESI from the native store')
    ap.add_argument('--merge', action='store_true')
    ap.add_argument('--out', default='qso_stacks', help='--merge output name (qso_stacks_dark for the provisional merge)')
    a = ap.parse_args()
    global V2, OUT, SNAME
    if a.v2:
        V2, OUT, SNAME = True, f'{PROD}/v2', 'step_qso_stack_v2'
        os.makedirs(OUT, exist_ok=True)
    if a.accdir is None:
        a.accdir = f'{OUT}/stack_acc'
    if a.store is None:
        a.store = STORE + ('_v2' if V2 else '')
    if a.merge:
        return merge(a.accdir, a.out)
    si, sn = map(int, a.shard.split('/'))
    t0 = time.time()
    mem, cells = members()
    os.makedirs(a.accdir, exist_ok=True)
    pd.DataFrame({'CELL_I': np.arange(len(cells)), 'STACK_CELL': cells}).to_csv(f'{a.accdir}/cells_index.csv', index=False)
    if a.survey == 'desi':
        m = pd.read_parquet(f'{PROD}/qso_master_desi.parquet', columns=DESI_COLS)
        sub = mem[mem.SURVEY == 'desi'].merge(m, on=['TARGETID', 'PROGRAM'], how='inner')
        for l in open(Q.passB.MANIFEST):
            p = l.split()[1]; PATHS[os.path.basename(p)] = p
    else:
        m = pd.read_parquet(f'{PROD}/qso_master_sdss.parquet', columns=SDSS_COLS)
        sub = mem[mem.SURVEY.isin(['dr17', 'dr8'])].drop(columns=['TARGETID', 'PROGRAM']).merge(m, on=['PLATE', 'MJD', 'FIBERID', 'SURVEY'], how='inner')
        sub['FILE'] = ['spPlate-%04d-%05d' % (p, j) for p, j in zip(sub.PLATE.to_numpy(), sub.MJD.to_numpy())]
        for sv in ('dr17', 'dr8'):
            for p in PQS.plate_files(sv):
                PATHS[os.path.basename(p).replace('.fits', '')] = p
    groups = sorted(sub.groupby('FILE'), key=lambda x: x[0])
    files = [g for i, g in enumerate(groups) if i % sn == si]      # shard membership fixed by the sorted index
    if a.survey == 'desi':
        # order (2026-10-28 01:10, shared with the Stellar thread): phase 0 = coadds holding the stars' 100-per-stack
        # set (first STAR_TIER1 lines of star_priority_files.txt), phase 1 = the remaining dark coadds (the bulk of the
        # quasars), phase 2 = the rest (bright), the stars' priority rank within each phase
        rank = {}
        if os.path.exists(STAR_PRIORITY):
            rank = {l.strip(): i for i, l in enumerate(open(STAR_PRIORITY)) if l.strip()}
        def phase(fn):
            return 0 if rank.get(fn, 10**9) < STAR_TIER1 else (1 if 'dark' in fn else 2)
        files.sort(key=lambda g: (phase(g[0]), rank.get(g[0], 10**9), g[0]))
        PHASE.update({g[0]: phase(g[0]) for g in files})
    if a.nfiles:
        files = files[:a.nfiles]
    acc = Acc(f'{a.accdir}/acc_{a.survey}_shard{si}of{sn}.npz', len(cells))
    if a.survey == 'desi':      # complete list of this shard's coadds (pending or done) and the done ones, for the Stellar thread
        if not V2: open(f'{STATUS_DIR}/qso_desi_files_shard{si}of{sn}.txt', 'w').write(''.join(g[0] + '\n' for g in files))
        if not V2: open(f'{STATUS_DIR}/qso_desi_done_shard{si}of{sn}.txt', 'w').write(''.join(f + '\n' for f in sorted(acc.done)))
    ntot = len(files); nsince = [0]
    stp = f'{STATUS_DIR}/{SNAME}_{a.survey}_shard{si}of{sn}.json'
    print(f'{a.survey} shard {si}/{sn}: {ntot} files, {sum(len(g) for _, g in files):,} quasars, {len(cells)} cells; '
          f'{len(acc.done)} files already done', flush=True)

    def status(fname):
        nsince[0] += 1
        nd = len(acc.done)
        if nd % 20 == 0:
            print(f'{nd}/{ntot} {fname} nobj={acc.nobj:,} src={acc.nsrc.tolist()} {time.time()-t0:.0f}s', flush=True)
        if nsince[0] >= CHECKPOINT_EVERY:
            acc.save(); nsince[0] = 0
        json.dump(dict(done=nd, total=ntot, status='running', started=t0, note=f'{acc.nobj:,} quasars'), open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
        agg_status()

    (run_desi if a.survey == 'desi' else run_sdss)(files, mem, acc, a.store or None, t0, status)
    acc.save()
    nd = len(acc.done)
    json.dump(dict(done=nd, total=ntot, status='done' if nd >= ntot else 'stopped', started=t0, note=f'{acc.nobj:,} quasars'),
              open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    agg_status()
    print(f"{'stopped' if STOP['flag'] else 'done'}: {nd}/{ntot} files, {acc.nobj:,} quasars, norm src {acc.nsrc.tolist()}, {time.time()-t0:.0f}s", flush=True)


def agg_status():
    """Sum of the shard status files -> step_qso_stack.json (the monitor's quasar 'stacks' gauge)."""
    js = []
    pat = f'{STATUS_DIR}/{SNAME}_*_shard*.json'
    paths = [p for p in glob.glob(pat) if os.path.basename(p).split('_shard')[0].rsplit('_', 1)[0] == SNAME]
    for p in paths:
        try:
            js.append(json.load(open(p)))
        except Exception:
            pass
    if not js:
        return
    st = 'done' if all(j['status'] == 'done' for j in js) else ('running' if any(j['status'] == 'running' for j in js) else 'waiting')
    out = dict(done=sum(j['done'] for j in js), total=sum(j['total'] for j in js), status=st, started=min(j['started'] for j in js),
               note=', '.join(sorted(os.path.basename(p)[len(SNAME) + 1:-5] + ' ' + j['note'] for p, j in zip(paths, js))))
    p = f'{STATUS_DIR}/{SNAME}.json'; tmp = f'{p}.{os.getpid()}.tmp'      # one tmp per process (shards race)
    json.dump(out, open(tmp, 'w')); os.replace(tmp, p)


def merge(accdir, out='qso_stacks'):
    cells = pd.read_csv(f'{accdir}/cells_index.csv')
    tot = None; nobj = 0; nsrc = np.zeros(3, np.int64); parts = []
    for p in sorted(glob.glob(f'{accdir}/acc_*_shard*.npz')):
        z = np.load(p); parts.append(os.path.basename(p))
        if tot is None:
            tot = {k: z[k].astype(np.float64) for k in ('s1', 's2', 'n', 'wf', 'w')}
        else:
            for k in tot:
                tot[k] += z[k]
        nobj += int(z['nobj']); nsrc += z['nsrc']
    n = tot['n']
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, tot['s1'] / n, np.nan); sd = np.sqrt(np.maximum(np.where(n > 1, tot['s2'] / n - mean**2, np.nan), 0))
        wmean = np.where(tot['w'] > 0, tot['wf'] / tot['w'], np.nan); werr = np.where(tot['w'] > 0, 1 / np.sqrt(tot['w']), np.nan)
    np.savez(f'{OUT}/{out}.npz', WAVE=Q.QWAVE, MEAN=mean.astype(np.float32), SD=sd.astype(np.float32), N=n.astype(np.int32),
             WMEAN=wmean.astype(np.float32), WERR=werr.astype(np.float32), STACK_CELL=cells.STACK_CELL.to_numpy())
    info = pd.read_csv(f'{PROD}/qso_stack_cells_all.csv') if os.path.exists(f'{PROD}/qso_stack_cells_all.csv') else None
    cells['N_MAX'] = n.max(axis=1)
    cells = cells.merge(info, on='STACK_CELL', how='left') if info is not None else cells
    cells.to_csv(f'{OUT}/{out}_cells.csv', index=False)
    print(f'merged {parts}: {nobj:,} quasars, norm src {nsrc.tolist()}, {len(cells)} cells -> {OUT}/{out}.npz')


if __name__ == '__main__':
    main()
