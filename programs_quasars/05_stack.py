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
v3 (--v3, Nao 2026-10-31): v2 rest frame, cells from qso_cells_v3_all (04_cells.py --chain: no redshift axis). The window is
chosen by each quasar's own Z, and the normalized spectrum is divided by 10**CHAIN_LOGR (the fitted R = F3025/F(window)),
so every member reads ~1 at 3000-3050 A and a cell joins all redshifts. No store by default. Outputs PROD/v3.
v3 measurement pass (--v3 --measure, Nao 2026-10-31 'Yes, re-stack'): the same read, but per training quasar only the
ivar-weighted means of its own window (NORM_OWN, SRC_OWN) and of 3000-3050 A (M3025, COV3025 = fraction of good pixels)
are recorded (per-shard *_meas.npz, checkpointed); --merge --measure -> PROD/v3/qso_norm_table.parquet (exact IDs).
04_cells.py --chain then chains with the MEASURED ratio M3025/NORM_OWN wherever 3000-3050 A is covered (CHAIN_SRC 0) and
with the fitted R(M) only where it is not (z < ~0.19, z > ~2.25; CHAIN_SRC 1), so the stack normalization is a direct
3000-3050 A normalization for every quasar that reaches the window.
v3 template normalization (--v3 --template STACKS, Nao 2026-11-01 'Yes, template norm'): the 50 A window mean is noisy for the
faint low-z quasars (window at the DESI blue edge, S/N of a few), and dividing by a noisy mean biases the stack outside the
window (faint stacks jumped 20-130 % at the window edge, which the EM-PCA e1/e3/e4 inherited). Instead each member is
normalized by the ivar-weighted least-squares amplitude of its own stack's mean spectrum (WMEAN of STACKS, scaled to unit
mean at 3000-3050 A where it covers the window) over all its good pixels: the whole spectrum sets the norm, every redshift
alike, no R(M) fallback. Falls back to the window x chain norm where the template has < TPL_MIN_PIX pixels in common
(NORM_SRC 1). The per-member norms go to the accumulator files (ID, NORM, NORM_SRC) and --merge writes them to
v3/qso_norm_table_template.parquet.
Usage:
  python3 05_stack.py --survey sdss --shard 0/4          # spPlates, internal disk
  python3 05_stack.py --survey desi --shard 0/2          # DESI coadds, exdisk2 (read-only)
  python3 05_stack.py --merge                            # all shard accumulators -> qso_stacks.npz + qso_stacks_cells.csv
  python3 05_stack.py --v3 --measure --survey desi --shard 0/4 --accdir DIR; ... --v3 --measure --merge --accdir DIR
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
V3 = False
MEASURE = False
MEAS_WIN = (3000.0, 3050.0)
TPL = None            # --template: [ncell, QNREST] template per cell (0 where undefined)
TELL = None           # --telluric: programs_telluric.tellapply module (Stellar thread); C = 1 + T0 + sum a_k T_k per spectrum, observed frame
TELL_LOG0 = 3.5502; TELL_NPIX = 4447                                # the telluric lattice: log10 lambda_obs = 3.5502 + 1e-4 i
BROAD = [1033.8, 1215.67, 1240.8, 1397.0, 1549.06, 1640.4, 1908.73, 2798.75, 3426.9, 3727.1, 3869.9, 4102.9, 4341.7,
         4862.7, 4960.3, 5008.2, 6564.6]                             # rest A, ivar = 0 within +-5000 km/s in the telluric fit
TPL_MIN_PIX = 300
CELLS = 'qso_cells_all'
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
    c = pd.read_parquet(f'{PROD}/{CELLS}.parquet', columns=['PLATE', 'MJD', 'FIBERID', 'TARGETID', 'PROGRAM', 'SURVEY', 'TRAIN',
                                                            'STACK_CELL', 'STACK_ZBIN', 'STACK_N'] + (['Z', 'CHAIN_LOGR'] if V3 else []))
    c = c[c.TRAIN & ((c.STACK_CELL != '') | MEASURE)].copy()      # --measure: every training quasar, pooled ones too
    cells = np.array(sorted(c.STACK_CELL.unique()))
    c['CELL_I'] = np.searchsorted(cells, c.STACK_CELL.to_numpy()).astype(np.int32)
    if V3:           # own redshift picks the window (no z axis); Z renamed so the master-table merge keeps its own Z
        zc = c.pop('Z').to_numpy(dtype=float)
    else:
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


def window_mean(rflux, rivar, lo, hi):
    """ivar-weighted mean of each row in [lo, hi] and the fraction of good pixels there (no fallback)."""
    sel = np.where((Q.QWAVE >= lo) & (Q.QWAVE <= hi))[0]
    f = rflux[:, sel]; iv = rivar[:, sel]; good = iv > 0
    with np.errstate(all='ignore'):
        m = (f * iv).sum(axis=1) / np.where(good.any(axis=1), (iv * good).sum(axis=1), np.nan)
    return m, good.mean(axis=1)


class Meas:
    """--measure: per-quasar window means instead of stack accumulators (same done / save / nobj / nsrc interface)."""
    COLS = ('ID', 'Z', 'NORM_OWN', 'SRC_OWN', 'COV_OWN', 'M3025', 'COV3025')

    def __init__(self, path, ncell=None):
        self.path = path; self.rows = {k: [] for k in self.COLS}; self.prog = []
        if os.path.exists(path):
            z = np.load(path, allow_pickle=False)
            for k in self.COLS:
                self.rows[k].append(z[k])
            self.prog.append(z['PROG']); self.done = set(z['done'].tolist()); self.nobj = int(z['nobj']); self.nsrc = z['nsrc']
        else:
            self.done = set(); self.nobj = 0; self.nsrc = np.zeros(3, np.int64)

    def add(self, ids, prog, z, norm, src, cov, m3025, cov3025):
        for k, v in zip(self.COLS, (ids, z, norm, src, cov, m3025, cov3025)):
            self.rows[k].append(np.asarray(v))
        self.prog.append(np.asarray(prog).astype('U8'))

    def save(self):
        tmp = self.path + '.tmp.npz'
        arr = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in self.rows.items()}
        np.savez(tmp, **arr, PROG=np.concatenate(self.prog) if self.prog else np.zeros(0, 'U8'),
                 done=np.array(sorted(self.done)), nobj=self.nobj, nsrc=self.nsrc)
        os.replace(tmp, self.path)


class Acc:
    def __init__(self, path, ncell):
        self.path = path
        if os.path.exists(path):
            z = np.load(path, allow_pickle=False)
            self.a = {k: z[k] for k in ('s1', 's2', 'n', 'wf', 'w')}; self.done = set(z['done'].tolist())
            self.nobj = int(z['nobj']); self.nsrc = z['nsrc']
            self.ids = [z['ID']] if 'ID' in z.files else []; self.norms = [z['NORM']] if 'ID' in z.files else []; self.srcs = [z['NORM_SRC']] if 'ID' in z.files else []
        else:
            shp = (ncell, Q.QNREST)
            self.a = dict(s1=np.zeros(shp, np.float32), s2=np.zeros(shp, np.float32), n=np.zeros(shp, np.int32),
                          wf=np.zeros(shp, np.float32), w=np.zeros(shp, np.float32))
            self.done = set(); self.nobj = 0; self.nsrc = np.zeros(3, np.int64); self.ids = []; self.norms = []; self.srcs = []
        self.ntell = int(z['ntell']) if os.path.exists(path) and 'ntell' in z.files else 0

    def add(self, cell_i, fn, ivn):
        good = ivn > 0
        for k in np.unique(cell_i):
            r = cell_i == k; f = np.where(good[r], fn[r], 0.0); g = good[r]; iv = np.where(g, ivn[r], 0.0)
            self.a['s1'][k] += f.sum(0); self.a['s2'][k] += (f * f).sum(0); self.a['n'][k] += g.sum(0)
            self.a['wf'][k] += (iv * f).sum(0); self.a['w'][k] += iv.sum(0)

    def save(self):
        tmp = self.path + '.tmp.npz'
        extra = {}
        if self.ids:
            extra = dict(ID=np.concatenate(self.ids), NORM=np.concatenate(self.norms), NORM_SRC=np.concatenate(self.srcs))
        np.savez(tmp, **self.a, done=np.array(sorted(self.done)), nobj=self.nobj, nsrc=self.nsrc, ntell=self.ntell, **extra)
        os.replace(tmp, self.path)


def template_norm(rflux, rivar, cell_i):
    """ivar-weighted least-squares amplitude of each row's cell template over its good pixels: norm = sum(iv f T) / sum(iv T^2).
    NaN where fewer than TPL_MIN_PIX pixels are shared or the amplitude is not positive."""
    T = TPL[cell_i].astype(np.float64); good = (rivar > 0) & (T != 0)
    iv = np.where(good, rivar, 0.0).astype(np.float64)              # float64: cgs ivar (~1e34-1e36) x T^2 overflows float32
    with np.errstate(all='ignore'):
        num = (iv * rflux * T).sum(axis=1); den = (iv * T * T).sum(axis=1)
        norm = num / den
    ok = (good.sum(axis=1) >= TPL_MIN_PIX) & np.isfinite(norm) & (norm > 0)
    return np.where(ok, norm, np.nan)


def _shift_rows(X, sh, n_out, fill=0.0, nearest=False):
    """X [n, npix] sampled at integer positions j; returns Y [n, n_out] with Y[:, i] = X at j = i + sh[row] (fractional,
    linear interpolation; nearest=True picks the nearer pixel, for ivar). fill outside."""
    n = len(X); i = np.arange(n_out)[None, :]; jf = i + sh[:, None]
    j0 = np.floor(jf).astype(np.int64); fr = jf - j0
    if nearest:
        j0 = np.rint(jf).astype(np.int64); fr = np.zeros_like(fr)
    ok0 = (j0 >= 0) & (j0 < X.shape[1]); ok1 = (j0 + 1 >= 0) & (j0 + 1 < X.shape[1])
    r = np.arange(n)[:, None]
    x0 = np.where(ok0, X[r, np.clip(j0, 0, X.shape[1] - 1)], fill); x1 = np.where(ok1, X[r, np.clip(j0 + 1, 0, X.shape[1] - 1)], fill)
    return np.where(fr > 0, x0 * (1 - fr) + x1 * fr, x0)


def telluric_correct(rflux, rivar, z, model, survey):
    """Stellar thread's telluric correction (programs_telluric/tellapply.py) applied in the rest frame: the rest spectrum and its
    model are put on the observed lattice (log10 lambda = 3.5502 + 1e-4 i, the same 1e-4 step, a fractional shift
    log10(1+z)/1e-4), tellapply gives C = 1 + T0 + sum a_k T_k per spectrum there (1 outside its pixels), and C is shifted
    back: flux / C, ivar x C^2. ivar = 0 on the broad lines and the forest for the fit. Rows without a model (0) are untouched."""
    n = len(rflux); sh = (TELL_LOG0 / Q.DLOG - Q.QRID_MIN) - np.log10(1.0 + z) / Q.DLOG       # rest j = obs i + sh
    has = (model != 0).any(axis=1)
    if not has.any():
        return rflux, rivar, 0
    f = _shift_rows(rflux, sh, TELL_NPIX); m = _shift_rows(model, sh, TELL_NPIX)
    iv = _shift_rows(np.where((rivar > 0) & (model != 0), rivar, 0.0), sh, TELL_NPIX, nearest=True)
    lrest = TELL_LOG0 + Q.DLOG * np.arange(TELL_NPIX)[None, :] - np.log10(1.0 + z)[:, None]    # rest log10 lambda of each obs pixel
    bad = lrest < np.log10(1215.67); v = np.log10(1 + 5000.0 / 2.998e5)
    for l in BROAD:
        bad |= np.abs(lrest - np.log10(l)) < v
    iv[bad] = 0.0
    wave = 10.0**(TELL_LOG0 + Q.DLOG * np.arange(TELL_NPIX))
    C = np.ones((n, TELL_NPIX))
    C[has], _ = TELL.factor(survey, wave, f[has], iv[has], m[has])      # (C, a [n, 7]); a enters C only where tellapply.APPLY_A
    Cr = _shift_rows(C, -sh, rflux.shape[1], fill=1.0)
    Cr = np.where(np.isfinite(Cr) & (Cr > 0), Cr, 1.0)
    return rflux / Cr, rivar * Cr**2, int(has.sum())


def store_file(store, survey, base, ids, cell_i, norm, src, fn, ivn):
    """Normalized spectra, float16, cropped to each row's coverage (START = first QWAVE pixel, OFFSET into FLUX/IVAR)."""
    d = os.path.join(store, survey); os.makedirs(d, exist_ok=True)
    starts, flx, ivs = [], [], []
    for i in range(len(fn)):
        g = np.where(ivn[i] > 0)[0]
        a0, a1 = (g[0], g[-1] + 1) if len(g) else (0, 0)
        starts.append(a0); flx.append(fn[i, a0:a1].astype(np.float16)); ivs.append(np.minimum(ivn[i, a0:a1], 65000).astype(np.float16))   # float16 max 65504
    lens = np.array([len(x) for x in flx]); off = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
    out = os.path.join(d, base + '.npz'); tmp = out + '.tmp.npz'
    np.savez(tmp, ID=ids, CELL_I=cell_i, NORM=norm, NORM_SRC=src, START=np.array(starts, np.int32), OFFSET=off,
             FLUX=np.concatenate(flx) if flx else np.zeros(0, np.float16), IVAR=np.concatenate(ivs) if ivs else np.zeros(0, np.float16),
             QRID_MIN=Q.QRID_MIN)
    os.replace(tmp, out)


def finish_rows(acc, store, survey, base, ids, cell_i, rflux, rivar, sub):
    norm, src = normalize(rflux, rivar, sub.WIN_LO.to_numpy(), sub.WIN_HI.to_numpy(),
                          sub.CONT_LOGF3000.to_numpy(dtype=float), sub.ALPHA_LAM.to_numpy(dtype=float))
    if MEASURE:                   # v3 measurement pass: own window and 3000-3050 A per quasar, no stacking
        wlo, whi = sub.WIN_LO.to_numpy(), sub.WIN_HI.to_numpy(); cov_own = np.full(len(sub), np.nan)
        for lo, hi in set(zip(wlo, whi)):
            rows = np.where((wlo == lo) & (whi == hi))[0]; cov_own[rows] = window_mean(rflux[rows], rivar[rows], lo, hi)[1]
        m3025, cov3025 = window_mean(rflux, rivar, *MEAS_WIN)
        prog = sub.PROGRAM.to_numpy() if 'PROGRAM' in sub else sub.SURVEY.to_numpy()
        acc.add(np.asarray(ids, np.int64), prog, sub.Z.to_numpy(dtype=float), norm, src, cov_own, m3025, cov3025)
        acc.nobj += len(ids); acc.nsrc += np.bincount(src, minlength=3)
        return
    if 'CHAIN_LOGR' in sub:       # v3: chain to 3025 A
        norm = norm * 10**sub.CHAIN_LOGR.to_numpy(dtype=float)
    if TPL is not None:           # v3 template normalization: window x chain only as the fallback (src 1)
        tn = template_norm(rflux, rivar, cell_i)
        use = np.isfinite(tn)
        src = np.where(use, 0, np.where(src < 2, 1, 2)).astype(np.int8); norm = np.where(use, tn, norm)
        acc.ids.append(np.asarray(ids, np.int64)); acc.norms.append(norm.astype(np.float64)); acc.srcs.append(src.copy())
        if TELL is not None:      # run 5: telluric correction with the member's template at its amplitude as the model
            model = np.where(use[:, None], TPL[cell_i].astype(np.float64) * np.where(use, tn, 0.0)[:, None], 0.0)
            rflux, rivar, nt = telluric_correct(rflux, rivar, sub.Z.to_numpy(dtype=float), model, survey)
            acc.ntell += nt
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
                                           coeff0=float(h['COEFF0']) if V2 else None, survey=sub.SURVEY.iloc[0])
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
    ap.add_argument('--v3', action='store_true', help='v3: v2 rest frame, no z axis, chained to 3025 A (qso_cells_v3_all)')
    ap.add_argument('--merge', action='store_true')
    ap.add_argument('--measure', action='store_true', help='v3: record per-quasar window means (own window, 3000-3050 A) instead of stacking')
    ap.add_argument('--template', default=None, help='v3: stacks npz whose WMEAN (per cell, unit mean at 3000-3050 A) normalizes each member by least squares')
    ap.add_argument('--telluric', action='store_true', help='run 5: divide every member by the Stellar thread\'s telluric factor (programs_telluric/tellapply.py; needs --template)')
    ap.add_argument('--out', default='qso_stacks', help='--merge output name (qso_stacks_dark for the provisional merge)')
    a = ap.parse_args()
    global V2, V3, CELLS, OUT, SNAME, MEASURE, TPL
    MEASURE = a.measure
    if a.v3:
        V2, V3, CELLS, OUT, SNAME = True, True, 'qso_cells_v3_all', f'{PROD}/v3', 'step_qso_stack_v3' + ('_meas' if MEASURE else '')
        os.makedirs(OUT, exist_ok=True)
        if a.store is None:
            a.store = ''
    elif a.v2:
        V2, OUT, SNAME = True, f'{PROD}/v2', 'step_qso_stack_v2'
        os.makedirs(OUT, exist_ok=True)
    if a.accdir is None:
        a.accdir = f'{OUT}/stack_acc'
    if a.store is None:
        a.store = STORE + ('_v2' if V2 else '')
    if a.merge:
        return merge_measure(a.accdir) if MEASURE else merge(a.accdir, a.out)
    si, sn = map(int, a.shard.split('/'))
    t0 = time.time()
    mem, cells = members()
    if a.template:
        TPL = load_templates(a.template, cells)
    if a.telluric:
        global TELL
        if TPL is None:
            sys.exit('--telluric needs --template (the model)')
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'programs_telluric'))
        import tellapply; TELL = tellapply
        print('telluric correction: programs_telluric/tellapply.py', getattr(tellapply, '__version__', ''), flush=True)
    os.makedirs(a.accdir, exist_ok=True)
    pd.DataFrame({'CELL_I': np.arange(len(cells)), 'STACK_CELL': cells}).to_csv(f'{a.accdir}/cells_index.csv', index=False)
    if a.survey == 'desi':
        m = pd.read_parquet(f'{PROD}/qso_master_desi.parquet', columns=DESI_COLS)
        dm = mem[mem.SURVEY == 'desi'].copy()
        assert not pd.api.types.is_float_dtype(dm.TARGETID), 'float TARGETID in qso_cells_all: rebuild it with 04_cells.py'
        dm['TARGETID'] = dm.TARGETID.astype(np.int64)     # exact ids (04_cells Int64), never float
        sub = dm.merge(m, on=['TARGETID', 'PROGRAM'], how='inner')
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
    acc = (Meas(f'{a.accdir}/meas_{a.survey}_shard{si}of{sn}.npz') if MEASURE else
           Acc(f'{a.accdir}/acc_{a.survey}_shard{si}of{sn}.npz', len(cells)))
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


def load_templates(path, cells):
    """[ncell, QNREST] templates from a stacks npz: WMEAN where >= max(30, 10 % of the stack's best coverage) members, scaled to
    unit mean at 3000-3050 A where the stack covers the window (else left at its chained level); 0 (undefined) elsewhere."""
    z = np.load(path, allow_pickle=True); idx = {c: i for i, c in enumerate(z['STACK_CELL'])}
    wm, n = z['WMEAN'], z['N']; T = np.zeros((len(cells), Q.QNREST), np.float32)
    piv = (Q.QWAVE >= MEAS_WIN[0]) & (Q.QWAVE <= MEAS_WIN[1]); nscaled = 0
    for k, c in enumerate(cells):
        if c not in idx:
            continue
        i = idx[c]; ok = np.isfinite(wm[i]) & (n[i] >= max(30, 0.1 * n[i].max()))
        t = np.where(ok, wm[i], 0.0)
        if (ok & piv).sum() >= 0.5 * piv.sum() and t[ok & piv].mean() > 0:
            t /= t[ok & piv].mean(); nscaled += 1
        T[k] = t
    print(f'templates: {sum(c in idx for c in cells)} of {len(cells)} cells from {path}, {nscaled} scaled to unit mean at 3000-3050 A', flush=True)
    return T


def merge_measure(accdir):
    """Shard measurement files -> PROD/v3/qso_norm_table.parquet: one row per training quasar, exact IDs
    (TARGETID + PROGRAM for DESI, PLATE/MJD/FIBERID for SDSS), own-window norm and the 3000-3050 A mean."""
    parts = []
    for p in sorted(glob.glob(f'{accdir}/meas_*_shard*.npz')):
        z = np.load(p); sv = os.path.basename(p).split('_')[1]
        d = pd.DataFrame({k: z[k] for k in Meas.COLS}); d['PROG'] = z['PROG']; d['STORE_SURVEY'] = sv; parts.append(d)
        print(p, len(d), 'rows, done', len(z['done']))
    d = pd.concat(parts, ignore_index=True)
    desi = d.STORE_SURVEY == 'desi'
    d['TARGETID'] = pd.array(np.where(desi, d.ID, 0), dtype='Int64'); d.loc[~desi, 'TARGETID'] = pd.NA
    d['PROGRAM'] = np.where(desi, d.PROG, ''); d['SURVEY'] = np.where(desi, 'desi', d.PROG)
    sid = d.ID.to_numpy(np.int64)
    d['PLATE'] = pd.array(np.where(desi, 0, sid // 10**9), dtype='Int64'); d['MJD'] = pd.array(np.where(desi, 0, sid // 10**4 % 10**5), dtype='Int64')
    d['FIBERID'] = pd.array(np.where(desi, 0, sid % 10**4), dtype='Int64')
    for c in ('PLATE', 'MJD', 'FIBERID'):
        d.loc[desi, c] = pd.NA
    with np.errstate(all='ignore'):
        d['W3025'] = d.M3025 / d.NORM_OWN                      # measured R = F(3025) / F(own window); 1 for 0.8 <= z < 2
    d = d.drop(columns=['PROG']).rename(columns={'SRC_OWN': 'NORM_SRC'})
    os.makedirs(f'{PROD}/v3', exist_ok=True); out = f'{PROD}/v3/qso_norm_table.parquet'
    d.to_parquet(out, index=False)
    ok = (d.COV3025 >= MIN_WIN_GOOD) & np.isfinite(d.W3025) & (d.W3025 > 0)
    print(f'{len(d):,} quasars -> {out}; 3000-3050 A covered (>= {MIN_WIN_GOOD:.0%} good): {ok.sum():,} ({100*ok.mean():.1f} %); '
          f'own-window src {np.bincount(d.NORM_SRC, minlength=3).tolist()}')


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
    tn = [(np.load(p)['ID'], np.load(p)['NORM'], np.load(p)['NORM_SRC'], os.path.basename(p).split('_')[1]) for p in sorted(glob.glob(f'{accdir}/acc_*_shard*.npz')) if 'ID' in np.load(p).files]
    if tn:
        d = pd.concat([pd.DataFrame(dict(ID=i, NORM=v, NORM_SRC=s, STORE_SURVEY=sv)) for i, v, s, sv in tn], ignore_index=True)
        d.to_parquet(f'{OUT}/qso_norm_table_template.parquet', index=False)
        print(f'{len(d):,} member norms (template src {(d.NORM_SRC == 0).sum():,}, fallback {(d.NORM_SRC == 1).sum():,}) -> {OUT}/qso_norm_table_template.parquet')
    n = tot['n']
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, tot['s1'] / n, np.nan); sd = np.sqrt(np.maximum(np.where(n > 1, tot['s2'] / n - mean**2, np.nan), 0))
        wmean = np.where(tot['w'] > 0, tot['wf'] / tot['w'], np.nan); werr = np.where(tot['w'] > 0, 1 / np.sqrt(tot['w']), np.nan)
    np.savez(f'{OUT}/{out}.npz', WAVE=Q.QWAVE, MEAN=mean.astype(np.float32), SD=sd.astype(np.float32), N=n.astype(np.int32),
             WMEAN=wmean.astype(np.float32), WERR=werr.astype(np.float32), STACK_CELL=cells.STACK_CELL.to_numpy())
    ip = f"{PROD}/qso_stack_cells_{'v3_' if V3 else ''}all.csv"
    info = pd.read_csv(ip) if os.path.exists(ip) else None
    cells['N_MAX'] = n.max(axis=1)
    cells = cells.merge(info, on='STACK_CELL', how='left') if info is not None else cells
    cells.to_csv(f'{OUT}/{out}_cells.csv', index=False)
    print(f'merged {parts}: {nobj:,} quasars, norm src {nsrc.tolist()}, {len(cells)} cells -> {OUT}/{out}.npz')


if __name__ == '__main__':
    main()
