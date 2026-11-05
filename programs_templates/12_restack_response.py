#!/usr/bin/env python3
"""Re-stack the DESI galaxies with the per-coadd spectral response of Calibration 2 (stellar thread, step8_response.py),
from the rest-frame store written by pass C (no coadd reads on exdisk2).

Pass C calibrated each galaxy with a fiber-flux SCALE only (the per-object tilt was noisy; CAL_B_EST). Calibration 2 solves
scale + tilt per (EXPID, PETAL_LOC) from the stars and averages it per coadd:
  response_coadd_desi.parquet  key (TARGETID, PROGRAM); CORR(l_obs) = 10^(-0.4 (A + B_MAG_DEX log10(l_obs/6000)))
Only the shape B_MAG_DEX is applied (the scale is already set by the fiber flux, and every spectrum is renormalized
in its window anyway). For each stored galaxy (normalized rest-frame flux on the pass-C grid):
  f_rest(l) *= 10^(-0.4 B log10(l (1+z) / 6000)), ivar /= that factor^2, then renormalize in the z-slice window
and accumulate in the SAME cell (ACC_INDEX of the store: the tilt moves Dn4000, the EWs and M by much less than a bin).
Galaxies whose coadd has no solved tilt (B NaN) enter uncorrected and are counted (N_NOB).

Check: --response none applies no correction and must reproduce the pass-C accumulators (up to float16 store precision).
Output: <out>/acc/acc_shard<i>of<N>.npz (same format as pass C, so 09_merge_stacks.py --passc <out> works).
Per-galaxy normalization record (docs/HANDOVER_normalization.md step 2), <out>/normtab/norm_shard<i>of<N>.parquet: one row per
stacked galaxy with NORM_SCALE (pass C) x RESCALE (this renormalization after the tilt) = the factor the rest-frame
F_lambda (pb.FLUX_UNIT) was multiplied by, and the mean of the final normalized spectrum (W*) and the good-pixel
fraction (COV*) in every galaxy window (WINDOWS below); the own window reads 1.
The window is the one normwin gives now: v2 = two points (5400-5500 below z = 0.30, 3000-3100 above; Nao 2026-10-31), so
this also moves the store's v1 five-window normalization to the two windows (out: ~/data/sdssdesi_passC_v2_2win).
Usage: python3 12_restack_response.py --workers 4 [--response FILE|none] [--out ~/data/sdssdesi_passC_v2_resp] [--limit N]
"""
import os, sys, glob, json, time, argparse
import numpy as np
import pandas as pd
from multiprocessing import Process

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pc = importlib.import_module('08_passC_stack')
pb = pc.pb
import gridcells as gc, normwin

STORE = '/Volumes/exdisk3/sdssdesi_restframe_v2'        # v1 store (whole-pixel shift): /Volumes/exdisk3/sdssdesi_restframe
RESPONSE = os.path.expanduser('~/data/sdssdesi_products/stars/response_coadd_desi.parquet')
STATUS = os.path.expanduser('~/data/sdssdesi_status')
WINDOWS = {'6125': (6000, 6250), '5450': (5400, 5500), '4700': (4600, 4800), '3050': (3000, 3100), '2150': (2100, 2200)}
TABCOLS = ['TARGETID', 'Z', 'ZB', 'NORM_SCALE', 'MABS', 'MBAND', 'MABS_R', 'MABS_G', 'MABS_U', 'MABS_NUV', 'COLOR', 'EW', 'DN4000',
           'MB', 'CB', 'EB', 'DB', 'SNR_SDSSPIX']


def load_response(path):
    """{program: (sorted TARGETID, B_MAG_DEX)}; looked up with searchsorted (53M rows: a dict would cost GBs per worker)."""
    if path == 'none':
        return None
    r = pd.read_parquet(path, columns=['TARGETID', 'PROGRAM', 'B_MAG_DEX'])
    out = {}
    for p, g in r.groupby('PROGRAM'):
        t = g.TARGETID.to_numpy(np.int64); o = np.argsort(t)
        out[p] = (t[o], g.B_MAG_DEX.to_numpy(np.float32)[o])
    return out


def lookup_b(rp, tid):
    t, b = rp
    if len(t) == 0:
        return np.full(len(tid), np.nan, np.float32)
    k = np.clip(np.searchsorted(t, tid), 0, len(t) - 1)
    return np.where(t[k] == tid, b[k], np.nan)


def worker(i, n, files, resp_path, out, nc, passc):
    resp = load_response(resp_path)
    acc = dict(wf=np.zeros((nc, pb.NREST), pc.F32), w=np.zeros((nc, pb.NREST), pc.F32), wf2=np.zeros((nc, pb.NREST), pc.F32),
               n=np.zeros((nc, pb.NREST), np.int32), ncell=np.zeros(nc, np.int64))
    done = []; ngal = nacc = nnob = 0; t0 = time.time(); allrec = []
    wsel = {k: (int(np.searchsorted(pb.RWAVE, lo)), int(np.searchsorted(pb.RWAVE, hi, side='right'))) for k, (lo, hi) in WINDOWS.items()}
    lw = np.log10(pb.RWAVE)
    wins = {zb: normwin.slice_window(zb * gc.DZ, (zb + 1) * gc.DZ, 'GALAXY') for zb in range(int(gc.ZMAX / gc.DZ))}
    stat = os.path.join(out, f'status_shard{i}of{n}.json') if 'test' in out else os.path.join(STATUS, f'step_gal_restack_shard{i}of{n}.json')
    for k, (prog, fn) in enumerate(files):
        hpx = int(os.path.basename(fn)[3:-4])
        s = np.load(fn)
        tab = pd.read_parquet(os.path.join(passc, 'tables', prog, f'hpx{hpx}.parquet'), columns=TABCOLS)
        zmap = dict(zip(tab.TARGETID.to_numpy(), zip(tab.Z.to_numpy(float), tab.ZB.to_numpy(int))))
        trow = dict(zip(tab.TARGETID.to_numpy(), range(len(tab))))
        rp = resp.get(prog, (np.zeros(0, np.int64), np.zeros(0, np.float32))) if resp is not None else None
        tid, ci, st, off = s['TARGETID'], s['ACC_INDEX'], s['START'], s['OFFSET']
        bb = lookup_b(rp, tid.astype(np.int64)) if rp is not None else None
        F, I = s['FLUX'].astype(np.float32), s['IVAR'].astype(np.float32); recs = []
        for j in range(len(tid)):
            ngal += 1
            c = int(ci[j])
            if c < 0:
                continue
            z, zb = zmap[tid[j]]
            a0 = int(st[j]); f = F[off[j]:off[j + 1]].copy(); iv = I[off[j]:off[j + 1]].copy(); a1 = a0 + len(f)
            if rp is not None:
                b = bb[j]
                if np.isfinite(b):
                    corr = 10.0**(-0.4 * b * (lw[a0:a1] + np.log10(1.0 + z) - np.log10(6000.0)))
                    f *= corr; iv /= corr**2
                else:
                    nnob += 1
            ff = np.zeros(pb.NREST, np.float32); ii = np.zeros(pb.NREST, np.float32); ff[a0:a1] = f; ii[a0:a1] = iv
            sc, ok = normwin.normalize(pb.RWAVE, ff, ii, wins[zb], min_good=normwin.MIN_GOOD_GALAXY)
            if not ok:
                continue
            ff *= sc; ii /= sc**2
            rec = dict(TARGETID=int(tid[j]), PROGRAM=prog, HPX=hpx, ACC_INDEX=c, RESCALE=float(sc),
                       B_MAG_DEX=float(bb[j]) if bb is not None else np.nan)
            for wk, (s0, s1) in wsel.items():
                g = ii[s0:s1] > 0; ng = int(g.sum())
                rec['COV' + wk] = ng / (s1 - s0); rec['W' + wk] = float(ff[s0:s1][g].mean()) if ng else np.nan
            rec['_row'] = trow[tid[j]]; recs.append(rec)
            acc['wf'][c] += ii * ff; acc['w'][c] += ii; acc['wf2'][c] += ii * ff**2; acc['n'][c] += ii > 0; acc['ncell'][c] += 1
            nacc += 1
        if recs:
            R = pd.DataFrame(recs); T = tab.iloc[R.pop('_row').to_numpy()].reset_index(drop=True).drop(columns='TARGETID')
            allrec.append(pd.concat([R, T], axis=1))
        done.append(f'{prog}/{hpx}')
        if k % 50 == 0 or k == len(files) - 1:
            with open(stat + '.tmp', 'w') as fo:
                json.dump(dict(done=k + 1, total=len(files), status='running' if k < len(files) - 1 else 'done',
                               note=f'galaxies {ngal:,}, stacked {nacc:,}, no tilt {nnob:,}', started=t0), fo)
            os.replace(stat + '.tmp', stat)
            print(f'shard {i}: {k+1}/{len(files)} galaxies {ngal:,} stacked {nacc:,} no-tilt {nnob:,} {time.time()-t0:.0f}s', flush=True)
    meta = dict(cells=pc.CELLTABLE, acc_min_n=pc.ACC_MIN_N, nrest=pb.NREST, rid_min=pb.RID_MIN, dlog=pc.DLOG, lines='vacuum',
                ncells=nc, writer='12_restack_response', response=resp_path, n_no_tilt=nnob)
    pc.save_acc(os.path.join(out, 'acc', f'acc_shard{i}of{n}.npz'), acc, set(done), meta)
    os.makedirs(os.path.join(out, 'normtab'), exist_ok=True)
    if allrec:
        nt = pd.concat(allrec, ignore_index=True); nt['SURVEY'] = 'DESI'
        nt.to_parquet(os.path.join(out, 'normtab', f'norm_shard{i}of{n}.parquet'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--response', default=RESPONSE)
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC_v2_resp')); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--store', default=STORE); ap.add_argument('--passc', default=pc.OUTDIR, help='pass-C dir with tables/')
    a = ap.parse_args()
    files = [(p, f) for p in ('bright', 'dark') for f in sorted(glob.glob(os.path.join(a.store, p, 'hpx*.npz')))]
    if a.limit:
        files = files[:a.limit]
    cells, _ = pc.load_cells(); nc = len(cells)
    os.makedirs(os.path.join(a.out, 'acc'), exist_ok=True)
    print(f'{len(files)} store files, {a.workers} workers, response {a.response}', flush=True)
    ps = [Process(target=worker, args=(i, a.workers, files[i::a.workers], a.response, a.out, nc, a.passc)) for i in range(a.workers)]
    [p.start() for p in ps]; [p.join() for p in ps]
    sys.exit(max(p.exitcode for p in ps))


if __name__ == '__main__':
    main()
