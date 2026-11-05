#!/usr/bin/env python3
"""SDSS galaxies: renormalize in the v2 two-point windows and re-stack, from the SDSS v2 rest-frame store (11_passG_sdss.py
--store), with the per-galaxy normalization record (docs/HANDOVER_normalization.md step 2). The DESI counterpart is
12_restack_response.py.

The store holds each stacked galaxy's rest-frame F_lambda x NORM_SCALE (unit mean in its v1 window). Here every galaxy is
renormalized to unit mean in the window normwin gives now (5400-5500 below z = 0.30, 3000-3100 above; Nao 2026-10-31),
RESCALE = that factor, and accumulated into its cell (ACC_INDEX) as in pass G. One row per stacked galaxy: NORM_SCALE x
RESCALE (the factor the rest-frame F_lambda was multiplied by), the mean of the final normalized spectrum (W*) and its
good-pixel fraction (COV*) in every galaxy window (12_restack_response.WINDOWS), plus the table columns.
Outputs: <out>/acc/acc_<survey>_shard<i>of<n>.npz (09_merge_stacks.py --passc <out> --acc 'acc_dr*.npz'),
         <out>/normtab/norm_<survey>.parquet
Usage: python3 16_norm_table_sdss.py [--nproc 4] [--out ~/data/sdssdesi_passG_sdss_v2_2win] [--limit N]
"""
import os, sys, glob, argparse, time
import numpy as np
import pandas as pd
from multiprocessing import Process

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pc = importlib.import_module('08_passC_stack'); pb = pc.pb
rs = importlib.import_module('12_restack_response')
import normwin, gridcells as gc

STORE = '/Volumes/exdisk3/sdssdesi_restframe_sdss_v2'
PASSG = os.path.expanduser('~/data/sdssdesi_passG_sdss_v2')
COLS = ['PLATE', 'MJD', 'FIBERID', 'Z', 'ZB', 'NORM_SCALE', 'MABS', 'MBAND', 'MABS_R', 'MABS_G', 'MABS_U', 'MABS_NUV', 'COLOR',
        'EW', 'DN4000', 'MB', 'CB', 'EB', 'DB', 'SNR_SDSSPIX']
WSEL = {k: (int(np.searchsorted(pb.RWAVE, lo)), int(np.searchsorted(pb.RWAVE, hi, side='right'))) for k, (lo, hi) in rs.WINDOWS.items()}


def one(survey, fn, acc, wins):
    base = os.path.basename(fn)[:-4]
    tp = os.path.join(PASSG, 'tables', survey, base + '.parquet')
    if not os.path.exists(tp):
        return None
    s = np.load(fn); tab = pd.read_parquet(tp, columns=COLS)
    row = dict(zip(tab.FIBERID.to_numpy(), range(len(tab))))
    fib, ci, st, off = s['FIBERID'], s['ACC_INDEX'], s['START'], s['OFFSET']
    F, I = s['FLUX'].astype(np.float32), s['IVAR'].astype(np.float32)
    recs = []
    for j in range(len(fib)):
        c = int(ci[j])
        if c < 0:
            continue
        r = row[fib[j]]
        a0 = int(st[j]); f = F[off[j]:off[j + 1]]; iv = I[off[j]:off[j + 1]]; a1 = a0 + len(f)
        ff = np.zeros(pb.NREST, np.float32); ii = np.zeros(pb.NREST, np.float32); ff[a0:a1] = f; ii[a0:a1] = iv
        sc, ok = normwin.normalize(pb.RWAVE, ff, ii, wins[int(tab.ZB.iat[r])], min_good=normwin.MIN_GOOD_GALAXY)
        if not ok:
            continue
        ff *= sc; ii /= sc**2
        rec = dict(ACC_INDEX=c, RESCALE=float(sc), _row=r)
        for wk, (s0, s1) in WSEL.items():
            g = ii[s0:s1] > 0; ng = int(g.sum())
            rec['COV' + wk] = ng / (s1 - s0); rec['W' + wk] = float(ff[s0:s1][g].mean()) if ng else np.nan
        recs.append(rec)
        acc['wf'][c] += ii * ff; acc['w'][c] += ii; acc['wf2'][c] += ii * ff**2; acc['n'][c] += ii > 0; acc['ncell'][c] += 1
    if not recs:
        return None
    R = pd.DataFrame(recs); T = tab.iloc[R.pop('_row').to_numpy()].reset_index(drop=True)
    out = pd.concat([R, T], axis=1); out['SURVEY'] = survey.upper()
    return out


def worker(i, n, survey, files, out, nc):
    acc = dict(wf=np.zeros((nc, pb.NREST), pc.F32), w=np.zeros((nc, pb.NREST), pc.F32), wf2=np.zeros((nc, pb.NREST), pc.F32),
               n=np.zeros((nc, pb.NREST), np.int32), ncell=np.zeros(nc, np.int64))
    wins = {zb: normwin.slice_window(zb * gc.DZ, (zb + 1) * gc.DZ, 'GALAXY') for zb in range(int(gc.ZMAX / gc.DZ))}
    parts, done, t0 = [], [], time.time()
    for k, fn in enumerate(files):
        p = one(survey, fn, acc, wins)
        if p is not None:
            parts.append(p)
        done.append(os.path.basename(fn).replace('.npz', '.fits'))
        if k % 200 == 0 or k == len(files) - 1:
            print(f'{survey} shard {i}: {k+1}/{len(files)} plates, {sum(len(x) for x in parts):,} stacked, {time.time()-t0:.0f}s', flush=True)
    meta = dict(cells=pc.CELLTABLE, acc_min_n=pc.ACC_MIN_N, nrest=pb.NREST, rid_min=pb.RID_MIN, dlog=pc.DLOG, lines='vacuum',
                ncells=nc, survey=survey, writer='16_norm_table_sdss', windows=str(normwin.GALAXY_V2), switch_z=normwin.GALAXY_SWITCH_Z)
    pc.save_acc(os.path.join(out, 'acc', f'acc_{survey}_shard{i}of{n}.npz'), acc, set(done), meta)
    if parts:
        pd.concat(parts, ignore_index=True).to_parquet(os.path.join(out, 'normtab', f'norm_{survey}_part{i}.parquet'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passG_sdss_v2_2win'))
    a = ap.parse_args()
    for d in ('acc', 'normtab'):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    cells, _ = pc.load_cells(); nc = len(cells)
    for survey in ('dr8', 'dr17'):
        files = sorted(glob.glob(os.path.join(STORE, survey, '*.npz')))
        if a.limit:
            files = files[:a.limit]
        ps = [Process(target=worker, args=(i, a.nproc, survey, files[i::a.nproc], a.out, nc)) for i in range(a.nproc)]
        [p.start() for p in ps]; [p.join() for p in ps]
        assert all(p.exitcode == 0 for p in ps), f'{survey}: a worker failed'
        pp = sorted(glob.glob(os.path.join(a.out, 'normtab', f'norm_{survey}_part*.parquet')))
        t = pd.concat([pd.read_parquet(f) for f in pp], ignore_index=True)
        t.to_parquet(os.path.join(a.out, 'normtab', f'norm_{survey}.parquet')); [os.remove(f) for f in pp]
        print(f'{survey}: {len(files)} plates, {len(t):,} stacked galaxies', flush=True)


if __name__ == '__main__':
    main()
