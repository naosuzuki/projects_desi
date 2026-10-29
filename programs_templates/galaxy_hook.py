"""Galaxy hook for the shared DESI walk (programs_quasars/02_passQ_desi.py --with-galaxies): pass C's per-coadd galaxy
work (08_passC_stack.process_arrays) on the spectra the walk has already read, so each coadd is read once for quasars,
stars and galaxies.

State is the SAME as the standalone pass C with N shards: accumulators ~/data/sdssdesi_passC/acc/acc_shard<i>of<N>.npz
(done set + per-cell sums) and per-galaxy tables ~/data/sdssdesi_passC/tables/<program>/hpx<N>.parquet. The walk's
shard i/N and run_passC.sh N use the same manifest split, so whatever the walk does, a later `run_passC.sh N` skips.

Use in the walk:
    G = galaxy_hook.GalaxyHook(shard_i, shard_n, store='/Volumes/exdisk3/sdssdesi_restframe_v2', zl=zl, offsets=offsets)
    if G.needs(prog, hpx): ... read FIBERMAP (galaxy_hook.FM_COLS) before closing the FITS, then after combine_cameras:
        G.process(path, fm, w, flux, ivar, mask, prog, hpx)
    G.close()      # at the end and on SIGTERM: writes the accumulator checkpoint
"""
import os, sys, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pc = importlib.import_module('08_passC_stack')
FM_COLS = pc.FM_COLS
CHECKPOINT_EVERY = pc.CHECKPOINT_EVERY


class GalaxyHook:
    def __init__(self, shard_i, shard_n, store=None, zl=None, offsets=None, outdir=pc.OUTDIR):
        self.si, self.sn, self.store, self.outdir = shard_i, shard_n, store, outdir
        self.zl = zl if zl is not None else dict(np.load(pc.pb.ZLOOKUP_COMPACT))
        self.offsets = offsets if offsets is not None else {int(k): (int(s), int(e)) for k, s, e in
                                                            zip(self.zl['OFF_KEY'], self.zl['OFF_START'], self.zl['OFF_STOP'])}
        cells, self.cellmap = pc.load_cells(); nc = len(cells)
        self.accp = os.path.join(outdir, 'acc', f'acc_shard{shard_i}of{shard_n}.npz')
        os.makedirs(os.path.dirname(self.accp), exist_ok=True)
        if os.path.exists(self.accp):
            z = np.load(self.accp, allow_pickle=False)
            self.acc = {k: z[k] for k in ('wf', 'w', 'wf2', 'n', 'ncell')}; self.done = set(z['done'].tolist())
        else:
            self.acc = dict(wf=np.zeros((nc, pc.NREST), pc.F32), w=np.zeros((nc, pc.NREST), pc.F32), wf2=np.zeros((nc, pc.NREST), pc.F32),
                            n=np.zeros((nc, pc.NREST), np.int32), ncell=np.zeros(nc, np.int64)); self.done = set()
        self.meta = dict(cells=pc.CELLTABLE, acc_min_n=pc.ACC_MIN_N, nrest=pc.NREST, rid_min=pc.pb.RID_MIN, dlog=pc.DLOG,
                         lines='vacuum', ncells=nc, writer='galaxy_hook (shared DESI walk)')
        self.since = 0; self.ngal = 0; self.nacc = 0

    def needs(self, prog, hpx):
        return f'{prog}/{int(hpx)}' not in self.done

    def galaxy_rows(self, fm, path, prog, hpx):
        return pc.galaxy_rows(fm, path, prog, int(hpx), self.zl, self.offsets)

    def process(self, path, fm, w, flux, ivar, mask, prog, hpx):
        """fm: FIBERMAP record array with FM_COLS (all rows); w/flux/ivar/mask: combine_cameras output for all rows."""
        hpx = int(hpx); key = f'{prog}/{hpx}'
        if key in self.done:
            return 0
        df = self.galaxy_rows(fm, path, prog, hpx)
        if len(df):
            tab, ng, na = pc.process_arrays(df, w, flux, ivar, mask, prog, hpx, self.cellmap, self.acc, self.store)
            td = os.path.join(self.outdir, 'tables', prog); os.makedirs(td, exist_ok=True)
            out = os.path.join(td, f'hpx{hpx}.parquet'); tab.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
            self.ngal += ng; self.nacc += na
        else:                                         # empty table too: the monitor counts one table per coadd
            td = os.path.join(self.outdir, 'tables', prog); os.makedirs(td, exist_ok=True)
            out = os.path.join(td, f'hpx{hpx}.parquet'); df.to_parquet(out + '.tmp'); os.replace(out + '.tmp', out)
        self.done.add(key); self.since += 1
        if self.since >= CHECKPOINT_EVERY:
            self.checkpoint()
        return len(df)

    def checkpoint(self):
        pc.save_acc(self.accp, self.acc, self.done, self.meta); self.since = 0

    def close(self):
        self.checkpoint()
