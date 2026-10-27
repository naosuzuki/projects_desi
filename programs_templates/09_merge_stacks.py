#!/usr/bin/env python3
"""Merge the pass-C shard accumulators into the cell stacks and write them per redshift slice.

Input : <passC>/acc/acc_shard*of*.npz  (from 08_passC_stack.py)
Output: <passC>/stacks/stacks_z<zb>.fits  one per Delta z = 0.05 slice, rows = cells of that slice
          HDU FLUX (weighted mean, normalized F_lambda), IVAR (sum of weights), SCATTER (weighted rms about the
          mean, per pixel), NSPEC (spectra contributing per pixel), WAVE (rest-frame vacuum A), CELLS (cell table rows)
        <passC>/stacks/cell_stacks.csv   cell table with N_STACKED and the median stack S/N
Cells with fewer than --floor stacked spectra are left out.
"""
import os, glob, json, argparse
import numpy as np
import pandas as pd
import fitsio

ap = argparse.ArgumentParser()
ap.add_argument('--passc', default=os.path.expanduser('~/data/sdssdesi_passC'))
ap.add_argument('--floor', type=int, default=200)
a = ap.parse_args()
files = sorted(glob.glob(os.path.join(a.passc, 'acc', 'acc_shard*of*.npz')))
assert files, 'no accumulators'
tot = None; done = []
for f in files:
    z = np.load(f, allow_pickle=False)
    meta = json.loads(str(z['meta'])); done += z['done'].tolist()
    if tot is None:
        tot = {k: z[k].astype(np.float64) for k in ('wf', 'w', 'wf2', 'n', 'ncell')}
    else:
        for k in tot:
            tot[k] += z[k]
print(f'{len(files)} shard accumulators, {len(set(done))} files')
cells = pd.read_csv(meta['cells']); cells = cells[cells.N >= meta['acc_min_n']].reset_index(drop=True)
cells['N_STACKED'] = tot['ncell'].astype(np.int64)
with np.errstate(all='ignore'):
    mean = np.where(tot['w'] > 0, tot['wf'] / tot['w'], np.nan)
    var = np.where(tot['w'] > 0, tot['wf2'] / tot['w'] - mean**2, np.nan)
    scat = np.sqrt(np.clip(var, 0, None))
    snr = np.nanmedian(np.where(tot['w'] > 0, mean * np.sqrt(tot['w']), np.nan), axis=1)
cells['STACK_SNR_MED'] = snr
wave = 10.0**((meta['rid_min'] + np.arange(meta['nrest'])) * meta['dlog'])
od = os.path.join(a.passc, 'stacks'); os.makedirs(od, exist_ok=True)
keep = cells.N_STACKED >= a.floor
for zb in sorted(cells[keep].ZB.unique()):
    r = np.where(keep & (cells.ZB == zb))[0]
    fn = os.path.join(od, f'stacks_z{zb:02d}.fits')
    with fitsio.FITS(fn + '.tmp', 'rw', clobber=True) as fo:
        fo.write(mean[r].astype(np.float32), extname='FLUX'); fo.write(tot['w'][r].astype(np.float32), extname='IVAR')
        fo.write(scat[r].astype(np.float32), extname='SCATTER'); fo.write(tot['n'][r].astype(np.int32), extname='NSPEC')
        fo.write(wave, extname='WAVE'); fo.write(cells.iloc[r].to_records(index=False), extname='CELLS')
    os.replace(fn + '.tmp', fn)
cells.to_csv(os.path.join(od, 'cell_stacks.csv'), index=False)
print(f'stacks written: {keep.sum()} cells (floor {a.floor}) in {cells[keep].ZB.nunique()} z slices; '
      f'median stack S/N {np.nanmedian(snr[keep]):.1f}; spectra stacked {int(tot["ncell"].sum()):,}')
