#!/usr/bin/env python3
"""Build the compact redshift lookup used by pass B (one-time, ~1 min).

Input : ~/data/sdssdesi_cache/zlookup_main.parquet  (main survey rows of zall-pix-loa, built from the zall cache)
Output: ~/data/sdssdesi_cache/zlookup_compact.npz with
        TARGETID int64, Z float32, ZERR float32, ZWARN int32, DELTACHI2 float32, ZCAT_PRIMARY bool,
        SPECTYPE_CODE int8 (0 GALAXY, 1 QSO, 2 STAR, 3 other), rows sorted by (program, healpix),
        and OFFSETS: for each (program code 0=dark 1=bright, healpix) the [start, stop) row range.
A worker loads this as numpy arrays (~1.6 GB, no conversion transient) and slices by offsets.
"""
import os, time
import numpy as np
import pyarrow.parquet as pq

CACHE = os.path.expanduser('~/data/sdssdesi_cache')
t0 = time.time()
t = pq.read_table(f'{CACHE}/zlookup_main.parquet', columns=['TARGETID', 'PROGRAM', 'HEALPIX', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'ZCAT_PRIMARY'])
prog = t.column('PROGRAM').to_numpy(zero_copy_only=False)
pcode = np.where(prog == 'dark', 0, np.where(prog == 'bright', 1, 2)).astype(np.int8)
hpx = t.column('HEALPIX').to_numpy().astype(np.int32)
stype = t.column('SPECTYPE').to_numpy(zero_copy_only=False)
scode = np.full(len(stype), 3, np.int8)
for k, s in enumerate(['GALAXY', 'QSO', 'STAR']):
    scode[stype == s] = k
order = np.lexsort((hpx, pcode))
arrays = dict(
    TARGETID=t.column('TARGETID').to_numpy()[order].astype(np.int64),
    Z=t.column('Z').to_numpy()[order].astype(np.float32),
    ZERR=t.column('ZERR').to_numpy()[order].astype(np.float32),
    ZWARN=t.column('ZWARN').to_numpy()[order].astype(np.int32),
    DELTACHI2=t.column('DELTACHI2').to_numpy()[order].astype(np.float32),
    ZCAT_PRIMARY=t.column('ZCAT_PRIMARY').to_numpy(zero_copy_only=False)[order].astype(bool),
    SPECTYPE_CODE=scode[order], PCODE=pcode[order], HEALPIX=hpx[order])
key = arrays['PCODE'].astype(np.int64) * 10**6 + arrays['HEALPIX']
uniq, start = np.unique(key, return_index=True)
stop = np.append(start[1:], len(key))
arrays['OFF_KEY'] = uniq; arrays['OFF_START'] = start.astype(np.int64); arrays['OFF_STOP'] = stop.astype(np.int64)
np.savez(f'{CACHE}/zlookup_compact.npz', **arrays)
print(f'{len(key):,} rows, {len(uniq):,} (program,healpix) groups, {os.path.getsize(f"{CACHE}/zlookup_compact.npz")/1e9:.2f} GB, {time.time()-t0:.0f}s')
