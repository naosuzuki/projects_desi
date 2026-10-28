#!/usr/bin/env python3
"""Copy the step-6 stars of coadds already in star_hook's exdisk3 ARCHIVE into the internal STORE (for coadds archived before the
hook wrote the internal store, 2026-10-28). Resumable; loops until no archived coadd is left unconverted."""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import star_hook
SH = star_hook.StarHook()
t0 = time.time(); n = 0
while True:
    todo = [e.name for e in os.scandir(SH.archive) if e.name.endswith('.npz') and '.tmp' not in e.name
            and e.name[:-4] + '.fits' in SH.by_file and not os.path.exists(os.path.join(SH.store, e.name))]
    if not todo:
        break
    for f in todo:
        z = np.load(os.path.join(SH.archive, f)); keep = np.isin(z['TARGETID'], SH.by_file[f[:-4] + '.fits'])
        p = os.path.join(SH.store, f); tmp = p + '.tmp.npz'
        np.savez(tmp, TARGETID=z['TARGETID'][keep], RID0=z['RID0'], FLUX=z['FLUX'][keep], IVAR=z['IVAR'][keep]); os.replace(tmp, p)
        n += 1
        if n % 100 == 0:
            print(f'{n:,} converted, {time.time()-t0:.0f}s', flush=True); SH.gauge()
SH.gauge(); print(f'done: {n:,} converted, {time.time()-t0:.0f}s')
