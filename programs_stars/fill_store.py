#!/usr/bin/env python3
"""Fill star_hook's store for the coadds the quasar-stack read does not visit (or finished before the hook was in): one
niced reader, priority order, same arrays as the hook. Resumable (a coadd already in the store is skipped).

Usage: python3 fill_store.py --list FILES.txt      # coadd basenames, in the order to read
"""
import os, sys, time, argparse, signal
import numpy as np, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import star_hook, starlib
passB = starlib.passB
STOP = {'flag': False}
signal.signal(signal.SIGTERM, lambda *a: STOP.update(flag=True))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--list', required=True)
    a = ap.parse_args(); t0 = time.time()
    SH = star_hook.StarHook()
    todo = [l.strip() for l in open(a.list) if l.strip()]
    todo = [f for f in todo if SH.needs(f)]
    print(f'{len(todo):,} coadds to read', flush=True)
    for k, fname in enumerate(todo):
        if STOP['flag']:
            break
        if not SH.needs(fname):
            continue
        prog = 'dark' if '-dark-' in fname else 'bright'; hpx = fname.split('-')[-1].replace('.fits', '')
        p = passB.local_path(f'/Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main/{prog}/{int(hpx)//100}/{hpx}/{fname}')
        if p is None:
            print('MISSING', fname, flush=True); continue
        F = fitsio.FITS(p); tid = F['FIBERMAP'].read(columns=['TARGETID'])['TARGETID']
        w, flux, ivar, mask, _ = passB.combine_cameras(F); F.close()
        SH.process(fname, tid, w, flux, ivar, mask)
        if (k + 1) % 50 == 0:
            print(f'{k+1:,}/{len(todo):,} {fname} {time.time()-t0:.0f}s', flush=True)
    SH.gauge()
    print(f"{'stopped' if STOP['flag'] else 'done'}: {SH.nfile:,} coadds, {SH.nstar:,} stars, {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    main()
