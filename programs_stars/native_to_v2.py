#!/usr/bin/env python3
"""Build star_hook's v2 store (rest frame, one resampling at the pipeline RV) from the NATIVE spectra that the shared v2
walk (programs_quasars/run_walk_v2.sh, native_hook.py) writes to exdisk3: /Volumes/exdisk3/sdssdesi_native/desi/<coadd>.npz
with TARGETID, KIND, W0 = 3600, DW = 0.8, FLUX/IVAR [n, 7781] (combine_cameras output; ivar = 0 where masked).
No exdisk2 read. Follows the walk: converts each needed coadd as soon as its npz exists (the walk writes atomically),
until --until-done sees all walk shards finished and nothing left. Resumable (stored coadds are skipped).

Usage: python3 native_to_v2.py [--follow] [--dlog 1e-4]
"""
import os, sys, time, glob, argparse, signal
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import star_hook
NATIVE = '/Volumes/exdisk3/sdssdesi_native/desi'
STOP = {'flag': False}
signal.signal(signal.SIGTERM, lambda *a: STOP.update(flag=True))


def walk_running():
    import subprocess
    return subprocess.run(['pgrep', '-f', '02_passQ_desi.py --all'], capture_output=True).returncode == 0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--follow', action='store_true', help='keep following the walk until it ends')
    a = ap.parse_args(); t0 = time.time()
    SH = star_hook.StarHook(archive='')            # internal v2 store only: the native store on exdisk3 is the archive
    todo = sorted(SH.by_file); n = 0
    print(f'{len(todo):,} coadds hold step-6 stars; v2 store {SH.store}', flush=True)
    while not STOP['flag']:
        pend = [f for f in todo if SH.needs(f)]
        ready = [f for f in pend if os.path.exists(os.path.join(NATIVE, f.replace('.fits', '.npz')))]
        for f in ready:
            if STOP['flag']:
                break
            z = np.load(os.path.join(NATIVE, f.replace('.fits', '.npz')))
            fl = z['FLUX']; iv = z['IVAR']; w = float(z['W0']) + float(z['DW']) * np.arange(fl.shape[1])
            SH.process(f, z['TARGETID'], w, fl, iv, (iv <= 0).astype(np.int8)); n += 1
            if n % 200 == 0:
                print(f'  {n:,} converted, {len(pend) - len(ready):,} waiting for the walk; {time.time()-t0:.0f}s', flush=True)
        if not a.follow or (not walk_running() and not [f for f in todo if SH.needs(f) and os.path.exists(os.path.join(NATIVE, f.replace('.fits', '.npz')))]):
            break
        time.sleep(120)
    SH.gauge()
    left = [f for f in todo if SH.needs(f)]
    print(f"{'stopped' if STOP['flag'] else 'done'}: {n:,} converted this run; {len(todo)-len(left):,}/{len(todo):,} stored; {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    main()
