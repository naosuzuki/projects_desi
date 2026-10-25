#!/usr/bin/env python3
"""Steps monitor for the stars thread (gauge page with progress, rate and ETA per step), built on the pass-B monitor
(programs_templates/04_monitor.py: resources, workers, disks; the guard is off here - the star passes are light).

Steps shown: 1 catalogs (files present) · 2-3 SDSS pass (plates written / all plates) · 3b 3D extinction (progress
line of its log or step_ebv3d.json) · pass B (prerequisite of the DESI pass) · 2-3 DESI pass (coadds written / 34,508)
· 5 HR-cell check (cell tables) · 6 stacks (stack files). Status dir: ~/data/sdssdesi_status_stars (status.json,
monitor.html); page: http://127.0.0.1:8767/monitor.html (or preview "monitor-stars" from .claude/launch.json).
Usage: python3 monitor_stars.py [--port 8767] [--interval 10] [--once]
"""
import os, sys, re, glob, json, time, argparse, shutil, importlib.util
HOME = os.path.expanduser('~'); HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('monitor', os.path.join(HERE, '..', 'programs_templates', '04_monitor.py'))
M = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(M)

STATUS_DIR = f'{HOME}/data/sdssdesi_status_stars'
M.STATUS_DIR = STATUS_DIR
CACHE = f'{HOME}/data/sdssdesi_cache'


def n_plates():
    S = f'{HOME}/sdssredux'
    return sum(len(glob.glob(f'{d}/*/spPlate-*.fits')) for d in (f'{S}/v5_13_2', f'{S}/DR8/26', f'{S}/DR8/103', f'{S}/DR8/104'))


M.STEPS = [
    dict(name='1 catalogs (DESI lookup, DESI HR, SDSS+Gaia)', kind='exists', files=[f'{CACHE}/zall_stars_main.parquet', f'{CACHE}/gaia_desi_stars_hr.parquet', f'{CACHE}/sdss_stars_gaia.parquet'], proc='step1_'),
    dict(name='2-3 SDSS pass (spPlates)', kind='files', pattern=f'{HOME}/data/sdssdesi_passS_sdss/*/*.parquet', total=n_plates(), proc='step2_3_pass_sdss'),
    dict(name='2b joint spectrophotometric fit', kind='json', path=f'{HOME}/data/sdssdesi_status/step_star_spcal.json', proc='step2b_spcal'),
    dict(name='3b 3D extinction (Edenhofer + Bayestar)', kind='log', path=f'{HOME}/data/sdssdesi_status/ebv3d.log', json=f'{HOME}/data/sdssdesi_status/step_star_ebv3d.json', proc='step3b_ebv3d'),
    dict(name='pass B (prerequisite of the DESI pass)', kind='files', pattern=f'{HOME}/data/sdssdesi_passB/*/*.parquet', total=34508, proc='02_passB_measure'),
    dict(name='2-3 DESI pass (coadds)', kind='files', pattern=f'{HOME}/data/sdssdesi_passS/*/*.parquet', total=34508, proc='step2_3_pass_desi'),
    dict(name='5 HR-cell check (SDSS, DESI)', kind='exists', files=[f'{CACHE}/stars/cells_sdss_*.csv', f'{CACHE}/stars/cells_desi_*.csv'], proc='step5_hr_cells'),
    dict(name='6 stacks (SDSS, DESI)', kind='exists', files=[f'{HOME}/data/sdssdesi_stars/stacks_sdss_*.fits', f'{HOME}/data/sdssdesi_stars/stacks_desi_*.fits'], proc='step6_stack'),
]
M.DISKS = {'internal (output)': HOME, 'exdisk2 (coadds)': '/Volumes/exdisk2'}
_orig = M.step_progress
_LOGRE = re.compile(r'^\s*([\d,]+)/([\d,]+)\s*$', re.M)


def step_progress(step, hist):
    now = time.time()
    if step['kind'] == 'exists':
        hits = [glob.glob(f) for f in step['files']]
        done = sum(1 for h in hits if h); total = len(step['files'])
        mt = [max(os.path.getmtime(x) for x in h) for h in hits if h]
        return dict(name=step['name'], done=done, total=total, pct=100.0 * done / total, status=('done' if done == total else ('running' if done else 'idle')),
                    rate_per_min=0.0, eta_min=None, eta_at=None, started=None, elapsed_min=None,
                    last_update=(time.strftime('%H:%M:%S', time.localtime(max(mt))) if mt else None))
    if step['kind'] == 'log':
        j = step.get('json')
        if j and os.path.exists(j):
            return _orig(dict(step, kind='json', path=j), hist)
        p = step['path']
        if not os.path.exists(p):
            return dict(name=step['name'], done=0, total=0, pct=0.0, status='idle', rate_per_min=0.0, eta_min=None, eta_at=None, started=None, elapsed_min=None, last_update=None)
        txt = open(p).read()[-20000:]
        m = _LOGRE.findall(txt)
        done, total = (int(m[-1][0].replace(',', '')), int(m[-1][1].replace(',', ''))) if m else (0, 0)
        finished = 'written' in txt.split('\n')[-3:][0] if txt else False
        h = hist.setdefault(step['name'], []); h.append((now, done)); del h[:-180]
        old = [x for x in h if x[0] < now - 600] or h[:1]
        rate = (done - old[0][1]) / max((now - old[0][0]) / 60, 1e-9) if len(h) > 2 else 0.0
        last = os.path.getmtime(p); alive = last > now - 300
        status = 'done' if (total and done >= total and not alive) or ('written' in txt[-2000:] and not alive) else ('running' if alive else ('stalled' if done else 'idle'))
        eta = (total - done) / rate if rate > 0 and total > done else None
        return dict(name=step['name'], done=done, total=total, pct=(100.0 * done / total if total else 0.0), status=status, rate_per_min=rate, eta_min=eta,
                    eta_at=(time.strftime('%a %H:%M', time.localtime(now + eta * 60)) if eta else None), started=None, elapsed_min=None,
                    last_update=time.strftime('%H:%M:%S', time.localtime(last)))
    return _orig(step, hist)


M.step_progress = step_progress


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--port', type=int, default=8767); ap.add_argument('--interval', type=float, default=10); ap.add_argument('--once', action='store_true')
    a = ap.parse_args(); M.INTERVAL = a.interval
    os.makedirs(STATUS_DIR, exist_ok=True)
    shutil.copy(os.path.join(HERE, 'monitor_stars.html'), os.path.join(STATUS_DIR, 'monitor.html'))
    import psutil
    psutil.cpu_percent(interval=None)
    for p in psutil.process_iter():
        try: p.cpu_percent(interval=None)
        except Exception: pass
    if not a.once:
        M.serve(a.port)
    hist = {}
    while True:
        try:
            st = M.collect(hist, False)
            tmp = os.path.join(STATUS_DIR, 'status.json.tmp'); json.dump(st, open(tmp, 'w'), indent=1); os.replace(tmp, os.path.join(STATUS_DIR, 'status.json'))
        except Exception as e:
            M.log(f'monitor error: {e!r}')
        if a.once:
            print(json.dumps([{k: s[k] for k in ('name', 'done', 'total', 'pct', 'status', 'rate_per_min', 'eta_at')} for s in st['steps']], indent=1)); break
        time.sleep(a.interval)


if __name__ == '__main__':
    main()
