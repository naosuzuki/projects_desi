#!/usr/bin/env python3
"""Pipeline monitor for the SDSS-DESI template project: one page, every step, with ETA, plus CPU / memory / disk and a
resource guard. Supersedes programs_templates/04_monitor.py (same protocol, same status dir, same port).

Every INTERVAL seconds:
  * progress of each step -> done/total, rate, ETA, finish time; workers of the step (CPU %, RSS, state)
      kind 'files': count output files matching a glob against a known total (passes B, Q-SDSS, Q-DESI, ...)
      kind 'json' : a step writes ~/data/sdssdesi_status/step_<name>.json
                    {"done": n, "total": N, "status": "running|done|idle", "note": "...", "started": unix_time}
  * system: CPU total and per core, load, memory, swap, disk free, disk read rate
  * GUARD (unchanged from 04_monitor.py): SIGSTOP all workers when available memory < MEM_PAUSE_GB or swap >
    SWAP_PAUSE_GB, SIGCONT at MEM_RESUME_GB; SIGSTOP the highest shard when CPU > CPU_PAUSE_PCT for CPU_PAUSE_SEC
    while our workers are the load; all actions logged
  * writes ~/data/sdssdesi_status/status.json and serves that directory: http://127.0.0.1:PORT/monitor.html

Usage:  python3 monitor.py [--port 8765] [--interval 10] [--no-guard] [--once]
Run detached:  nohup python3 monitor.py --port 8765 > ~/data/sdssdesi_status/monitor.log 2>&1 &
Only one guard may run at a time (two guards would fight over the same workers).
"""
import os, sys, time, json, glob, signal, argparse, threading, shutil, socket
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
import psutil

HOME = os.path.expanduser('~')
STATUS_DIR = f'{HOME}/data/sdssdesi_status'
HERE = os.path.dirname(os.path.abspath(__file__))
N_DESI_COADDS = 34508
N_SDSS_PLATES = 4131 + 3953          # DR17 v5_13_2 + DR8 (26/103/104) spPlate files on disk
# The shared DESI walk (02_passQ_desi.py --with-stars --expfibermap --with-galaxies) writes several outputs ('walk=True'
# below). After a restart it re-reads coadds that some outputs already have, so those outputs make no files for hours
# while the walk is busy: they show as 'waiting' with the walk's ETA, not 'stalled'.
WALK_PROC = '02_passQ_desi'

# The pipeline, by track. 'proc' is the substring of the worker command line (python scripts), used for the worker
# table, the per-step CPU/RSS and the guard. Steps that do not exist yet use the JSON protocol and show as 'idle'.
STEPS = [
    # --- common
    dict(track='common', name='Pass B: DESI per-spectrum measurements (all classes)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passB/*/*.parquet', total=N_DESI_COADDS, proc='02_passB_measure'),
    dict(track='common', name='v2 shared DESI walk: one read for galaxies, quasars, stars (coadds with native quasar/star spectra stored)', kind='files',
         pattern='/Volumes/exdisk3/sdssdesi_native/desi/*.npz', total=34450, proc='02_passQ_desi.py --all --with-native', walk=True),
    # --- galaxies
    dict(track='galaxies', name='Galaxy cell assignment', kind='json', path=f'{STATUS_DIR}/step_gal_assign.json', proc='07_assign_cells'),
    dict(track='galaxies', name='v2 SDSS galaxy pass DR8 (spPlates, sub-pixel; outputs on exdisk3)', kind='json', path=f'{STATUS_DIR}/step_gal_sdss_dr8_shard0of1.json', proc='11_passG_sdss.py --survey dr8'),
    dict(track='galaxies', name='v2 SDSS galaxy pass DR17 (spPlates, sub-pixel; outputs on exdisk3)', kind='json', path=f'{STATUS_DIR}/step_gal_sdss_dr17_shard0of1.json', proc='11_passG_sdss.py --survey dr17'),
    dict(track='galaxies', name='v1 Pass C: galaxy stacks (DESI coadds, whole-pixel shift)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passC/tables/*/*.parquet', total=N_DESI_COADDS, proc='08_passC_stack'),
    dict(track='galaxies', name='v2 Pass C: galaxy stacks (galaxy hook in the v2 shared DESI walk, sub-pixel)', kind='files',
         pattern='/Volumes/exdisk3/sdssdesi_passC_v2/tables/*/*.parquet', total=N_DESI_COADDS, proc='02_passQ_desi.py --all --with-native', walk=True),
    dict(track='galaxies', name='Galaxy EM-PCA (latest run; v2 after the walk)', kind='json', path=f'{STATUS_DIR}/step_gal_empca.json', proc='10_galaxy_empca'),
    # --- quasars
    dict(track='quasars', name='Pass Q: SDSS quasars (DR17 + DR8 spPlates)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passQ_sdss/*/*.parquet', total=N_SDSS_PLATES, proc='01_passQ_sdss'),
    dict(track='quasars', name='Pass Q: DESI quasars (coadds)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passQ/*/*.parquet', total=N_DESI_COADDS, proc='02_passQ_desi', walk=True),
    dict(track='quasars', name='Quasar master table, SDSS (DR16Q match, checks)', kind='json', path=f'{STATUS_DIR}/step_qso_table_sdss.json', proc='03_qso_table_sdss'),
    dict(track='quasars', name='Quasar master table, DESI', kind='json', path=f'{STATUS_DIR}/step_qso_table_desi.json', proc='03_qso_table_desi'),
    dict(track='quasars', name='Quasar cell table (6 axes)', kind='json', path=f'{STATUS_DIR}/step_qso_cells.json', proc='04_cells'),
    dict(track='quasars', name='Quasar stacks v1 (whole-pixel shift)', kind='json', path=f'{STATUS_DIR}/step_qso_stack.json', proc='05_stack.py --survey'),
    dict(track='quasars', name='Quasar EM-PCA v1', kind='json', path=f'{STATUS_DIR}/step_qso_empca.json', proc='06_empca.py --stacks /Users/suzuki/data/sdssdesi_products/quasar/qso_'),
    dict(track='quasars', name='Quasar stacks v2 (sub-pixel; SDSS spPlates + DESI native store, after the v2 walk)', kind='json', path=f'{STATUS_DIR}/step_qso_stack_v2.json', proc='05_stack.py --v2'),
    dict(track='quasars', name='Quasar EM-PCA v2', kind='json', path=f'{STATUS_DIR}/step_qso_empca_v2.json', proc='quasar/v2/qso_stacks'),
    # --- stars (Nao's steps, 2026-10-24 night; scripts in programs_stars/, runbook there; detail page http://127.0.0.1:8767/monitor.html)
    dict(track='stars', name='Stars 1: catalogs with Gaia parallax (DESI lookup, DESI HR, SDSS)', kind='json', path=f'{STATUS_DIR}/step_star_catalog.json', proc='step1_'),
    dict(track='stars', name='Stars 2-3: SDSS pass (PSF/Gaia calibration, RV frame, indices)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passS_sdss/*/*.parquet', total=N_SDSS_PLATES, proc='step2_3_pass_sdss'),
    dict(track='stars', name='Stars 2b: joint spectrophotometric fit (Gaia scale, SDSS/LS tilt, zero points)', kind='json', path=f'{STATUS_DIR}/step_star_spcal.json', proc='step2b_spcal'),
    dict(track='stars', name='Stars 3b: 3D extinction at the parallax distance (Edenhofer + Bayestar)', kind='json', path=f'{STATUS_DIR}/step_star_ebv3d.json', proc='step3b_ebv3d'),
    dict(track='stars', name='Stars 2-3: DESI pass (coadds, after pass B)', kind='files',
         pattern=f'{HOME}/data/sdssdesi_passS/*/*.parquet', total=N_DESI_COADDS, proc='step2_3_pass_desi', walk=True),
    dict(track='stars', name='Stars 5: HR grid, same-cell similarity test (SDSS, DESI)', kind='json', path=f'{STATUS_DIR}/step_star_cells.json', proc='step5_hr_cells'),
    dict(track='stars', name='Stars 6: stacks of the alike cells (SDSS, DESI)', kind='json', path=f'{STATUS_DIR}/step_star_stack.json', proc='step6_stack'),
    dict(track='stars', name='Stars 7: EM-PCA', kind='json', path=f'{STATUS_DIR}/step_star_empca.json', proc='step7_empca'),
    dict(track='stars', name='Stars 8: Calibration 2 (DESI response per exposure and petal, per coadd)', kind='json', path=f'{STATUS_DIR}/step_star_response.json', proc='step8_response'),
]
DISKS = {'internal (output)': HOME, 'exdisk2 (coadds)': '/Volumes/exdisk2', 'exdisk1 / My Book (overflow)': ['/Volumes/exdisk1', '/Volumes/My Book']}

MEM_PAUSE_GB = 8.0; MEM_RESUME_GB = 16.0; SWAP_PAUSE_GB = 6.0
CPU_PAUSE_PCT = 92.0; CPU_PAUSE_SEC = 90; CPU_RESUME_PCT = 70.0

events = []; paused = {}; cpu_hot_since = None; io_prev = None
INTERVAL = 10.0


def log(msg):
    events.append(dict(t=time.strftime('%Y-%m-%d %H:%M:%S'), msg=msg)); del events[:-200]
    print(time.strftime('%H:%M:%S'), msg, flush=True)


def step_progress(step, hist, workers):
    now = time.time()
    note = ''
    if step['kind'] == 'files':
        files = glob.glob(step['pattern'])
        done = len(files); total = step['total']
        mt = sorted(os.path.getmtime(f) for f in files) if files else []
        started = mt[0] if mt else None
        recent = [t for t in mt if t > now - 1800]
        window_min = min(30.0, max((now - started) / 60, 0.5)) if started else 30.0
        rate = len(recent) / window_min if len(recent) >= 3 else (done / max((mt[-1] - mt[0]) / 60, 1e-9) if len(mt) > 10 else 0.0)
        status = 'done' if done >= total else ('running' if recent and recent[-1] > now - 600 else ('stalled' if done else 'idle'))
        last = mt[-1] if mt else None
    else:
        p = step.get('path'); j = {}
        if p and os.path.exists(p):
            try:
                j = json.load(open(p))
            except Exception:
                j = {}
            done, total = j.get('done', 0), j.get('total', 0); status = j.get('status', 'idle'); started = j.get('started'); note = j.get('note', '')
            last = os.path.getmtime(p)
            h = hist.setdefault(step['name'], []); h.append((now, done)); del h[:-180]
            old = [x for x in h if x[0] < now - 1800] or h[:1]
            rate = (done - old[0][1]) / max((now - old[0][0]) / 60, 1e-9) if len(h) > 2 else 0.0
            if status != 'done' and total and done >= total: status = 'done'
            if status == 'running' and now - last > 900: status = 'stalled'
        else:
            done, total, status, started, last, rate = 0, 0, 'idle', None, None, 0.0
    ws = [w for w in workers if w['step'] == step['name'] or (step.get('walk') and w.get('walk'))]
    if status == 'idle' and ws: status = 'running'
    remaining = max(total - done, 0)
    eta_min = remaining / rate if (rate > 0 and status != 'done') else (0.0 if status == 'done' else None)
    if step['kind'] == 'json' and status != 'done' and j.get('eta_epoch'):      # the step knows better (e.g. Stars 6: its inputs arrive in a known read order)
        eta_min = max((j['eta_epoch'] - now) / 60, 0.0)
    return dict(track=step['track'], name=step['name'], done=done, total=total, pct=(100.0 * done / total if total else 0.0), status=status,
                rate_per_min=rate, eta_min=eta_min, eta_at=(time.strftime('%a %H:%M', time.localtime(now + eta_min * 60)) if eta_min else None),
                started=(time.strftime('%m-%d %H:%M', time.localtime(started)) if started else None),
                elapsed_min=((now - started) / 60 if started and status != 'done' else ((last - started) / 60 if started and last else None)),
                last_update=(time.strftime('%H:%M:%S', time.localtime(last)) if last else None), note=note,
                n_workers=len(ws), workers_cpu=sum(w['cpu'] for w in ws), workers_rss_gb=sum(w['rss_gb'] for w in ws),
                n_paused=sum(1 for w in ws if w['state'] == 'paused'))


def find_workers():
    out = []
    for p in psutil.process_iter(['pid', 'cmdline', 'create_time']):
        try:
            cmd = ' '.join(p.info['cmdline'] or [])
            if 'python' not in cmd or 'monitor.py' in cmd:
                continue
            for s in STEPS:
                if s['proc'] in cmd:
                    shard = cmd.split('--shard')[1].split()[0] if '--shard' in cmd else '-'
                    with p.oneshot():
                        out.append(dict(pid=p.pid, step=s['name'], track=s['track'], shard=shard, cpu=p.cpu_percent(interval=None),
                                        rss_gb=p.memory_info().rss / 1e9, state=('paused' if p.pid in paused else p.status()),
                                        age_min=(time.time() - p.info['create_time']) / 60, walk=WALK_PROC in cmd))
                    break
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return out


def guard(workers, mem, swap, cpu_total):
    global cpu_hot_since
    avail_gb = mem.available / 1e9; swap_gb = swap.used / 1e9
    running = [w for w in workers if w['pid'] not in paused]
    if (avail_gb < MEM_PAUSE_GB or swap_gb > SWAP_PAUSE_GB) and running:
        for w in running:
            try:
                os.kill(w['pid'], signal.SIGSTOP); paused[w['pid']] = f'memory: {avail_gb:.1f} GB available, swap {swap_gb:.1f} GB'
            except ProcessLookupError:
                pass
        log(f'GUARD: paused {len(running)} worker(s) - available memory {avail_gb:.1f} GB, swap {swap_gb:.1f} GB')
    elif paused and avail_gb > MEM_RESUME_GB and swap_gb < SWAP_PAUSE_GB:
        mem_paused = [pid for pid, r in paused.items() if r.startswith('memory')]
        for pid in mem_paused:
            try:
                os.kill(pid, signal.SIGCONT)
            except ProcessLookupError:
                pass
            paused.pop(pid, None)
        if mem_paused:
            log(f'GUARD: resumed {len(mem_paused)} worker(s) - available memory {avail_gb:.1f} GB')
    now = time.time(); ncpu = psutil.cpu_count() or 1
    worker_share = sum(w['cpu'] for w in running) / (100.0 * ncpu) * 100.0
    if cpu_total > CPU_PAUSE_PCT and worker_share > 40.0:
        cpu_hot_since = cpu_hot_since or now
        if now - cpu_hot_since > CPU_PAUSE_SEC and len(running) > 1:
            w = sorted(running, key=lambda x: (x['step'], x['shard']))[-1]
            try:
                os.kill(w['pid'], signal.SIGSTOP); paused[w['pid']] = f'cpu: {cpu_total:.0f}%'
                log(f'GUARD: paused pid {w["pid"]} ({w["step"]}, shard {w["shard"]}) - CPU {cpu_total:.0f}% for {CPU_PAUSE_SEC}s')
            except ProcessLookupError:
                pass
            cpu_hot_since = now
    else:
        cpu_hot_since = None
        if cpu_total < CPU_RESUME_PCT:
            for pid, r in list(paused.items()):
                if r.startswith('cpu'):
                    try:
                        os.kill(pid, signal.SIGCONT)
                    except ProcessLookupError:
                        pass
                    paused.pop(pid, None); log(f'GUARD: resumed pid {pid} - CPU {cpu_total:.0f}%')
    for pid in list(paused):
        if not psutil.pid_exists(pid):
            paused.pop(pid)


def collect(hist, do_guard):
    global io_prev
    now = time.time()
    cpu_pc = psutil.cpu_percent(interval=None, percpu=True); cpu_total = sum(cpu_pc) / len(cpu_pc)
    mem = psutil.virtual_memory(); swap = psutil.swap_memory(); load = os.getloadavg()
    io = psutil.disk_io_counters()
    read_mbs = ((io.read_bytes - io_prev[0]) / max(now - io_prev[1], 1e-9) / 1e6) if io_prev else 0.0
    io_prev = (io.read_bytes, now)
    disks = {}
    for k, p in DISKS.items():
        cands = p if isinstance(p, list) else [p]
        disks[k] = dict(free_gb=None, total_gb=None, used_pct=None, missing=True)
        for c in cands:
            if os.path.ismount(c) or c == HOME:
                u = shutil.disk_usage(c); disks[k] = dict(free_gb=u.free / 1e9, total_gb=u.total / 1e9, used_pct=100.0 * u.used / u.total, path=c)
                break
    workers = find_workers()
    if do_guard:
        guard(workers, mem, swap, cpu_total)
    steps = [step_progress(s, hist, workers) for s in STEPS]
    # walk outputs: one that is behind while another advances and the walk is alive is waiting for the walk to reach it
    walk = [st for st, s in zip(steps, STEPS) if s.get('walk')]
    lead = [st for st in walk if st['status'] == 'running' and st['eta_min']]
    if lead and any(w.get('walk') for w in workers):
        lead = max(lead, key=lambda x: x['rate_per_min'])     # the output the walk is producing now
        for st in walk:      # stalled, or trickling (a few retried coadds) far below the walk's rate
            if st is not lead and (st['status'] == 'stalled' or (st['status'] == 'running' and st['rate_per_min'] < 0.25 * lead['rate_per_min'])):
                st.update(status='waiting', rate_per_min=0.0, eta_min=lead['eta_min'], eta_at=lead['eta_at'],
                          note=f"walk alive; resumes when the walk passes coadd {st['done']:,} (now at {lead['done']:,} for "
                               f"'{lead['name'].split(':')[0]}'); ETA = walk ETA, a lower limit")
    # pipeline summary per track: running step, next ETA
    tracks = {}
    for st in steps:
        t = tracks.setdefault(st['track'], dict(track=st['track'], n_steps=0, n_done=0, running=None, eta_at=None, eta_min=None))
        t['n_steps'] += 1; t['n_done'] += st['status'] == 'done'
        rank = {'running': 0, 'stalled': 1, 'waiting': 2}.get(st['status'], 9) + (0 if st['eta_min'] else 0.5)
        if rank < 9 and rank < t.get('_rank', 9):     # running first, then stalled, then waiting; with an ETA first
            t['_rank'] = rank; t['running'] = st['name']; t['eta_at'] = st['eta_at']; t['eta_min'] = st['eta_min']
    hist_sys = hist.setdefault('__sys__', []); hist_sys.append((now, cpu_total, (mem.total - mem.available) / 1e9)); del hist_sys[:-360]   # 1 h at 10 s
    return dict(
        time=time.strftime('%Y-%m-%d %H:%M:%S'), host=socket.gethostname(), interval=INTERVAL,
        steps=steps, tracks=[{k: v for k, v in t.items() if k != '_rank'} for t in tracks.values()], workers=workers, n_workers=len(workers), n_paused=len(paused),
        system=dict(cpu_total=cpu_total, cpu_per_core=cpu_pc, ncpu=len(cpu_pc), load1=load[0], load5=load[1],
                    mem_total_gb=mem.total / 1e9, mem_used_gb=(mem.total - mem.available) / 1e9, mem_available_gb=mem.available / 1e9,
                    mem_pct=100.0 * (mem.total - mem.available) / mem.total, swap_used_gb=swap.used / 1e9, disk_read_mbs=read_mbs,
                    history=[dict(t=t, cpu=c, mem=m) for t, c, m in hist_sys[::3]]),
        disks=disks, guard=dict(enabled=do_guard, mem_pause_gb=MEM_PAUSE_GB, mem_resume_gb=MEM_RESUME_GB, swap_pause_gb=SWAP_PAUSE_GB,
                                 cpu_pause_pct=CPU_PAUSE_PCT, cpu_resume_pct=CPU_RESUME_PCT, paused={str(k): v for k, v in paused.items()}),
        events=events[-30:][::-1])


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *a): pass
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store'); super().end_headers()


def serve(port):
    os.chdir(STATUS_DIR)
    httpd = ThreadingHTTPServer(('127.0.0.1', port), Quiet)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    log(f'serving {STATUS_DIR} on http://127.0.0.1:{port}/monitor.html')


def main():
    global INTERVAL
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=8765); ap.add_argument('--interval', type=float, default=10)
    ap.add_argument('--no-guard', action='store_true'); ap.add_argument('--once', action='store_true')
    a = ap.parse_args(); INTERVAL = a.interval
    os.makedirs(STATUS_DIR, exist_ok=True)
    shutil.copy(os.path.join(HERE, 'monitor.html'), os.path.join(STATUS_DIR, 'monitor.html'))
    psutil.cpu_percent(interval=None)
    for p in psutil.process_iter():
        try: p.cpu_percent(interval=None)
        except Exception: pass
    if not a.once:
        serve(a.port)
    hist = {}
    while True:
        try:
            st = collect(hist, not a.no_guard)
            tmp = os.path.join(STATUS_DIR, 'status.json.tmp')
            json.dump(st, open(tmp, 'w'), indent=1); os.replace(tmp, os.path.join(STATUS_DIR, 'status.json'))
        except Exception as e:
            log(f'monitor error: {e!r}')
        if a.once:
            print(json.dumps(st, indent=1)[:4000]); break
        time.sleep(INTERVAL)


if __name__ == '__main__':
    main()
