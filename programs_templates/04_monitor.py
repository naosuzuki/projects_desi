#!/usr/bin/env python3
"""Progress / resource monitor for the SDSS-DESI pipeline, with a resource guard.

Every INTERVAL seconds:
  * progress of each step (files written or a step's progress.json) -> rate, ETA
  * worker processes (CPU %, RSS, state) found by command-line pattern
  * system CPU (total, per core), load, memory, swap, disk free, disk read rate
  * GUARD: pauses workers (SIGSTOP) when memory gets low or swap grows, resumes (SIGCONT)
    when it recovers; pauses the highest shard when CPU stays saturated; all actions logged
  * writes ~/data/sdssdesi_status/status.json and serves the directory on http://127.0.0.1:PORT
    (monitor.html next to it is the gauge page)

Usage:  python3 04_monitor.py [--port 8765] [--interval 10] [--no-guard]
Step progress protocol (for steps other than pass B): write JSON to
~/data/sdssdesi_status/step_<name>.json  {"done": n, "total": N, "status": "running|done|idle", "note": "..."}
"""
import os, sys, time, json, glob, signal, argparse, threading, shutil, socket
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
import psutil

HOME = os.path.expanduser('~')
STATUS_DIR = f'{HOME}/data/sdssdesi_status'
HERE = os.path.dirname(os.path.abspath(__file__))

STEPS = [  # name, kind, source, total, worker pattern
    dict(name='Pass B: DESI measure', kind='files', pattern=f'{HOME}/data/sdssdesi_passB/*/*.parquet', total=34508, proc='02_passB_measure'),
    dict(name='SDSS pass', kind='json', path=f'{STATUS_DIR}/step_sdss.json', proc='05_passB_sdss'),
    dict(name='Cell assignment', kind='json', path=f'{STATUS_DIR}/step_assign.json', proc='06_assign_cells'),
    dict(name='Pass C: stacks', kind='json', path=f'{STATUS_DIR}/step_passC.json', proc='07_passC_stack'),
]
DISKS = {'internal (output)': HOME, 'exdisk2 (coadds)': '/Volumes/exdisk2', 'exdisk1 / My Book (overflow)': ['/Volumes/exdisk1', '/Volumes/My Book']}

# guard thresholds
MEM_PAUSE_GB = 8.0      # pause all workers when available memory falls below this
MEM_RESUME_GB = 16.0    # resume when available memory is back above this
SWAP_PAUSE_GB = 6.0     # or when swap in use exceeds this
CPU_PAUSE_PCT = 92.0    # pause the highest shard when total CPU stays above this ...
CPU_PAUSE_SEC = 90      # ... for this long
CPU_RESUME_PCT = 70.0

events = []
paused = {}             # pid -> reason
cpu_hot_since = None
io_prev = None

def log(msg):
    events.append(dict(t=time.strftime('%Y-%m-%d %H:%M:%S'), msg=msg))
    del events[:-200]
    print(time.strftime('%H:%M:%S'), msg, flush=True)

def step_progress(step, hist):
    now = time.time()
    if step['kind'] == 'files':
        files = glob.glob(step['pattern'])
        done = len(files); total = step['total']
        mt = sorted(os.path.getmtime(f) for f in files) if files else []
        started = mt[0] if mt else None
        recent = [t for t in mt if t > now - 1800]
        window_min = min(30.0, max((now - started) / 60, 0.5)) if started else 30.0     # the job may be younger than the window
        rate = len(recent) / window_min if len(recent) >= 3 else (done / max((mt[-1] - mt[0]) / 60, 1e-9) if len(mt) > 10 else 0.0)  # files/min
        status = 'done' if done >= total else ('running' if recent and recent[-1] > now - 600 else ('stalled' if done else 'idle'))
        last = mt[-1] if mt else None
    else:
        p = step.get('path')
        if p and os.path.exists(p):
            try:
                j = json.load(open(p))
            except Exception:
                j = {}
            done, total = j.get('done', 0), j.get('total', 0); status = j.get('status', 'idle'); started = j.get('started'); last = os.path.getmtime(p)
            h = hist.setdefault(step['name'], [])
            h.append((now, done)); del h[:-180]
            old = [x for x in h if x[0] < now - 1800] or h[:1]
            rate = (done - old[0][1]) / max((now - old[0][0]) / 60, 1e-9) if len(h) > 2 else 0.0
            if status not in ('done',) and total and done >= total: status = 'done'
        else:
            done, total, status, started, last, rate = 0, 0, 'idle', None, None, 0.0
    remaining = max(total - done, 0)
    eta_min = remaining / rate if rate > 0 else None
    return dict(name=step['name'], done=done, total=total, pct=(100.0 * done / total if total else 0.0), status=status,
                rate_per_min=rate, eta_min=eta_min, eta_at=(time.strftime('%a %H:%M', time.localtime(now + eta_min * 60)) if eta_min else None),
                started=(time.strftime('%Y-%m-%d %H:%M', time.localtime(started)) if started else None),
                elapsed_min=((now - started) / 60 if started else None), last_update=(time.strftime('%H:%M:%S', time.localtime(last)) if last else None))

def find_workers():
    out = []
    for p in psutil.process_iter(['pid', 'cmdline', 'create_time']):
        try:
            cmd = ' '.join(p.info['cmdline'] or [])
            for s in STEPS:
                if s['proc'] in cmd and 'python' in cmd:
                    shard = cmd.split('--shard')[1].split()[0] if '--shard' in cmd else '-'
                    with p.oneshot():
                        out.append(dict(pid=p.pid, step=s['name'], shard=shard, cpu=p.cpu_percent(interval=None), rss_gb=p.memory_info().rss / 1e9,
                                        state=('paused' if p.pid in paused else p.status()), age_min=(time.time() - p.info['create_time']) / 60))
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
    now = time.time()
    ncpu = psutil.cpu_count() or 1
    worker_share = sum(w['cpu'] for w in running) / (100.0 * ncpu) * 100.0     # % of the machine used by our workers
    # pause for CPU only when our own workers are the load (another user's job saturating the CPU is not our problem:
    # the workers are niced and the scheduler handles it)
    if cpu_total > CPU_PAUSE_PCT and worker_share > 40.0:
        cpu_hot_since = cpu_hot_since or now
        if now - cpu_hot_since > CPU_PAUSE_SEC and len(running) > 1:
            w = sorted(running, key=lambda x: x['shard'])[-1]
            try:
                os.kill(w['pid'], signal.SIGSTOP); paused[w['pid']] = f'cpu: {cpu_total:.0f}%'
                log(f'GUARD: paused pid {w["pid"]} (shard {w["shard"]}) - CPU {cpu_total:.0f}% for {CPU_PAUSE_SEC}s')
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
    return dict(
        time=time.strftime('%Y-%m-%d %H:%M:%S'), host=socket.gethostname(), interval=INTERVAL,
        steps=[step_progress(s, hist) for s in STEPS],
        workers=workers, n_workers=len(workers), n_paused=len(paused),
        system=dict(cpu_total=cpu_total, cpu_per_core=cpu_pc, ncpu=len(cpu_pc), load1=load[0], load5=load[1],
                    mem_total_gb=mem.total / 1e9, mem_used_gb=(mem.total - mem.available) / 1e9, mem_available_gb=mem.available / 1e9,
                    mem_pct=100.0 * (mem.total - mem.available) / mem.total, swap_used_gb=swap.used / 1e9, disk_read_mbs=read_mbs),
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
    ap.add_argument('--passb-dir', default=None, help='watch another pass-B output dir (e.g. the 20-file test)')
    a = ap.parse_args(); INTERVAL = a.interval
    if a.passb_dir:
        STEPS[0]['pattern'] = os.path.join(os.path.expanduser(a.passb_dir), '*', '*.parquet')
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
            print(json.dumps(st, indent=1)[:3000]); break
        time.sleep(INTERVAL)

if __name__ == '__main__':
    main()
