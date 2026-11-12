#!/usr/bin/env python3
"""igm_status.py — print igmdust pipeline status as JSON (runs on astro4).

Parses the SLURM job logs in ~/slurmlogs and the download logs; the gauge
page (index.html) renders the result. Add gauges here as jobs change.
"""
import json
import os
import re
import time

LOGS = os.path.expanduser("~/slurmlogs")


def read(path):
    try:
        with open(os.path.expanduser(path), errors="ignore") as f:
            return f.read()
    except OSError:
        return ""


def mtime(path):
    try:
        return os.path.getmtime(os.path.expanduser(path))
    except OSError:
        return 0


def fmt_eta(s):
    s = max(0, int(s))
    return f"{s // 3600}h{s % 3600 // 60:02d}m" if s >= 3600 \
        else f"{s // 60}m"


def frac(text, tag):
    """Last 'fg done X/Y ... Ts' after a [tag] marker, or None."""
    hits = re.findall(
        rf"\[{tag}\].*?fg done ([\d,]+)/([\d,]+).*?(\d+)s", text)
    if not hits:
        return None
    x, y = (int(v.replace(",", "")) for v in hits[-1][:2])
    return x, y, int(hits[-1][2])


def fmt_t(e):
    return time.strftime("%m-%d %H:%M", time.localtime(e))


def runtimes(elapsed, ts, done):
    """'started ~A' (log mtime minus in-log elapsed; catalog-loading time
    is not in the counter, so '~') plus ', ended B' once finished — the
    mtime of a finished log is the moment of its last write."""
    if not elapsed or not ts:
        return ""
    s = f"started ~{fmt_t(ts - elapsed)}"
    return s + (f", ended {fmt_t(ts)}" if done else "")


def pairhits(text, tag):
    """All '[tag] ... fg done X/Y pairs P Ts' tuples (X, Y, P, T)."""
    return [(int(x.replace(",", "")), int(y.replace(",", "")),
             float(p.replace(",", "")), int(t)) for x, y, p, t in
            re.findall(rf"\[{tag}\].*?fg done ([\d,]+)/([\d,]+)"
                       rf"\s+pairs ([\d,]+)\s+(\d+)s", text)]


def pairstage(text, w_data=45, w_ran=50):
    """Progress through a 02-style run: data pass then randoms pass.

    The naive per-object ETA is badly wrong on z-ordered passes: the
    opening low-z slices have huge search apertures (13 deg at z=0.075
    for 45 Mpc), so extrapolating the opening rate over-estimates by
    ~an order of magnitude (the lcr '264 h' of 2026-11-09). For the
    randoms pass the total pair count is predictable from the finished
    data pass (x ran/fg), and the PAIRS rate is nearly constant in
    time, so use pairs when possible; data-pass ETAs are upper bounds.
    """
    if "saved ->" in text or "contours ->" in text:
        return 100.0, "done"
    r = frac(text, "random")
    if r and r[0] > 0:
        eta = (r[1] - r[0]) * r[2] / r[0]        # per-object fallback
        rh, dh = pairhits(text, "random"), pairhits(text, "data")
        if rh and dh and dh[-1][0] >= dh[-1][1]:   # data pass complete
            pred = dh[-1][2] * r[1] / max(dh[-1][1], 1)  # total ran pairs
            pnow, tnow = rh[-1][2], rh[-1][3]
            old = next((h for h in reversed(rh) if tnow - h[3] >= 1800),
                       rh[0])
            dt = tnow - old[3]
            rate = (pnow - old[2]) / dt if dt > 0 else 0.0
            if rate > 0 and pred > pnow:
                eta = (pred - pnow) / rate
        return (w_data + w_ran * r[0] / max(r[1], 1),
                f"randoms pass, ETA ~{fmt_eta(eta)}")
    d = frac(text, "data")
    if d and d[0] > 0:
        eta = (d[1] - d[0]) * d[2] / d[0]
        return (w_data * d[0] / max(d[1], 1),
                f"data pass, ETA ≲{fmt_eta(eta)} + randoms")
    return (2.0, "loading catalogs") if text else (0.0, "queued")


def state_of(text, pct):
    if re.search(r"Traceback|oom_kill|OOM Killed", text):
        return "fail"
    if pct >= 100:
        return "done"
    return "run" if text else "wait"


gauges = [dict(name="DR11 sweeps (1.9 TB)", state="done", pct=100,
               unit="2126/2126 files", ts=0,
               detail="south+north complete, sha256 verified 04:10 Sat")]

t60 = read(f"{LOGS}/igmdust_60.log")
m = re.findall(r"\[(\d+)/2126\] matched so far ([\d,]+)\s+(\d+)s", t60)
if "saved ->" in t60:
    p60, d60 = 100.0, "qso_merged_dr11.fits written"
    mm = re.search(r"matched ([\d,]+) / ([\d,]+) \(([\d.]+%)\)", t60)
    if mm:
        d60 = f"{mm.group(1)} of {mm.group(2)} QSOs matched ({mm.group(3)})"
elif m:
    n, t = int(m[-1][0]), int(m[-1][2])
    p60 = 100.0 * n / 2126
    d60 = (f"sweep {n}/2126; {m[-1][1]} QSOs matched; "
           f"ETA ~{fmt_eta((2126 - n) * t / n)}")
else:
    p60, d60 = (2.0, "starting") if t60 else (0.0, "queued")
ti = runtimes(int(m[-1][2]) if m else 0, mtime(f"{LOGS}/igmdust_60.log"),
              p60 >= 100)
gauges.append(dict(name="Legacy: QSO × DR11 photometry match",
                   state=state_of(t60, p60), pct=p60,
                   unit="of 2126 sweeps",
                   detail=d60 + (f"; {ti}" if ti else ""),
                   ts=mtime(f"{LOGS}/igmdust_60.log")))

def pair_elapsed(text):
    """Seconds of counting so far: last data T + last random T."""
    return sum(h[-1][3] for h in (pairhits(text, "data"),
                                  pairhits(text, "random")) if h)


t61 = read(f"{LOGS}/igmdust_61.log")
p61, s61 = pairstage(t61)
ti = runtimes(pair_elapsed(t61), mtime(f"{LOGS}/igmdust_61.log"),
              p61 >= 100)
gauges.append(dict(name="Legacy: quasars (Ménard recipe) × BGS",
                   state=state_of(t61, p61), pct=p61, unit="pair counting",
                   detail=s61 + ("" if t61 else " (after photometry match)")
                   + (f"; {ti}" if ti else ""),
                   ts=mtime(f"{LOGS}/igmdust_61.log")))

import glob as _g
for pat, label, det in (
        ("igmdust_rmg_*.log", "DES: redMaGiC × crayons (wide)",
         "685k highdens fg × 0.65–0.90 crayons, 10 kpc–45 Mpc"),
        ("igmdust_rmw_*.log", "DES: redMaPPer × crayons (wide)",
         "19k λ>20 clusters, 10 kpc–45 Mpc"),
        ("igmdust_rmc_*.log", "DES: redMaGiC × crayons (CSFD)",
         "iteration 2: Milky Way re-correction, mid-scale test"),
        ("igmdust_mqd_*.log", "DES: quasars (Ménard recipe)",
         "redMaGiC fg × merged QSOs, DECam griz"),
        ("igmdust_lmq_*.log", "Legacy: quasars × LRG fg",
         "merged QSOs, joint dust+μ fit"),
        ("igmdust_snl_*.log", "Legacy: star null test",
         "DR11 PSF stars × BGS — systematics template"),
        ("igmdust_lcr_*.log", "Legacy: LRG crayons × BGS",
         "spec-z locus (f_abs = 1), 10 kpc–45 Mpc"),
        ("igmdust_bqw_*.log", "Legacy: BGS × QSOs wide (lw recount)",
         "45 Mpc, CSFD, locus residuals + ivar weights, 6 shards"),
        ("igmdust_lqw_*.log", "Legacy: LRG × QSOs wide (lw recount)",
         "45 Mpc, CSFD, locus residuals + ivar weights"),
        ("igmdust_pgw_*.log", "Legacy: phot. galaxies × QSOs (replica, lw)",
         "r<21 sweep galaxies, 45 Mpc, CSFD, locus+ivar wt, 4 shards"),
        ("igmdust_hcr_*.log", "HSC: BGS × quiescent crayons",
         "DR4 Wide i<24, photo-z 0.7–1.7, 45 Mpc, CSFD, 2 shards"),
        ("igmdust_hsn_*.log", "HSC: star null test",
         "i<22.5 point sources × BGS — systematics template"),
        ("igmdust_hlr_*.log", "HSC: LRG × quiescent crayons",
         "higher-z fg slice (0.4–1.1), i<24, 45 Mpc, CSFD"),
        ("igmdust_hls_*.log", "HSC: LRG star null",
         "i<22.5 point sources × LRG — systematics template")):
    # newest log drives the card: a dead shard's lexically-last log must
    # not mask the shard actually running (2026-11-11 night: s3's stale
    # MemoryError hid s1's rerun all night)
    files = sorted(_g.glob(f"{LOGS}/{pat}"), key=os.path.getmtime)
    txt = read(files[-1]) if files else ""
    p, s = pairstage(txt, w_data=10, w_ran=85)
    st = state_of(txt, p)
    # per-task logs never see the kill — sniff the SLURM job log too
    # (the 2026-11-08 OOMs sat at "loading catalogs" for hours)
    jpat = {"igmdust_hcr_*.log": "igmdust_hsc_[0-9]*.log",
            "igmdust_hsn_*.log": "igmdust_hsc_[0-9]*.log",
            "igmdust_hlr_*.log": "igmdust_hsclrg_[0-9]*.log",
            "igmdust_hls_*.log": "igmdust_hsclrg_[0-9]*.log"}.get(pat)
    if st != "done" and jpat:
        jf = sorted(_g.glob(f"{LOGS}/{jpat}"), key=os.path.getmtime)
        if jf and mtime(jf[-1]) >= (mtime(files[-1]) if files else 0) \
                and re.search(r"oom_kill|Killed", read(jf[-1])):
            st, s = "fail", "OOM-killed — needs resubmit"
    # sharded lw recounts: saved npz shards are ground truth for progress
    # (each shard log only knows itself)
    PCDIR = "/home/suzuki/data/igmdust/paircount"
    shn = {"igmdust_bqw_*.log":
           (f"{PCDIR}/paircount_bgs_mqcsfdlw_s6.*_qso.npz", 6),
           "igmdust_pgw_*.log":
           (f"{PCDIR}/paircount_photgal_mqcsfdlw_s4.*_qso.npz", 4)}.get(pat)
    if shn:
        nd = len(_g.glob(shn[0]))
        if nd >= shn[1]:
            st, p, s = "done", 100.0, f"all {shn[1]} shards saved — merge next"
        elif nd:
            p = (100.0 * nd + min(p, 99.0)) / shn[1]
            s = f"{nd}/{shn[1]} shards saved; current: {s}"
    ti = runtimes(pair_elapsed(txt), mtime(files[-1]) if files else 0,
                  st == "done")
    gauges.append(dict(name=label, state=st, pct=p,
                       unit="pair counting",
                       detail=f"{det} — {s}" + (f"; {ti}" if ti else ""),
                       ts=mtime(files[-1]) if files else 0))

# sweep-match style gauges: "[i/N] matched so far X  Ts" progress lines
for pat, label, what, obj in (
        ("igmdust_lmt_*.log", "Legacy: LRG × DR11 sweep match",
         "lrg_dr11.fits written", "LRGs"),
        ("igmdust_pgs_*.log", "Legacy: sweep scan — photometric foregrounds",
         "photgal_dr11.fits written", "r<21 galaxies"),
        ("igmdust_pgr_*.log", "Legacy: DR11 imaging randoms cut",
         "photran_dr11.fits written", "randoms")):
    files = sorted(_g.glob(f"{LOGS}/{pat}"), key=os.path.getmtime)
    tl = read(files[-1]) if files else ""
    ml = re.findall(r"\[(\d+)/(\d+)\] matched so far ([\d,]+)\s+(\d+)s", tl)
    if "saved ->" in tl:
        pl, dl = 100.0, what
    elif ml:
        n, tot, t = int(ml[-1][0]), int(ml[-1][1]), int(ml[-1][3])
        pl = 100.0 * n / tot
        dl = (f"file {n}/{tot}; {ml[-1][2]} {obj}; "
              f"ETA ~{fmt_eta((tot - n) * t / max(n, 1))}")
    else:
        pl, dl = (2.0, "starting") if tl else (0.0, "queued")
    ti = runtimes(int(ml[-1][3]) if ml else 0,
                  mtime(files[-1]) if files else 0, pl >= 100)
    gauges.append(dict(name=label, state=state_of(tl, pl), pct=pl,
                       unit="of files", detail=dl + (f"; {ti}" if ti else ""),
                       ts=mtime(files[-1]) if files else 0))

# HSC catalog downloads (3 CAS jobs; .part suffix = still streaming)
hsc_dir = "/home/suzuki/data/HSC/catalog"
jobs_hsc = ["669956", "669957", "669937"]
try:
    names = os.listdir(hsc_dir)
    done = [j for j in jobs_hsc if f"hsc_dr4_igmdust_{j}" in names]
    gb = sum(os.path.getsize(os.path.join(hsc_dir, n))
             for n in names if n.startswith("hsc_dr4")) / 1e9
    ph = 100.0 * len(done) / len(jobs_hsc)
    sh = "done" if len(done) == 3 else "run"
    gauges.append(dict(name="HSC DR4 catalog download",
                       state=sh if sh == "run" else "done", pct=ph,
                       unit=f"{len(done)}/3 files",
                       detail=f"{gb:.1f} GB on disk (CAS jobs "
                              "669956/669957/669937)",
                       ts=max((mtime(os.path.join(hsc_dir, n))
                               for n in names if n.startswith("hsc_dr4")),
                              default=0)))
except OSError:
    pass

# McCleary reproduction verdict (07_assess, SLURM job or astro4 nohup run)
files = sorted(_g.glob(f"{LOGS}/igmdust_ass_*.log"), key=os.path.getmtime)
ta = read(files[-1]) if files else ""
verdict = re.search(r"MCCLEARY_REPRODUCTION_(PASS|FAIL)", ta)
zones = re.findall(r"zone\s+(\S+): (mean pull [^\n]+)", ta)
if verdict:
    vres = ["McCleary 2026 reproduction — CSFD iteration VERDICT",
            verdict.group(1) + "; " +
            "; ".join(f"{z}: {d}" for z, d in zones)]
elif ta:
    vres = ["McCleary 2026 reproduction — CSFD iteration",
            "assessment running..."]
else:
    vres = ["McCleary 2026 reproduction — CSFD iteration",
            "waiting for CSFD pair counts (iteration 2)"]

results = [
    vres,
    ["BGS × QSO colors (joint dust+μ fit)",
     "E(B−V) = +5.7±2.3×10⁻³ at 27 kpc; "
     "pure dust, δμ ≈ 0"],
    ["BGS × DES crayons",
     "null ≥ 82 kpc (−0.04±0.56×10⁻³); "
     "photo-z absorption f = 0.806 calibrated with DESI spec-z"],
    ["LRG × QSO colors",
     "foreground-light contamination inside 500 kpc "
     "(neither dust nor magnification)"],
    ["redMaPPer λ>20 × crayons",
     "null 150–500 kpc; 787 kpc bin = R_λ selection artifact; "
     "A_V ≲ 10⁻⁴ beyond 10 Mpc"],
    ["QSO merge", "2,836,852 QSOs (DESI + SDSS DR16Q + Quaia, 1\" dedup)"],
    ["DES crayon colors", "846,578/850,028 matched to Y3 Gold (99.6%)"],
]

# order: failures first (need attention), then running, then queued,
# then done — done sorted newest-finished first, oldest at the bottom
_rank = {"fail": 0, "run": 1, "wait": 2, "done": 3}
gauges = [g for _, g in sorted(
    enumerate(gauges),
    key=lambda ig: (_rank.get(ig[1]["state"], 2),
                    -ig[1].get("ts", 0) if ig[1]["state"] == "done"
                    else ig[0]))]
for g in gauges:
    g.pop("ts", None)

print(json.dumps(dict(updated=time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                      gauges=gauges, results=results)))
