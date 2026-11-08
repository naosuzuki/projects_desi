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


def pairstage(text, w_data=45, w_ran=50):
    """Progress through a 02-style run: data pass then randoms pass."""
    if "saved ->" in text or "contours ->" in text:
        return 100.0, "done"
    r = frac(text, "random")
    if r and r[0] > 0:
        eta = (r[1] - r[0]) * r[2] / r[0]
        return (w_data + w_ran * r[0] / max(r[1], 1),
                f"randoms pass, ETA ~{fmt_eta(eta)}")
    d = frac(text, "data")
    if d and d[0] > 0:
        eta = (d[1] - d[0]) * d[2] / d[0]
        return (w_data * d[0] / max(d[1], 1),
                f"data pass, ETA ~{fmt_eta(eta)} + randoms")
    return (2.0, "loading catalogs") if text else (0.0, "queued")


def state_of(text, pct):
    if re.search(r"Traceback|oom_kill|OOM Killed", text):
        return "fail"
    if pct >= 100:
        return "done"
    return "run" if text else "wait"


gauges = [dict(name="DR11 sweeps (1.9 TB)", state="done", pct=100,
               unit="2126/2126 files",
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
gauges.append(dict(name="Legacy: QSO × DR11 photometry match",
                   state=state_of(t60, p60), pct=p60,
                   unit="of 2126 sweeps", detail=d60))

t61 = read(f"{LOGS}/igmdust_61.log")
p61, s61 = pairstage(t61)
gauges.append(dict(name="Legacy: quasars (Ménard recipe) × BGS",
                   state=state_of(t61, p61), pct=p61, unit="pair counting",
                   detail=s61 + ("" if t61 else " (after photometry match)")))

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
        ("igmdust_bqw_*.log", "Legacy: BGS × QSOs wide (CSFD)",
         "fix wave: 45 Mpc, CSFD, 6 shards + cov-ML"),
        ("igmdust_lqw_*.log", "Legacy: LRG × QSOs wide (CSFD)",
         "fix wave: 45 Mpc, CSFD + cov-ML"),
        ("igmdust_pgw_*.log", "Legacy: phot. galaxies × QSOs (Ménard replica)",
         "r<21 sweep galaxies (~17× BGS pairs), 45 Mpc, CSFD, astrofs4")):
    files = sorted(_g.glob(f"{LOGS}/{pat}"))
    txt = read(files[-1]) if files else ""
    p, s = pairstage(txt, w_data=10, w_ran=85)
    gauges.append(dict(name=label, state=state_of(txt, p), pct=p,
                       unit="pair counting", detail=f"{det} — {s}"))

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
    gauges.append(dict(name=label, state=state_of(tl, pl), pct=pl,
                       unit="of files", detail=dl))

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
                              "669956/669957/669937)"))
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

print(json.dumps(dict(updated=time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                      gauges=gauges, results=results)))
