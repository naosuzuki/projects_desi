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


def frac(text, tag):
    """Last 'fg done X/Y' after a [tag] marker; returns (x, y) or None."""
    hits = re.findall(rf"\[{tag}\].*?fg done ([\d,]+)/([\d,]+)", text)
    if not hits:
        return None
    x, y = (int(v.replace(",", "")) for v in hits[-1])
    return x, y


def pairstage(text, w_data=45, w_ran=50):
    """Progress through a 02-style run: data pass then randoms pass."""
    if "saved ->" in text or "contours ->" in text:
        return 100.0, "done"
    r = frac(text, "random")
    if r:
        return w_data + w_ran * r[0] / max(r[1], 1), "randoms pass"
    d = frac(text, "data")
    if d:
        return w_data * d[0] / max(d[1], 1), "data pass"
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
m = re.findall(r"\[(\d+)/2126\] matched so far ([\d,]+)", t60)
if "saved ->" in t60:
    p60, d60 = 100.0, "qso_merged_dr11.fits written"
    mm = re.search(r"matched ([\d,]+) / ([\d,]+) \(([\d.]+%)\)", t60)
    if mm:
        d60 = f"{mm.group(1)} of {mm.group(2)} QSOs matched ({mm.group(3)})"
elif m:
    p60 = 100.0 * int(m[-1][0]) / 2126
    d60 = f"sweep {m[-1][0]}/2126; {m[-1][1]} QSOs matched so far"
else:
    p60, d60 = (2.0, "starting") if t60 else (0.0, "queued")
gauges.append(dict(name="DR11 photometry match (job 60)",
                   state=state_of(t60, p60), pct=p60,
                   unit="of 2126 sweeps", detail=d60))

t61 = read(f"{LOGS}/igmdust_61.log")
p61, s61 = pairstage(t61)
gauges.append(dict(name="BGS × merged QSOs (job 61)",
                   state=state_of(t61, p61), pct=p61, unit="pair counting",
                   detail=s61 + ("" if t61 else " (after job 60)")))

import glob as _g
for pat, label, det in (
        ("igmdust_rmg_*.log", "redMaGiC × crayons (wide)",
         "685k highdens fg × 0.65–0.90 crayons, 10 kpc–45 Mpc"),
        ("igmdust_rmw_*.log", "redMaPPer × crayons (wide)",
         "19k λ>20 clusters, 10 kpc–45 Mpc")):
    files = sorted(_g.glob(f"{LOGS}/{pat}"))
    txt = read(files[-1]) if files else ""
    p, s = pairstage(txt, w_data=10, w_ran=85)
    gauges.append(dict(name=label, state=state_of(txt, p), pct=p,
                       unit="pair counting", detail=f"{det} — {s}"))

results = [
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
