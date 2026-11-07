#!/usr/bin/env python3
"""merge_status.py — merge base pipeline JSON with per-host stat lines.

Usage: merge_status.py base.json hostname statfile [hostname statfile ...]
Each statfile holds the output of:
    cat /proc/loadavg; nproc; free -b | awk 'NR==2{print $2,$7}'; \
    squeue -h -o %C -t R 2>/dev/null | awk '{s+=$1} END {print "ALLOC", s+0}'
(the squeue line is absent on hosts without SLURM).
"""
import json
import sys

base = json.load(open(sys.argv[1]))
hosts = []
args = sys.argv[2:]
for name, path in zip(args[::2], args[1::2]):
    try:
        lines = open(path).read().split()
        load1 = float(lines[0])
        # tokens: loadavg(5) nproc total avail [ALLOC n]
        nproc = int(lines[5])
        total, avail = int(lines[6]), int(lines[7])
        h = dict(name=name, cores=nproc, load=load1,
                 mem_total=total, mem_used=total - avail)
        if "ALLOC" in lines:
            h["alloc"] = int(lines[lines.index("ALLOC") + 1])
            h["alloc_total"] = nproc
        hosts.append(h)
    except (OSError, ValueError, IndexError):
        pass
base["hosts"] = hosts
json.dump(base, open(sys.argv[1].replace("base", "status"), "w"))
print(f"merged {len(hosts)} hosts")
