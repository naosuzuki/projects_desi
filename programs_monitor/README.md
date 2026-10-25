# programs_monitor — pipeline monitor (all tracks, all steps)

`monitor.py` + `monitor.html`: one page at **http://127.0.0.1:8765/monitor.html** (`.claude/launch.json` → "monitor"
attaches to the running server; "monitor-start" starts one) showing, per track (common / galaxies / quasars / stars),
every step as a gauge with done/total, rate, **ETA and finish time**, its workers' CPU % and RSS, plus system CPU
(total, per core, 1-h history), memory, swap, disk read, load, disk free, the worker table and the guard events.
Supersedes `programs_templates/04_monitor.py` (same status dir `~/data/sdssdesi_status/`, same `status.json`, same guard).

Steps are declared in `STEPS` in monitor.py. Two kinds:
- `files`: count output files of a glob against a known total (pass B, pass Q SDSS, pass Q DESI).
- `json`: the step writes `~/data/sdssdesi_status/step_<name>.json`:
  `{"done": n, "total": N, "status": "running|done|idle", "note": "...", "started": <unix time>}`
  (cell tables, stacks, EM-PCA, star steps). A step with no file shows as idle; a worker process matching the step's
  `proc` string makes it "running" even without a file.
Guard (unchanged): SIGSTOP all workers if available memory < 8 GB or swap > 6 GB (resume at 16 GB); SIGSTOP the
highest shard if CPU > 92 % for 90 s while our workers are > 40 % of the machine (resume < 70 %). Run **one** monitor
with the guard; `--no-guard` for a second instance (e.g. another port).

```
cd ~/github/projects_desi/programs_monitor
nohup python3 monitor.py --port 8765 > ~/data/sdssdesi_status/monitor.log 2>&1 &
python3 monitor.py --once --no-guard      # print one status snapshot
```
Started 2026-10-24 21:03 (replacing 04_monitor.py). ETAs are from the last 30 min of output files, so they move when
the machine load changes (2026-10-24 evening: Prospector + pass B + pass Q SDSS → load 14/14, pass B ETA slipped from
Sunday 07:00 to Sunday ~23:00; the SDSS pass was reniced to 19 to favour pass B).
