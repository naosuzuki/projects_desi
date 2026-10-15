# DESI DR2 spectra collection — runbook (state 2026-10-15, before a planned shutdown)

## Where the data are
DESI DR2 (production **loa**) HEALPix coadds, `main` survey, `dark` + `bright`: **34,508 files, ~24.03 TB**.
Source: NERSC `/global/cfs/cdirs/desi/public/dr2/spectro/redux/loa/healpix/main/{bright,dark}`.

| disk | path | files |
|---|---|---|
| exdisk2 (24.0 TB, Mac OS Extended Journaled, ~30 GB free) | `/Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main/...` | 33,539 |
| exdisk1 (`/Volumes/exdisk1` -> `/Volumes/My Book`, exFAT) | `/Volumes/exdisk1/desiredux/DR2/spectro/redux/loa/healpix/main/...` | 969 (overflow) |

The overflow list is `desiredux/DR2_overflow_on_exdisk1.txt` (on both disks). All 34,508 present, none missing.
DR2 does not fit on one 24 TB disk; the old "17 TB" estimate (494 MB/file) was wrong, and the pre-Oct-9 exFAT copy
was never complete (782 files missing, ~16 % corrupt). The tile catalog `zall-tilecumulative-loa.fits` (59 GB) is on
the SSD in `~/data/exdisk2_keep_tmp/`; `zall-pix-loa.fits` is in `~/desiredux/`.

## TODO after the reboot: finish the integrity check
Verification = one Globus task per disk with `--sync-level checksum`: checksums at NERSC and here are compared; a
corrupted file is transferred again over the bad copy. Started 2026-10-15 10:37, interrupted by the shutdown
(481 files checked, all identical, 0 re-downloaded). The first tasks had Globus' default 24 h deadline, too short.

1. Mount both disks; check `ls /Volumes/exdisk2/desiredux /Volumes/exdisk1/desiredux`
   (`/Volumes/exdisk1` must be the link to `/Volumes/My Book`; recreate with `sudo ln -s "/Volumes/My Book" /Volumes/exdisk1`).
2. Start Globus Connect Personal (`open -a "Globus Connect Personal"`); `~/.local/bin/globus whoami`.
3. Rerun the verification (cancels the old verify tasks, resubmits both with a 10-day deadline, waits, reports):
   ```
   cd ~/github/projects_desi/programs_setups
   python3 -c "import os;os.fork() and os._exit(0);os.setsid();os.fork() and os._exit(0);os.execvp('zsh',['zsh','-c','exec zsh 06_verify_redownload.sh >> logs/verify.out 2>&1'])"
   ```
   ~2 days (all ~24 TB are read back). Result: `logs/verify.log` ("VERIFY DONE"; per disk: files, re-downloaded,
   identical, faults) and `logs/verify_redownloaded_{exdisk2,exdisk1}.txt` (files that were corrupted and replaced).

## Scripts (programs_setups/)
| script | role |
|---|---|
| `01_build_manifests.py` | manifests (`manifests/healpix_main_darkbright.batch` = the 34,508 files) |
| `02_transfer_healpix.sh` | the main Globus transfer to exdisk2 (task 07cf18f4…, stopped by the guard 2026-10-15 09:31) |
| `03_monitor_transfer.sh` | 10-min rate/space log (`logs/transfer_rate.log`); falls back to disk growth when the Globus API is down |
| `04_disk_guard.sh [GB]` | cancels the transfer below N GB free on exdisk2 (was run with 30) |
| `05_overflow_to_exdisk1.sh` | overflow of what did not fit -> exdisk1 (used manually 2026-10-15: task d460fe25…, 969 files, 116 GB) |
| `06_verify_redownload.sh` | checksum verification + re-download of corrupted files (rerunnable) |

## Gotchas
* Globus Connect Personal folder grant: `~/Library/Preferences/org.globusonline.Globus-Connect.plist`
  `GC_RESTRICTED_PATHS` (granted: `/Volumes/exdisk2/desiredux`, `/Volumes/My Book/desiredux`, `/Users/suzuki`).
  Batch files need quotes around paths with spaces ("My Book").
* Network: prefer Ethernet; on Wi-Fi with IPv6 Globus can fail "IPv4 to IPv6" (`sudo networksetup -setv6off "Wi-Fi"`).
* Rate falls as exdisk2 fills (inner tracks + read-back for checksum): 58 MB/s empty -> ~40 MB/s at 99 %.
* Each coadd: B/R/Z RESOLUTION HDUs are ~78 % of the bytes; FLUX+IVAR+MASK ~22 %.
* QSO catalog / NERSC-side QSO extraction (old plan steps 5-6) not started.
