#!/usr/bin/env bash
#
# 02_transfer_healpix.sh
# ----------------------
# Bulk-transfer the DESI DR2 (loa) HEALPIX coadd files (STAR+GALAXY science)
# from NERSC to exdisk2 via Globus, using the manifest built by
# 01_build_manifests.py.
#
# Scope (chosen): main dark+bright  -> manifests/healpix_main_darkbright.batch
#   34,508 files, ~17 TB.  (Swap BATCH below for healpix_coadd_files.batch to
#   pull everything = ~28 TB, which will NOT fit 22 TB.)
#
# Prereqs (one-time, already done this session):
#   - globus login ; session update nersc.gov ; NERSC data_access consent
#   - Globus Connect Personal RUNNING with /Volumes/exdisk2 granted (Access tab)
#
# SPEED NOTE: a single big file (like the catalog) is capped ~40 Mbps by the
# 81 ms NERSC RTT + one TCP stream. A many-file batch like this parallelises
# across files (Globus concurrency) and should be much faster. To push harder,
# set Globus Connect Personal -> Preferences -> Network use -> "Aggressive"
# (or Custom: more concurrency + parallel streams) BEFORE launching.
#
set -euo pipefail

GLOBUS="${GLOBUS:-$HOME/.local/bin/globus}"
SRC_EP="${SRC_EP:-9d6d994a-6d04-11e5-ba46-22000b92c6ec}"             # NERSC DTN
DEST_EP="${DEST_EP:-$(cat "$HOME/.globusonline/lta/client-id.txt")}" # this Mac (Globus Connect Personal)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH="${BATCH:-$HERE/manifests/healpix_main_darkbright.batch}"

[ -f "$BATCH" ] || { echo "ERROR: batch file not found: $BATCH"; exit 1; }
NFILES=$(wc -l < "$BATCH")

echo "=================================================================="
echo " DESI DR2 (loa) HEALPIX coadd bulk transfer  (Globus)"
echo "   batch : $BATCH"
echo "   files : $NFILES   (~$(awk -v n=$NFILES 'BEGIN{printf "%.1f", n*494e6/1e12}') TB)"
echo "   src   : NERSC DTN $SRC_EP"
echo "   dest  : local Globus Connect Personal $DEST_EP"
echo "=================================================================="

# Preflight
$GLOBUS whoami >/dev/null 2>&1 || { echo "Not logged in: $GLOBUS login"; exit 1; }
pgrep -fi "Globus Connect Personal" >/dev/null 2>&1 || \
    echo "WARNING: Globus Connect Personal not running — start it first."

# Submit one task for the whole batch.
#   --sync-level checksum : safe to re-run; only transfers what's missing/changed
#   --preserve-timestamp  : keep NERSC mtimes
#   --verify-checksum     : end-to-end integrity (default; explicit here)
#   --label               : easy to find in the web UI
TASK_ID="$(
  $GLOBUS transfer \
    --batch "$BATCH" \
    --label "DESI DR2 loa healpix main dark+bright" \
    --sync-level checksum \
    --preserve-timestamp \
    --verify-checksum \
    --jmespath 'task_id' --format unix \
    "$SRC_EP" "$DEST_EP"
)"

echo "Task submitted: $TASK_ID"
echo "Monitor : $GLOBUS task show $TASK_ID"
echo "Web UI  : https://app.globus.org/activity/$TASK_ID"
echo
echo "This runs server-side. It survives logout/reboot: when this Mac is off the"
echo "task pauses, and resumes automatically once Globus Connect Personal is back"
echo "(within the task deadline). Re-running this script is safe (checksum sync)."
