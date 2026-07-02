#!/usr/bin/env bash
#
# 00_download_zall_tilecumulative_globus.sh
# -----------------------------------------
# Download the DESI DR2 ("loa") tile-cumulative redshift (specz) catalog from
# NERSC using GLOBUS (fault-tolerant, parallel, resumable, checksum-verified).
#
#   zall-tilecumulative-loa.fits
#
# Source : NERSC DTN endpoint  ->  /global/cfs/cdirs/desi/public/dr2/...
# Dest   : your Globus Connect Personal collection on this Mac
#
# ---------------------------------------------------------------------------
# ONE-TIME SETUP (do these once, by hand):
#   1. globus-cli is installed at ~/.local/bin/globus  (already done).
#   2. Log in (opens a browser):            ~/.local/bin/globus login
#   3. In the Globus Connect Personal app:  Preferences -> Access ->
#        add  /Volumes/exdisk2/desiredux/DR2  (or /Volumes/exdisk2) as a
#        writable folder, then make sure Globus Connect Personal is RUNNING.
#   4. Accept the NERSC DTN endpoint's terms once via the Globus web app if
#      prompted (data access consent).
# ---------------------------------------------------------------------------
#
# Usage:
#   ./00_download_zall_tilecumulative_globus.sh           # submit + wait
#   WAIT=0 ./00_download_zall_tilecumulative_globus.sh     # submit, don't block
#
set -euo pipefail

# --- Tools -------------------------------------------------------------------
GLOBUS="${GLOBUS:-$HOME/.local/bin/globus}"
command -v "${GLOBUS}" >/dev/null 2>&1 || GLOBUS="globus"   # fall back to PATH

# --- Configuration -----------------------------------------------------------
RELEASE="dr2"
RELEASE_UC="$(printf '%s' "$RELEASE" | tr '[:lower:]' '[:upper:]')"
REDUX="loa"
FILENAME="zall-tilecumulative-${REDUX}.fits"
REMOTE_REL="${RELEASE}/spectro/redux/${REDUX}/zcatalog/v1/${FILENAME}"

SRC_PATH="/global/cfs/cdirs/desi/public/${REMOTE_REL}"

# NERSC DTN managed collection (long-standing, stable UUID).
# Verify any time with:  globus endpoint search 'NERSC DTN'
SRC_EP="${SRC_EP:-9d6d994a-6d04-11e5-ba46-22000b92c6ec}"

# This Mac's Globus Connect Personal collection (auto-read from Globus Connect Personal config).
DEST_EP="${DEST_EP:-$(cat "${HOME}/.globusonline/lta/client-id.txt" 2>/dev/null)}"

DEST_DIR="${DEST_DIR:-/Volumes/exdisk2/desiredux/DR2}"
DEST_DIR="${DEST_DIR%/}"
DEST_PATH="${DEST_DIR}/${FILENAME}"

WAIT="${WAIT:-1}"   # 1 = block until the transfer finishes; 0 = submit and exit

echo "=================================================================="
echo " DESI ${RELEASE_UC} (${REDUX}) specz catalog download  (Globus)"
echo "   file : ${FILENAME}"
echo "   src  : NERSC DTN  ${SRC_EP}"
echo "          ${SRC_PATH}"
echo "   dest : local Globus Connect Personal  ${DEST_EP}"
echo "          ${DEST_PATH}"
echo "=================================================================="

# --- Preflight ---------------------------------------------------------------
if [[ -z "${DEST_EP}" ]]; then
    echo "ERROR: could not determine local Globus collection UUID."
    echo "       Is Globus Connect Personal installed? Set DEST_EP=... to override."
    exit 1
fi

if ! "${GLOBUS}" whoami >/dev/null 2>&1; then
    echo "ERROR: not logged in to Globus. Run:"
    echo "         ${GLOBUS} login"
    exit 1
fi

if ! pgrep -fi "Globus Connect Personal" >/dev/null 2>&1; then
    echo "WARNING: Globus Connect Personal does not appear to be running."
    echo "         Start the Globus Connect Personal app (the local endpoint must be online to receive data)."
    echo
fi

# --- Submit the transfer -----------------------------------------------------
#   --sync-level checksum : skip/repair only what's needed; safe to re-run
#   --verify-checksum     : end-to-end integrity check (default on, explicit here)
#   --fail-on-quota-errors: surface dest-side problems instead of silently retrying
LABEL="DESI ${RELEASE} ${REDUX} zall-tilecumulative"
echo "Submitting transfer..."
TASK_ID="$(
    "${GLOBUS}" transfer \
        --label "${LABEL}" \
        --sync-level checksum \
        --verify-checksum \
        --fail-on-quota-errors \
        --jmespath 'task_id' --format unix \
        "${SRC_EP}:${SRC_PATH}" \
        "${DEST_EP}:${DEST_PATH}"
)"
echo "Task submitted: ${TASK_ID}"
echo "Monitor:  ${GLOBUS} task show ${TASK_ID}"
echo "Web UI :  https://app.globus.org/activity/${TASK_ID}"

# --- Optionally block until done ---------------------------------------------
if [[ "${WAIT}" == "1" ]]; then
    echo "Waiting for completion (Ctrl-C just stops waiting; the transfer keeps running)..."
    if "${GLOBUS}" task wait "${TASK_ID}" --polling-interval 15; then
        echo "------------------------------------------------------------------"
        echo "Transfer SUCCEEDED."
        echo "File should be at: ${DEST_PATH}"
        echo "(If the dest disk is local, verify: head -c 6 \"${DEST_PATH}\"  -> SIMPLE)"
    else
        echo "Transfer did not succeed; inspect with:"
        echo "  ${GLOBUS} task show ${TASK_ID}"
        echo "  ${GLOBUS} task event-list ${TASK_ID}"
        exit 1
    fi
fi
