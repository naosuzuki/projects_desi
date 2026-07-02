#!/usr/bin/env bash
#
# 00_download_zall_tilecumulative.sh
# ----------------------------------
# Download the DESI DR2 ("loa") tile-cumulative redshift (specz) catalog
# directly from NERSC via rsync/scp over SSH (MFA).
#
#   zall-tilecumulative-loa.fits
#
# NERSC source path (community filesystem):
#   /global/cfs/cdirs/desi/public/dr2/spectro/redux/loa/zcatalog/v1/zall-tilecumulative-loa.fits
#
# Usage:
#   NERSC_USER=yourname ./00_download_zall_tilecumulative.sh
#   NERSC_USER=yourname DEST_DIR=/path/to/dir ./00_download_zall_tilecumulative.sh
#
# Authentication (NERSC MFA):
#   - Plain rsync (default below) will prompt for: <password><OTP>  (your
#     password immediately followed by the 6-digit Google Authenticator code,
#     no space, on the SAME line).
#   - For repeated transfers without re-typing, get a 24-hour key first:
#         sshproxy.sh -u $NERSC_USER          # one MFA prompt, key good 24h
#     then this script's ssh/rsync run key-only. Get sshproxy from:
#         https://docs.nersc.gov/connect/mfa/#sshproxy
#
# The transfer is resumable: re-running continues a partial file (rsync -P).
#
set -euo pipefail

# --- Configuration ----------------------------------------------------------
RELEASE="dr2"
RELEASE_UC="$(printf '%s' "$RELEASE" | tr '[:lower:]' '[:upper:]')"
REDUX="loa"
FILENAME="zall-tilecumulative-${REDUX}.fits"

REMOTE_REL="${RELEASE}/spectro/redux/${REDUX}/zcatalog/v1/${FILENAME}"
NERSC_PATH="/global/cfs/cdirs/desi/public/${REMOTE_REL}"

NERSC_USER="${NERSC_USER:-${USER}}"            # set to your NERSC login if different
NERSC_HOST="${NERSC_HOST:-dtn01.nersc.gov}"    # NERSC data transfer node (best for bulk copies)

# Destination: external disk on this machine (override with DEST_DIR=...).
DEST_DIR="${DEST_DIR:-/Volumes/exdisk2/desiredux/DR2}"
DEST_DIR="${DEST_DIR%/}"
DEST="${DEST_DIR}/${FILENAME}"

mkdir -p "${DEST_DIR}"

echo "=================================================================="
echo " DESI ${RELEASE_UC} (${REDUX}) specz catalog download  (NERSC rsync)"
echo "   file   : ${FILENAME}"
echo "   source : ${NERSC_USER}@${NERSC_HOST}:${NERSC_PATH}"
echo "   dest   : ${DEST}"
echo "=================================================================="

if [[ "${NERSC_USER}" == "${USER}" ]]; then
    echo "NOTE: NERSC_USER defaulted to local user '${USER}'."
    echo "      If your NERSC login differs, re-run with NERSC_USER=<login>."
    echo
fi

# --- Download (rsync over SSH, resumable, with progress) --------------------
#   -a  archive (preserve perms/times)
#   -v  verbose
#   -P  --partial --progress  (resume + live progress bar)
#   --append-verify  resume a partial file then checksum-verify the whole thing
rsync -avP --append-verify \
      "${NERSC_USER}@${NERSC_HOST}:${NERSC_PATH}" \
      "${DEST_DIR}/"

# --- Verify ------------------------------------------------------------------
echo "------------------------------------------------------------------"
if [[ -s "${DEST}" ]]; then
    SIZE_H="$(du -h "${DEST}" | cut -f1)"
    echo "Done: ${DEST} (${SIZE_H})"
    # Sanity check: a real FITS file starts with the ASCII string 'SIMPLE'.
    if [[ "$(head -c 6 "${DEST}")" == "SIMPLE" ]]; then
        echo "FITS header looks valid (starts with SIMPLE)."
    else
        echo "WARNING: file does not start with 'SIMPLE' — inspect it:"
        echo "           head -c 300 \"${DEST}\""
        exit 1
    fi
else
    echo "ERROR: transfer produced an empty file."
    exit 1
fi

# =============================================================================
# Alternative one-liners (if you prefer not to run this script):
#
#   scp   "${NERSC_USER}@${NERSC_HOST}:${NERSC_PATH}"  "${DEST}"
#   rsync -avP "${NERSC_USER}@${NERSC_HOST}:${NERSC_PATH}" "${DEST_DIR}/"
#
# Globus is the most robust route for very large / many-file transfers:
#   endpoint  "NERSC DTN"  ->  your local Globus Connect Personal endpoint.
# =============================================================================
