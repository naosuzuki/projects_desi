#!/bin/bash
# 00c_download_qso_catalogs.sh — igmdust: external QSO catalogs for the
# merged background sample (05_qso_merge.py):
#   SDSS DR16Q v4 superset (~1.7 GB)  — spectroscopic QSOs, mostly north
#   Gaia-unWISE Quaia G<20.5 (~1 GB)  — all-sky spectrophotometric QSOs
#
# WHERE TO RUN: astro4 (direct internet).
#   nohup bash 00c_download_qso_catalogs.sh > ~/qso_ext_dl.log 2>&1 &
# Idempotent: wget -c resumes.
set -u
DEST=/data/astrofs2_1/suzuki/data/DESI/catalog/qso_ext
mkdir -p "$DEST"
cd "$DEST" || exit 1

wget -c -nv --tries=3 --retry-connrefused \
  https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits

wget -c -nv --tries=3 --retry-connrefused \
  -O quaia_G20.5.fits \
  "https://zenodo.org/records/8060755/files/quaia_G20.5.fits?download=1"

ls -ln "$DEST"
echo "QSO_EXT_DONE $(date)"
