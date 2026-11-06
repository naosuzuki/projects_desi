#!/bin/bash
# 00_download_dr1_catalogs.sh — igmdust phase 0: fetch the DESI DR1 pilot catalogs.
#
# Downloads the DR1 QSO catalog and the BGS_ANY/LRG LSS clustering catalogs
# (data + 2 randoms per cap + n(z)) into $DEST. Idempotent: wget -c resumes.
#
# WHERE TO RUN: astro4 (gateway; mounts /data/astrofs2_1, direct internet, no SLURM).
#   Do NOT run on astrofs4: its srun login wrapper turns jobs into login blockers.
#   Usage on astro4:  nohup bash 00_download_dr1_catalogs.sh > dl.log 2>&1 &
#
# Completed 2026-11-06 (8.8 GiB, 13 FITS files; QSO catalog sha256 verified).

set -u
DEST=/data/astrofs2_1/suzuki/data/DESI/catalog/dr1
LSS=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5
QSO=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron

mkdir -p "$DEST"
cd "$DEST" || exit 1

FILES=(
  BGS_ANY_NGC_clustering.dat.fits
  BGS_ANY_SGC_clustering.dat.fits
  LRG_NGC_clustering.dat.fits
  LRG_SGC_clustering.dat.fits
  BGS_ANY_NGC_0_clustering.ran.fits
  BGS_ANY_NGC_1_clustering.ran.fits
  BGS_ANY_SGC_0_clustering.ran.fits
  BGS_ANY_SGC_1_clustering.ran.fits
  LRG_NGC_0_clustering.ran.fits
  LRG_NGC_1_clustering.ran.fits
  LRG_SGC_0_clustering.ran.fits
  LRG_SGC_1_clustering.ran.fits
  BGS_ANY_NGC_nz.txt
  BGS_ANY_SGC_nz.txt
  LRG_NGC_nz.txt
  LRG_SGC_nz.txt
)

for f in "${FILES[@]}"; do wget -c -nv "$LSS/$f"; done
wget -c -nv "$QSO/QSO_cat_iron_cumulative_v0.fits"
wget -c -nv "$QSO/dr1_survey_catalogs_dr1_QSO_iron.sha256sum"

grep QSO_cat dr1_survey_catalogs_dr1_QSO_iron.sha256sum | sha256sum -c -
echo "DOWNLOAD_DONE $(date)"
