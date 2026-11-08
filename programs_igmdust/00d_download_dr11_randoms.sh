#!/bin/bash
# 00d_download_dr11_randoms.sh — igmdust: DR11 imaging randoms (2 files,
# ~27 GB each; full footprint, same MASKBITS/NOBS columns as the sweeps).
# Needed as the foreground randoms of the Menard-replica photometric run
# (the DESI LSS randoms carry fiber-assignment holes the photometric
# galaxy sample does not have).
#
# WHERE TO RUN: astro4 (direct internet, mounts the volume):
#   ( nohup bash 00d_download_dr11_randoms.sh > ~/igmdust_ranload.log 2>&1 & )
# Writes a .downloaded marker that 13_slurm_menard.sbatch (randoms job)
# waits for.
set -u
PORTAL=https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr11/randoms
DEST=/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11/randoms
PY=/home/suzuki/miniconda3/bin/python3
mkdir -p "$DEST"
cd "$DEST" || exit 1

curl -sO --retry 5 "$PORTAL/legacysurvey_dr11_randoms.sha256sum"
for f in randoms-1-0.fits randoms-1-1.fits; do
  echo "== $f $(date)"
  curl -O --retry 5 -C - "$PORTAL/$f" || exit 1
  grep " $f\$" legacysurvey_dr11_randoms.sha256sum | sha256sum -c - || exit 1
done
echo "== cutting to photran_dr11.fits $(date)"
cd /home/suzuki/github/projects_desi/programs_igmdust || exit 1
$PY 11_sweepgal.py --mode randoms --randoms "$DEST/randoms-1-*.fits" \
  > ~/slurmlogs/igmdust_pgr_astro4.log 2>&1 || exit 1
touch "$DEST/.downloaded"
echo "== done $(date)"
