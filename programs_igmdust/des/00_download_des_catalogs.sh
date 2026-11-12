#!/bin/bash
# 00_download_des_catalogs.sh — igmdust/DES: fetch the DES Y3A2 redMaGiC
# catalogs (the "standard crayon" background sample of McCleary+25) with
# their randoms and the joint mask.
#
# Files (desdr-server.ncsa.illinois.edu/despublic/y3a2_files/redmagic/):
#   y3a2_gold2.2.1_redmagic_highdens.fits            121 MB  z~0.15-0.65
#   y3a2_gold2.2.1_redmagic_highdens_randoms.fits    1.7 GB
#   y3a2_gold2.2.1_redmagic_highlum_highz.fits        56 MB  z~0.4-0.9 (crayons)
#   y3a2_gold2.2.1_redmagic_highlum_highz_randoms.fits 782 MB
#   y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits 327 MB
#   total ~3.0 GB
#
# WHERE TO RUN: astro4 (direct internet, no SLURM; never astrofs4 jobs).
#   nohup bash 00_download_des_catalogs.sh > ~/des_dl.log 2>&1 &
# DEST per Nao (2026-11-06): /home/suzuki/data/DES/catalog
#
# The DES DR2 coadd pull (des_dr2.main via Data Lab TAP: star null test +
# calibration cross-check) is a separate script: 00b (to be written).

set -u
DEST=/home/suzuki/data/DES/catalog/redmagic
BASE=https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/redmagic

mkdir -p "$DEST"
cd "$DEST" || exit 1

FILES=(
  y3a2_gold2.2.1_redmagic_highdens.fits
  y3a2_gold2.2.1_redmagic_highdens_randoms.fits
  y3a2_gold2.2.1_redmagic_highlum_highz.fits
  y3a2_gold2.2.1_redmagic_highlum_highz_randoms.fits
  y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits
)

printf '%s\n' "${FILES[@]}" | sed "s|^|$BASE/|" > urls_redmagic.txt
for pass in 1 2 3; do
  echo "== redmagic pass $pass $(date)"
  xargs -P 4 -n 1 nice -n 19 wget -c -nv --retry-connrefused --tries=3 \
    < urls_redmagic.txt && break
done
ls -la "$DEST"
echo "DES_REDMAGIC_DONE $(date)"
