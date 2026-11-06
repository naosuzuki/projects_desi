#!/bin/bash
# 00b_download_dr11_sweeps.sh — igmdust: fetch the Legacy Surveys DR11 MAIN
# sweep catalogs (the program photometry; approved by Nao 2026-11-06).
#
# Scope: dr11/{south,north}/sweep/11.0/ main sweeps only —
#   south: 1,600 files ~3.0 TB; north: 526 files ~0.1 TB; total ~3.1 TB.
#   The row-matched -lightcurves/-extra/-photo-z trees are NOT downloaded.
#
# WHERE TO RUN: astro4 (gateway; mounts /data/astrofs2_1, direct internet,
#   no SLURM). Do NOT run on astrofs4 (srun wrapper turns jobs into login
#   blockers). Usage:
#     nohup bash 00b_download_dr11_sweeps.sh > ~/dr11_dl.log 2>&1 &
#
# Idempotent: file lists are rebuilt from the portal index and wget -c
# resumes partial files; three passes catch transient failures.

set -u
DEST=/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11
PORTAL=https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr11

for PART in south north; do
  mkdir -p "$DEST/$PART"
  cd "$DEST/$PART" || exit 1
  BASE="$PORTAL/$PART/sweep/11.0"
  # build the file list from the index: main sweeps + checksum files
  curl -s --retry 5 "$BASE/" \
    | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//' \
    | grep -E '^sweep-[0-9mp-]+\.fits$|sha256sum' | sort -u \
    | sed "s|^|$BASE/|" > urls_$PART.txt
  n=$(wc -l < urls_$PART.txt)
  echo "== $PART: $n files queued -> $DEST/$PART"
  for pass in 1 2 3; do
    echo "== $PART pass $pass $(date)"
    wget -c -nv --retry-connrefused --tries=3 -i urls_$PART.txt && break
  done
done

echo "== verifying checksums (where provided)"
for PART in south north; do
  cd "$DEST/$PART" || continue
  for s in *.sha256sum; do
    [ -f "$s" ] && sha256sum -c "$s" 2>/dev/null | grep -vc ": OK$" \
      | xargs -I{} echo "$PART/$s: {} files NOT ok (0 = all good)"
  done
done
echo "DR11_SWEEPS_DONE $(date)"
