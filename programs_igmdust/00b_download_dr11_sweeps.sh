#!/bin/bash
# 00b_download_dr11_sweeps.sh — igmdust: fetch the Legacy Surveys DR11 MAIN
# sweep catalogs (the program photometry; approved by Nao 2026-11-06).
#
# Scope: dr11/{south,north}/sweep/11.0/ main sweeps only —
#   south: 1,600 files ~3.0 TB; north: 526 files ~0.1 TB; total ~3.1 TB.
#   The row-matched -lightcurves/-extra/-photo-z trees are NOT downloaded.
#
# WHERE TO RUN: astro4, and (since the 2026-11-06 admin fix) astrofs4 via
#   a small SLURM job. Two hosts run DISJOINT shards of the same list:
#     astro4  : nohup bash 00b_download_dr11_sweeps.sh 1 2 > ~/dr11_dl.log &
#     astrofs4: sbatch (2 cores) -> bash 00b_download_dr11_sweeps.sh 2 2
#   Args: SHARD NSHARD (default 1 1 = everything). Sharding is by line
#   number of the sorted url list, so shards never overlap; wget -c skips
#   files that are already complete. NPAR (env) sets streams per host.
#
# Idempotent: file lists are rebuilt from the portal index and aria2c -c
# resumes partial files; three passes catch transient failures.
# Parallel: aria2c -j6 (6 concurrent files, 1 connection each) — measured
# 9 MiB/s single-stream 2026-11-06; parallel streams cut ~4 days to ~1.

set -u
DEST=/data/astrofs2_1/suzuki/data/DESI/catalog/ls/dr11
PORTAL=https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr11
SHARD=${1:-1}
NSHARD=${2:-1}
NPAR=${NPAR:-8}
echo "== shard $SHARD of $NSHARD on $(hostname), $NPAR streams"

for PART in south north; do
  mkdir -p "$DEST/$PART"
  cd "$DEST/$PART" || exit 1
  BASE="$PORTAL/$PART/sweep/11.0"
  # build the file list from the index: main sweeps + checksum files
  curl -s --retry 5 "$BASE/" \
    | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//' \
    | grep -E '^sweep-[0-9mp-]+\.fits$|sha256sum' | sort -u \
    | sed "s|^|$BASE/|" > urls_$PART.txt
  SL=urls_${PART}_shard${SHARD}of${NSHARD}.txt       # unique per shard:
  awk -v s="$SHARD" -v n="$NSHARD" 'NR % n == s - 1' urls_$PART.txt > "$SL"
  n=$(wc -l < "$SL")                                 # both hosts share $DEST
  echo "== $PART: $n files in shard $SHARD/$NSHARD -> $DEST/$PART"
  # parallel wget streams (native tools only: conda aria2c lacks libstdc++
  # on astro4). Light on the host: ~10 MB per wget, nice -n 19.
  for pass in 1 2 3; do
    echo "== $PART pass $pass $(date)"
    xargs -P "$NPAR" -n 1 nice -n 19 wget -c -nv --retry-connrefused \
      --tries=3 < "$SL" && break
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
