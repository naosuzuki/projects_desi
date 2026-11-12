#!/bin/bash
# 00_submit_hsc_queries.sh — igmdust/HSC: submit the two catalog queries to
# the HSC DR4 (S23B) internal CAS and fetch the results.
#
# YOU (Nao) must run this — it needs your STARS credentials. It never
# stores your password: you are prompted by the official query tool.
#
# One-time setup:
#   1. Log in to https://hscdata.mtk.nao.ac.jp/datasearch/ and download the
#      batch tool hscReleaseQuery.py ("Command-line SQL tool" link), put it
#      in this directory.
#   2. export HSC_SSP_CAS_USER=<your STARS user>   (password is prompted)
#
# Run (on astro4, from this directory):
#   bash 00_submit_hsc_queries.sh
#
# Results land in /home/suzuki/data/HSC/catalog/ as FITS:
#   hsc_crayon_galaxies.fits  (~15-25 GB: i<25 galaxies, photoz 0.4-1.7)
#   hsc_null_stars.fits       (~1 GB: i<22.5 point sources)

set -eu
DEST=/home/suzuki/data/HSC/catalog
TOOL=./hscReleaseQuery.py
USER_OPT=${HSC_SSP_CAS_USER:?export HSC_SSP_CAS_USER=<STARS user> first}

[ -f "$TOOL" ] || { echo "Download hscReleaseQuery.py first (see header)"; exit 1; }
# the public tool hard-codes the public release names (pdr1..pdr3);
# drop the client-side choices list so the internal server validates
# the release name itself
sed "s/choices='pdr1 [^']*'\.split(), //" "$TOOL" > hscQuery_internal.py
APIURL=https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/
RELEASE=${HSC_RELEASE:-dr4}   # override with HSC_RELEASE=... if rejected
mkdir -p "$DEST"

for q in hsc_crayon_galaxies hsc_null_stars; do
  echo "== submitting $q  $(date)"
  python3 hscQuery_internal.py --user "$USER_OPT" --api-url "$APIURL" \
    --release-version "$RELEASE" --format fits --delete-job \
    sql/$q.sql > "$DEST/$q.fits"
  echo "== done: $(ls -la $DEST/$q.fits)"
done
echo "HSC_QUERIES_DONE $(date)"
