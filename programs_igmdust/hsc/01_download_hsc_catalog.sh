#!/bin/bash
# 01_download_hsc_catalog.sh — igmdust/HSC: fetch the finished DR4 CAS
# catalog jobs (i<26 crayon/background sample, split by redshift because
# of its size; submitted by 00_submit_hsc_queries.sh).
#
#   job 669956  dr4  done 2026-11-07 (low-z half)
#   job 669957  dr4  done 2026-11-07 (high-z half)
#
# Auth: ~/.netrc entry for hscdata.mtk.nao.ac.jp (curl -n; no password
# in command lines or logs). WHERE TO RUN: astro4.
#   nohup bash 01_download_hsc_catalog.sh > ~/hsc_dl.log 2>&1 &
# Idempotent: curl -C - resumes partial files.
set -u
DEST=/home/suzuki/data/HSC/catalog
BASE=https://hscdata.mtk.nao.ac.jp/datasearch/catalog_jobs/download
mkdir -p "$DEST"
cd "$DEST" || exit 1

fetch () {  # args: jobid token
  echo "== job $1 -> hsc_dr4_igmdust_$1  $(date)"
  curl -nsSL -C - -o "hsc_dr4_igmdust_$1.part" "$BASE/$2" \
    && mv "hsc_dr4_igmdust_$1.part" "hsc_dr4_igmdust_$1" \
    && ls -ln "hsc_dr4_igmdust_$1" \
    && file "hsc_dr4_igmdust_$1"
}

fetch 669956 2b3a87f16f886bd2df27a04b67ae7b565094bf47445ef1c7c84de267ef18b5d5
fetch 669957 e65308bfd3458ce5d2460bde3644a31e366a14c1622bd9f97f64723715120a78
# stars: i<22.5, i_extendedness=0 point sources (789 MB fits) - null test
fetch 669937 9f226e7f70e028af89aa4dc3ade68f950bdc70447cac8ec39f17e6d5cebe5fed

echo "HSC_CATALOG_DONE $(date)"
