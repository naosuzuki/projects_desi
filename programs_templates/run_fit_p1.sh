#!/bin/zsh
# Fit the p1 basis (passive composite + 10 residual vectors) to every DESI galaxy in the v2 rest-frame stores (Nao 2026-11-04).
# Output ~/data/sdssdesi_galfit_p1 -> exdisk3; resumable (files present are skipped). Status step_gal_fit_p1.json, log gal_fit_p1.log.
cd ~/github/projects_desi/programs_templates
OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 -W ignore 32_fit_galaxies.py --basis ~/data/sdssdesi_passC_v3t/empca/galaxy_empca_v3t_5450_sqrtn_p1.fits \
  --out ~/data/sdssdesi_galfit_p1 --tag p1 --workers 6 >> ~/data/sdssdesi_status/gal_fit_p1.log 2>&1
