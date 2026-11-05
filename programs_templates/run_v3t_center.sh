#!/bin/zsh
# Centered (classic) PCA of the v3t stacks, version 5450, sqrt(N) weights (Nao 2026-11-03): tag v3t_5450_sqrtn_c. Status: step_gal_v3t_center.json.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_center.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 1, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "centered EM-PCA, version 5450, sqrt N"
OPENBLAS_NUM_THREADS=8 nice -n 10 python3 -W ignore 10_galaxy_empca.py --stacks $D/stacks_5450_chain_stacks $G/stacks_5450_chain_stacks \
  --nvec 10 --niter 25 --min-member-frac 0.5 --min-cover 150 --min-cover-frac 0.2 --stack-weight sqrtn --center --tag v3t_5450_sqrtn_c --out $D/empca \
  > $S/gal_v3t_center.log 2>&1 && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_5450_sqrtn_c.fits --two-point >> $S/gal_v3t_center.log 2>&1
grep -q "^wrote" $S/gal_v3t_center.log && st 1 done "galaxy_empca_v3t_5450_sqrtn_c" || st 0 failed "see gal_v3t_center.log"
