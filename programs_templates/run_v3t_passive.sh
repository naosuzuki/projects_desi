#!/bin/zsh
# Passive-reference basis of the v3t stacks, version 5450, sqrt(N) weights, coverage cut 10 % (Nao 2026-11-03): tag v3t_5450_sqrtn_p.
# Status: step_gal_v3t_passive.json.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_passive.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 1, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "passive-e0 EM-PCA, version 5450, sqrt N, cover 10 %"
OPENBLAS_NUM_THREADS=8 nice -n 10 python3 -W ignore 31_passive_empca.py --stacks $D/stacks_5450_chain_stacks $G/stacks_5450_chain_stacks \
  --tag v3t_5450_sqrtn_p --out $D/empca --min-cover-frac 0.1 > $S/gal_v3t_passive.log 2>&1 \
  && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_5450_sqrtn_p.fits --two-point >> $S/gal_v3t_passive.log 2>&1
grep -q "^wrote" $S/gal_v3t_passive.log && st 1 done "galaxy_empca_v3t_5450_sqrtn_p" || st 0 failed "see gal_v3t_passive.log"
