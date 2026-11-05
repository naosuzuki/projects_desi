#!/bin/zsh
# Passive-reference basis with unit amplitude (residuals X - e0, as about a mean), version 5450, sqrt N, cover 10 %: tag v3t_5450_sqrtn_p1.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_passive1.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 1, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "passive-e0 EM-PCA, unit amplitude"
OPENBLAS_NUM_THREADS=6 nice -n 10 python3 -W ignore 31_passive_empca.py --stacks $D/stacks_5450_chain_stacks $G/stacks_5450_chain_stacks \
  --tag v3t_5450_sqrtn_p1 --out $D/empca --min-cover-frac 0.1 --amp unit > $S/gal_v3t_passive1.log 2>&1 \
  && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_5450_sqrtn_p1.fits --two-point >> $S/gal_v3t_passive1.log 2>&1
grep -q "^wrote" $S/gal_v3t_passive1.log && st 1 done "galaxy_empca_v3t_5450_sqrtn_p1" || st 0 failed "see gal_v3t_passive1.log"
