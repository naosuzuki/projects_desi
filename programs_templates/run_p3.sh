#!/bin/zsh
# p3: passive-reference basis with e0 = composite of the coefficient-selected passive galaxies ($1 = npz, $2 = tag, $3 = wave-min), Nao 2026-11-04.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_$2.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 1, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "passive-e0 EM-PCA $2 with external e0 $(basename $1)"
OPENBLAS_NUM_THREADS=4 nice -n 10 python3 -W ignore 31_passive_empca.py --stacks $D/stacks_5450_chain_stacks $G/stacks_5450_chain_stacks \
  --tag v3t_5450_sqrtn_$2 --out $D/empca --min-cover-frac 0.1 --amp unit --wave-min $3 --max-step 0.1 --e0 $1 > $S/gal_v3t_$2.log 2>&1 \
  && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_5450_sqrtn_$2.fits --two-point >> $S/gal_v3t_$2.log 2>&1
grep -q "^wrote" $S/gal_v3t_$2.log && st 1 done "galaxy_empca_v3t_5450_sqrtn_$2" || st 0 failed "see gal_v3t_$2.log"
