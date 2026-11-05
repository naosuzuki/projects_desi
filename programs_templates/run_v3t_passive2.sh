#!/bin/zsh
# p2: passive-reference basis as p1, dropping the stacks with a step on the 3000-3100 A window (|step| > 0.1 -> p2, > 0.05 -> p2b)
# and the pixels below 2000 A (Nao 2026-11-04). Tags v3t_5450_sqrtn_p2, _p2b; status step_gal_v3t_passive2.json.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_passive2.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 2, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "passive-e0 EM-PCA p2 (step 0.1) and p2b (step 0.05), wave-min 2000"
run() { OPENBLAS_NUM_THREADS=3 nice -n 10 python3 -W ignore 31_passive_empca.py --stacks $D/stacks_5450_chain_stacks $G/stacks_5450_chain_stacks \
  --tag v3t_5450_sqrtn_$1 --out $D/empca --min-cover-frac 0.1 --amp unit --wave-min 2000 --max-step $2 > $S/gal_v3t_passive_$1.log 2>&1 \
  && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_5450_sqrtn_$1.fits --two-point >> $S/gal_v3t_passive_$1.log 2>&1; }
run p2 0.1 & run p2b 0.05 & wait
N=$(grep -l "^wrote" $S/gal_v3t_passive_p2.log $S/gal_v3t_passive_p2b.log 2>/dev/null | wc -l | tr -d ' ')
[ "$N" = 2 ] && st 2 done "galaxy_empca_v3t_5450_sqrtn_p2, _p2b" || st $N failed "see gal_v3t_passive_p2*.log"
