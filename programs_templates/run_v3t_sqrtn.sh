#!/bin/zsh
# EM-PCA of the v3t stacks with sqrt(N) stack weighting (Quasar thread's recipe; Nao 2026-11-03), both versions G and 5450,
# tags v3t_G_sqrtn / v3t_5450_sqrtn, same settings otherwise as after_v3t.sh. Status: step_gal_v3t_sqrtn.json.
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_sqrtn.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t
st() { echo "{\"done\": $1, \"total\": 2, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "EM-PCA sqrt(N) weighting, G and 5450 in parallel"
for V in G 5450; do
  ( OPENBLAS_NUM_THREADS=4 nice -n 10 python3 -W ignore 10_galaxy_empca.py --stacks $D/stacks_${V}_chain_stacks $G/stacks_${V}_chain_stacks \
      --nvec 10 --niter 25 --min-member-frac 0.5 --min-cover 150 --min-cover-frac 0.2 --stack-weight sqrtn --tag v3t_${V}_sqrtn --out $D/empca \
      > $S/gal_v3t_sqrtn_$V.log 2>&1 && python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_${V}_sqrtn.fits --two-point >> $S/gal_v3t_sqrtn_$V.log 2>&1 ) &
done
wait
n=0; for V in G 5450; do grep -q "^wrote" $S/gal_v3t_sqrtn_$V.log && n=$((n+1)); done
[ $n -eq 2 ] && st 2 done "galaxy_empca_v3t_{G,5450}_sqrtn" || st $n failed "see gal_v3t_sqrtn_*.log"
