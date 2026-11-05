#!/bin/zsh
# After the v3t re-read (telluric-corrected spectra; DESI 4 shards + SDSS dr8/dr17, run_v3t.sh): merge the stacks of both versions and
# both surveys, chain with the stack-measured link (27_chain_v3_stacks.py), luminosity plots (mag and F_lambda), EM-PCA of both versions
# (same settings as after_chain_stacks.sh, tags v3t_G / v3t_5450) and the eigenspectrum plots. The uncorrected v3 products stay as the baseline.
# Progress: ~/data/sdssdesi_status/step_gal_v3t_after.json (monitor step "v3t stacks: merge, stack-link chaining, plots, EM-PCA").
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3t_after.json; D=~/data/sdssdesi_passC_v3t; G=~/data/sdssdesi_passG_sdss_v3t; P=~/github/projects_desi
st() { echo "{\"done\": $1, \"total\": 16, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); n=0; st 0 running "waiting for the v3t passes"
while pgrep -f "08_passC_stack.py --all --v3 --telluric|11_passG_sdss.py --v3 --telluric" >/dev/null; do sleep 60; done
for V in G 5450; do
  [ -s $D/stacks_$V/stacks_z00.fits ] || python3 09_merge_stacks.py --passc $D --acc "acc_${V}_shard*of4.npz" --out stacks_$V --floor 200 || { st $n failed "merge DESI $V"; exit 1; }
  [ -s $G/stacks_$V/stacks_z00.fits ] || python3 09_merge_stacks.py --passc $G --acc "acc_${V}_dr*.npz" --out stacks_$V --floor 200 || { st $n failed "merge SDSS $V"; exit 1; }
  n=$((n+1)); st $n running "merged $V"
  python3 -W ignore 27_chain_v3_stacks.py --version $V --desi $D/stacks_$V --desi-tables "$D/tables/*/*.parquet" --sdss $G/stacks_$V --sdss-tables "$G/tables/dr*/*.parquet" || { st $n failed "chain $V"; exit 1; }
  n=$((n+1)); st $n running "chained $V (stack link)"
  for BY in db ew; do for U in mag flam; do
    python3 -W ignore 21_plot_lum.py --axis $V --by $BY --unit $U --chain $D/stacks_${V}_chain_stacks $G/stacks_${V}_chain_stacks \
      --out $P/plots_luminosity/plots_v3t_${V}_lum_${BY}_${U}.png --title "v3t (telluric-corrected) DESI + SDSS, version $V: " || { st $n failed "plot $V $BY $U"; exit 1; }
    n=$((n+1)); st $n running "plotted $V $BY $U"
  done; done
done
for V in G 5450; do
  OPENBLAS_NUM_THREADS=4 nice -n 10 python3 -W ignore 10_galaxy_empca.py --stacks $D/stacks_${V}_chain_stacks $G/stacks_${V}_chain_stacks \
    --nvec 10 --niter 25 --min-member-frac 0.5 --min-cover 150 --min-cover-frac 0.2 --tag v3t_$V --out $D/empca || { st $n failed "EM-PCA $V"; exit 1; }
  n=$((n+1)); st $n running "EM-PCA $V done"
  python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_$V.fits; python3 24_plot_empca.py $D/empca/galaxy_empca_v3t_$V.fits --two-point
done
st 16 done "v3t: stacks merged, stack-link chaining, plots_v3t_*, EM-PCA galaxy_empca_v3t_{G,5450}"
