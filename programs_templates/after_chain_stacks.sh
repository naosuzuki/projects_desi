#!/bin/zsh
# v3 chaining with the stack-measured link (27_chain_v3_stacks.py, Nao 2026-11-01 option 1), then the luminosity plots (both
# versions, Dn4000 / EW classes, mag and F_lambda) and the EM-PCA of both versions (10_galaxy_empca.py, edge cuts, reference window).
# Progress: ~/data/sdssdesi_status/step_gal_v3_chain_stacks.json
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3_chain_stacks.json; D=~/data/sdssdesi_passC_v3; G=~/data/sdssdesi_passG_sdss_v3; P=~/github/projects_desi
st() { echo "{\"done\": $1, \"total\": 12, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); n=0; st 0 running "chaining G with the stack link"
for V in G 5450; do
  python3 -W ignore 27_chain_v3_stacks.py --version $V --desi $D/stacks_$V --desi-tables "$D/tables/*/*.parquet" --sdss $G/stacks_$V --sdss-tables "$G/tables/dr*/*.parquet" || { st $n failed "chain $V"; exit 1; }
  n=$((n+1)); st $n running "chained $V (stack link)"
  for BY in db ew; do for U in mag flam; do
    python3 -W ignore 21_plot_lum.py --axis $V --by $BY --unit $U --chain $D/stacks_${V}_chain_stacks $G/stacks_${V}_chain_stacks \
      --out $P/plots_luminosity/plots_v3s_${V}_lum_${BY}_${U}.png --title "v3 DESI + SDSS, version $V, stack link: " || { st $n failed "plot $V $BY $U"; exit 1; }
    n=$((n+1)); st $n running "plotted $V $BY $U"
  done; done
done
for V in G 5450; do
  OPENBLAS_NUM_THREADS=4 nice -n 10 python3 -W ignore 10_galaxy_empca.py --stacks $D/stacks_${V}_chain_stacks $G/stacks_${V}_chain_stacks \
    --nvec 10 --niter 25 --min-member-frac 0.5 --min-cover 150 --min-cover-frac 0.2 --tag v3s_$V --out $D/empca || { st $n failed "EM-PCA $V"; exit 1; }
  n=$((n+1)); st $n running "EM-PCA $V done"
  python3 24_plot_empca.py $D/empca/galaxy_empca_v3s_$V.fits; python3 24_plot_empca.py $D/empca/galaxy_empca_v3s_$V.fits --two-point
done
st 12 done "stack-link chaining, plots_v3s_*, EM-PCA galaxy_empca_v3s_{G,5450}"
