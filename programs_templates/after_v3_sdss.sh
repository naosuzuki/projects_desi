#!/bin/zsh
# After the v3 SDSS re-read (DR8 + DR17): merge the stacks of both versions, chain them, plot by luminosity (Dn4000 type, emission-line class).
# Progress: ~/data/sdssdesi_status/step_gal_v3_sdss_after.json (monitor step "v3 SDSS stacks: merge, chain, plots").
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3_sdss_after.json; D=~/data/sdssdesi_passG_sdss_v3; P=~/github/projects_desi
st() { echo "{\"done\": $1, \"total\": 9, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "waiting for the DR8 and DR17 passes"
while pgrep -f "11_passG_sdss.py --v3" >/dev/null; do sleep 30; done
n=0
for V in G 5450; do
  python3 09_merge_stacks.py --passc $D --acc "acc_${V}_dr*.npz" --out stacks_$V --floor 200 || { st $n failed "merge $V"; exit 1; }
  n=$((n+1)); st $n running "merged $V"
  python3 -W ignore 23_chain_v3.py --version $V --stacks $D/stacks_$V --tables "$D/tables/dr*/*.parquet" || { st $n failed "chain $V"; exit 1; }
  n=$((n+1)); st $n running "chained $V"
  for BY in db ew; do
    python3 -W ignore 21_plot_lum.py --axis $V --by $BY --chain $D/stacks_${V}_chain --out $P/plots_luminosity/plots_v3_sdss_${V}_lum_${BY}.png --title "v3 SDSS, version $V: " || { st $n failed "plot $V $BY"; exit 1; }
    n=$((n+1)); st $n running "plotted $V $BY"
  done
done
st 9 done "SDSS v3: stacks_G, stacks_5450 merged, chained, plotted (plots_v3_sdss_*)"
