#!/bin/zsh
# After the v3 DESI re-read (4 shards): merge the stacks of both versions, chain them, plot by luminosity together with the
# SDSS v3 chains (Dn4000 type, emission-line class; absolute AB mag and F_lambda at 10 pc).
# Progress: ~/data/sdssdesi_status/step_gal_v3_desi_after.json (monitor step "v3 DESI stacks: merge, chain, plots").
cd ~/github/projects_desi/programs_templates
S=~/data/sdssdesi_status; J=$S/step_gal_v3_desi_after.json; D=~/data/sdssdesi_passC_v3; G=~/data/sdssdesi_passG_sdss_v3; P=~/github/projects_desi
st() { echo "{\"done\": $1, \"total\": 12, \"status\": \"$2\", \"started\": $T0, \"note\": \"$3\"}" > $J; }
T0=$(date +%s); st 0 running "waiting for the 4 DESI shards"
while pgrep -f "08_passC_stack.py --all --v3" >/dev/null; do sleep 30; done
n=0
for V in G 5450; do
  [ -s $D/stacks_$V/stacks_z00.fits ] || python3 09_merge_stacks.py --passc $D --acc "acc_${V}_shard*of4.npz" --out stacks_$V --floor 200 || { st $n failed "merge $V"; exit 1; }
  n=$((n+1)); st $n running "merged $V"
  python3 -W ignore 23_chain_v3.py --version $V --stacks $D/stacks_$V --tables "$D/tables/*/*.parquet" || { st $n failed "chain $V"; exit 1; }
  n=$((n+1)); st $n running "chained $V"
  for BY in db ew; do
    for U in mag flam; do
      python3 -W ignore 21_plot_lum.py --axis $V --by $BY --unit $U --chain $D/stacks_${V}_chain $G/stacks_${V}_chain \
        --out $P/plots_luminosity/plots_v3_${V}_lum_${BY}_${U}.png --title "v3 DESI + SDSS, version $V: " || { st $n failed "plot $V $BY $U"; exit 1; }
      n=$((n+1)); st $n running "plotted $V $BY $U"
    done
  done
done
st 12 done "DESI v3: stacks_G, stacks_5450 merged, chained; plots_v3_{G,5450}_lum_{db,ew}_{mag,flam}.png with SDSS"
