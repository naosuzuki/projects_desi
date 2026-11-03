#!/bin/zsh
# After run_stack_v3.sh (2026-10-31): luminosity-stack plot (12_lum_stacks.py) and one EM-PCA basis for all redshifts
# (v3 has no z axis; stack pixels with < max(30, 10 % of the stack's best coverage) members dropped). Log: after_stack_v3.log.
# v4 (Nao 2026-11-01): the same stacks, but each stack's blue end trimmed where the S/N < 5 or the flux falls toward the blue
# (qsolib.blue_trim) and the stacks fainter than M_AB(3025) = -18 dropped (host galaxies, not quasars): defaults of 12_lum_stacks.py
# and 06_empca.py; basis qso_empca_v4.fits, plots plots_luminosity/plots_qso_stacks_by_lum_<tag>_*.png, tables csv_files/.
cd ~/github/projects_desi/programs_quasars
S=~/data/sdssdesi_status; P=~/data/sdssdesi_products/quasar
log() { echo "$(date '+%F %T') $*"; }
until grep -q 'v3 stacks finished\|not merging' $S/qso_stack_v3.log 2>/dev/null; do sleep 120; done
grep -q 'not merging' $S/qso_stack_v3.log && { log "stack run incomplete"; exit 1; }
log "luminosity stacks"; T=${QSO_EMTAG:-v4}; R=~/github/projects_desi    # Nao 2026-11-02: PNGs in plots_luminosity/, CSVs in csv_files/
for k in simple ew slope; do python3 12_lum_stacks.py --$( [ $k = simple ] && echo simple || echo "by $k" ) --out $R/plots_luminosity/plots_qso_stacks_by_lum_${T}_$k.png --csv $R/csv_files/plots_qso_stacks_by_lum_${T}_$k.csv; done
F=~/github/papers/26_sdssdesi_spectra/figures    # paper figures (Nao commits ms.tex)
cp $R/plots_luminosity/plots_qso_stacks_by_lum_${T}_simple.png $F/qso_lum_composites.png; cp $R/plots_luminosity/plots_qso_stacks_by_lum_${T}_ew.png $F/qso_lum_composites_ew.png; cp $R/plots_luminosity/plots_qso_stacks_by_lum_${T}_slope.png $F/qso_lum_composites_slope.png
T=${QSO_EMTAG:-v4}     # QSO_EMTAG=v5: run 5 (telluric-corrected members)
log "EM-PCA $T"; OPENBLAS_NUM_THREADS=8 python3 06_empca.py --stacks $P/v3/qso_stacks.npz --tag $T --min-member 30 --min-member-frac 0.1
log "after_stack_v3 finished"
