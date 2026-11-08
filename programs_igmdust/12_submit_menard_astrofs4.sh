#!/bin/bash
# 12_submit_menard_astrofs4.sh — submit the Menard-replica chain on
# astrofs4.  TYPE THIS YOURSELF in an interactive astrofs4 shell:
#     bash ~/github/projects_desi/programs_igmdust/12_submit_menard_astrofs4.sh
# (never via ssh remote command — the login wrapper duplicates those).
set -eu
SB=~/github/projects_desi/programs_igmdust/13_slurm_menard.sbatch
scan=$(sbatch --parsable -J igmdust_pgs "$SB" scan)
c0=$(sbatch --parsable -J igmdust_pgw_0 --dependency=afterok:$scan "$SB" count0)
c1=$(sbatch --parsable -J igmdust_pgw_1 --dependency=afterok:$scan "$SB" count1)
fin=$(sbatch --parsable -J igmdust_pgf --dependency=afterok:$c0:$c1 "$SB" finish)
echo "submitted: scan=$scan count0=$c0 count1=$c1 finish=$fin"
squeue -u suzuki
