#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=igmdust_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/home/suzuki/slurmlogs/%x_%j.log
#
# 03b: DR11 sweep download, shard 2 of 2, as a small SLURM job on astrofs4
# (downloads there OK per admin, 2026-11-06). astro4 runs shard 1 of 2
# outside SLURM at the same time; the shards are disjoint by construction.
# This is a network job: never request more cores or memory.
#
# Submit (typed on astrofs4):  sbatch -d singleton 03b_slurm_dr11_download.sbatch
set -u
echo "== $(hostname) $(date)  job $SLURM_JOB_ID"
df -h /data/astrofs2_1 | tail -1 || exit 1
export NPAR=8
exec bash "$HOME/github/projects_desi/programs_igmdust/00b_download_dr11_sweeps.sh" 2 2
