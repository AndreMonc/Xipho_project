#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=arg_branch_fst
#SBATCH --output=logs/arg_branch_fst_%j.out
#SBATCH --error=logs/arg_branch_fst_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail

cd /scratch/amonc/xipho_revision/arg_processing/arg_based_fst

mkdir -p logs per_region_branch_fst

module purge

source /home/amonc/miniforge3/etc/profile.d/conda.sh
conda activate arg_fst_env

unset PYTHONPATH
unset PYTHONHOME

python tre_to_branch_fst_xipho.py
