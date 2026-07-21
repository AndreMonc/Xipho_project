#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=arg_bed_region
#SBATCH --output=logs/arg_bed_region_%A_%a.out
#SBATCH --error=logs/arg_bed_region_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-575%575

set -euo pipefail

module load GCC/11.3.0

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel

JOB_TABLE="infoTables/region_job_table.txt"
SMC_DIR="/scratch/amonc/xipho_revision/arg_processing/verified_smcs"
LOG_DIR="/scratch/amonc/xipho_revision/arg_processing/verified_logs"
OUT_DIR="./argBedFiles"

line=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR == id + 1 {print}' "$JOB_TABLE")

region=$(echo "$line" | awk -F'\t' '{print $2}')

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "region=$region"

for iter in $(seq 1510 10 2000)
do
    Rscript scripts/smc_to_bed_xipho_parallel.R \
      "$region" \
      "$iter" \
      "$SMC_DIR" \
      "$LOG_DIR" \
      "$OUT_DIR"
done