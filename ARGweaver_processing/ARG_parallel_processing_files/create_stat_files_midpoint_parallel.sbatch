#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=arg_stat_mid
#SBATCH --output=logs/arg_stat_mid_%A_%a.out
#SBATCH --error=logs/arg_stat_mid_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-573%573

set -euo pipefail

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel

JOB_TABLE="infoTables/region_tree_job_table_midpoint.txt"
TREE_DIR="./argTreeFiles_midpoint"
STATS_DIR="./argStats_midpoint"

line=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR == id + 1 {print}' "$JOB_TABLE")

region=$(echo "$line" | awk -F'\t' '{print $2}')

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "region=$region"

for iter in $(seq 1510 10 2000)
do
    tree_file=$(find "${TREE_DIR}/${region}" -maxdepth 1 -type f -name "${region}.*.${iter}.tre.gz" | head -1)

    if [[ -z "$tree_file" ]]; then
        echo "Missing TRE for region=${region}, iter=${iter}"
        continue
    fi

    Rscript scripts/tre_to_stats_xipho_parallel.R \
        "$tree_file" \
        "${STATS_DIR}/${region}"
done

echo "Finished stats for $region"