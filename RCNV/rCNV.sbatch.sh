#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=rCNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --array=1-414%50
#SBATCH --chdir=/scratch/amonc/xipho_revision/rCNV
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

set -euo pipefail

module load R/4.3.2-gfbf-2023a

VCF_LIST=/scratch/amonc/xipho_revision/rCNV/vcf_list.txt
OUTDIR=/scratch/amonc/xipho_revision/rCNV/results
LOGDIR=/scratch/amonc/xipho_revision/rCNV/logs

mkdir -p "$OUTDIR" "$LOGDIR"

VCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$VCF_LIST")
SCAFFOLD=$(basename "$VCF" .vcf.gz)

exec > "${LOGDIR}/${SCAFFOLD}.out" 2> "${LOGDIR}/${SCAFFOLD}.err"

echo "================================"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scaffold: $SCAFFOLD"
echo "VCF: $VCF"
echo "================================"

Rscript \
/scratch/amonc/xipho_revision/rCNV/scripts/run_rCNV_one_chrom.R \
"$VCF" \
"$OUTDIR"
