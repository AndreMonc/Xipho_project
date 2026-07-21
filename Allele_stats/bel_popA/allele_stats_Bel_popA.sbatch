#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=175G
#SBATCH --output=allele_stats_Bel_popA_%J_stdout.txt
#SBATCH --error=allele_stats_Bel_popA_%J_stderr.txt
#SBATCH --time=01:00:00
#SBATCH --job-name=Bel_popA_allele_stats
#SBATCH --chdir=/scratch/amonc/xipho_revision/allele_stats/bel_popA
#
#################################################
source ~/miniforge3/etc/profile.d/conda.sh
conda activate allele_stats_env

python allele_stats.py --vcfFile /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf --skipRows 526 \
--windowFile windows.bed --popKey popKey.txt --popA belem \
--popB xingu --popC tapajos
