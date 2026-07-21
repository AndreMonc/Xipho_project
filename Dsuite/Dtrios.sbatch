#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Dtrios_%J_stdout.txt
#SBATCH --error=Dtrios_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Dtrios
#SBATCH --chdir=/scratch/amonc/xipho_revision/dsuite
#
#################################################

module load GCCcore/11.3.0

Dsuite Dtrios /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz /scratch/amonc/xipho_revision/dsuite/SETS.txt
