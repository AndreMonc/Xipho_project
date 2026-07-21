#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=Xingu_RAiSD_%J_stdout.txt
#SBATCH --error=Xingu_RAiSD_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=Xingu_RAiSD
#SBATCH --chdir=/scratch/amonc/xipho_revision/RAiSD/xingu
#
#################################################
module load GSL/1.16-goolf-1.4.10

/home/amonc/RAiSD/raisd-master/RAiSD -n xin_run1 -I /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD.vcf -s
