#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=Belem_RAiSD_%J_stdout.txt
#SBATCH --error=Belem_RAiSD_%J_stderr.txt
#SBATCH --time=3:00:00
#SBATCH --job-name=Belem_RAiSD
#SBATCH --chdir=/scratch/amonc/xipho_revision/RAiSD/belem
#
#################################################
module load GSL/1.16-goolf-1.4.10

/home/amonc/RAiSD/raisd-master/RAiSD -n bel_run1 -I /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf -s
