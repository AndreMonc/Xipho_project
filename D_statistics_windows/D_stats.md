# 5/25/2026

# Running to get D statistics across the genome

# VCF file
/scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz

# Get needed file (Simon Martin's pipeline)
wget -O parseVCF.py \
https://raw.githubusercontent.com/simonhmartin/genomics_general/master/VCF_processing/parseVCF.py

wget -O ABBABABAwindows.py \
https://raw.githubusercontent.com/simonhmartin/genomics_general/refs/heads/master/ABBABABAwindows.py

wget -O genomics.py \
https://raw.githubusercontent.com/simonhmartin/genomics_general/refs/heads/master/genomics.py

# convert my VCF files to geno files

#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Dstats_%J_stdout.txt
#SBATCH --error=Dstats_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Dstats
#SBATCH --chdir=/scratch/amonc/xipho_revision/Dstats
#
#################################################

source ~/.bashrc
conda activate genomics_general

python /scratch/amonc/xipho_revision/Dstats/parseVCF.py -i /scratch/amonc/xipho_revision/Dstats/Dstats.vcf.gz -o Dstats.geno.gz


# run ABBABABAwindows.py to get window data

#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=Dstats2_%J_stdout.txt
#SBATCH --error=Dstats2_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Dstats2
#SBATCH --chdir=/scratch/amonc/xipho_revision/Dstats
#
#################################################

source ~/.bashrc
conda activate genomics_general

python /scratch/amonc/xipho_revision/Dstats/ABBABABAwindows.py -w 10000 -m 100 -g Dstats.geno.gz -o xiph.Dstats.csv -f phased -T 4 -P1 Bel -P2 Xin -P3 Tap -O Out --popsFile xiph_pops.txt --writeFailedWindows


# Plotting summary stats
# https://speciationgenomics.github.io/sliding_windows/


s


cp /scratch/amonc/xipho_revision/Dstats