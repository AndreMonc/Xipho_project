e# June 14, 2026

LoStruct analysis to identify putative inversions 
Method developed by Li & Ralph 2019
Following pipeline in Huang et al. 2020, Molecular Ecology

cd /scratch/amonc/xipho_revision
mkdir lostruct

# Working directory for lostruct analysis
/scratch/amonc/xipho_revision/lostruct


# VCF input:
/scratch/amonc/xipho_revision/lostruct/allele_stats.vcf.gz
/scratch/amonc/xipho_revision/lostruct/allele_stats.vcf.gz.tbi

# Convert VCF to BCF format with BCFTOOLS
##bcftools_viewVersion=1.23.1+htslib-1.23.1

# Conversion script
#!/bin/bash
#SBATCH --job-name=vcf_to_bcf
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=vcf_to_bcf_%j.out
#SBATCH --error=vcf_to_bcf_%j.err

set -euo pipefail

WORKDIR="/scratch/amonc/xipho_revision/lostruct"

cd "$WORKDIR"

echo "Starting conversion at $(date)"

bcftools view \
    -O b \
    -o lostruct.bcf \
    lostruct.vcf.gz

echo "Indexing BCF file..."

bcftools index lostruct.bcf

echo "Finished at $(date)"

ls -lh lostruct.bcf lostruct.bcf.csi

# Load R module
module load R/4.2.1-foss-2022a
.libPaths() 
install.packages("data.table") # follow prompts to use a personal library

devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

# > packageVersion("lostruct")
[1] ‘0.0.0.9000’

# Next time just the following commands to start up R and lostruct
module load R/4.2.1-foss-2022a
library(lostruct)
