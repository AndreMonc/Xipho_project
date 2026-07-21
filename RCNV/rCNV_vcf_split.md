# May 21 2026, for Xiphorhynchus revision

# split the allele stats vcf for running rCNV by chromosome
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=as_vcfsplit_%J_stdout.txt
#SBATCH --error=as_vcfsplit_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=as_vcfsplit
#SBATCH --chdir=/scratch/amonc/xipho_revision/rCNV
#
#################################################
#!/bin/bash

# Input VCF
VCF=/scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz

# Output directory
OUTDIR=/scratch/amonc/xipho_revision/rCNV/vcfs

# Create output directory
mkdir -p "$OUTDIR"

# Loop through scaffolds and extract each into its own VCF
tabix -l "$VCF" | while read scaffold
do
    echo "Processing $scaffold"

    bcftools view \
        -r "$scaffold" \
        -Oz \
        -o "${OUTDIR}/${scaffold}.vcf.gz" \
        "$VCF"

    # index each output VCF
    tabix -p vcf "${OUTDIR}/${scaffold}.vcf.gz"

done
```