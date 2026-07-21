# 4/20/2026
# VCF filtering of Xipho for GADMA demographic analysis
## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

## must remove outgroup

#### Allele stat filtering
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=GADMA_filt_%J_stdout.txt
#SBATCH --error=GADMA_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=GADMA_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
    --remove-indv XELEGANS_MPEG75162_Mus \
    --exclude-bed /scratch/amonc/xipho_revision/ref/repeat_regions.bed \
    --remove-filtered FS_SOR_filter \
    --remove-filtered MQ_filter \
    --remove-filtered RPRS_filter \
    --remove-filtered QUAL_filter \
    --remove-indels \
    --minQ 30 \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --max-missing 0.95 \
    --thin 10000 \
    --recode --stdout | bgzip -c > GADMA.vcf.gz

# index VCF with tabix
tabix GADMA.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/GADMA.vcf.gz \
    --depth \
    --out GADMA_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/GADMA.vcf.gz \
    --missing-indv \
    --out GADMA_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/GADMA.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/GADMA.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/GADMA_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/GADMA_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 89028 out of a possible 1115307525 Sites

### Final GADMA VCF



# Quick scaffold count (requires index); 288 scaffolds
bcftools index -s GADMA.vcf.gz | wc -l


# This is a check to get the "raw" number of snps for 33 individuals in the GADMA dataset--Variable X in: L = (X - Y) / X * Nseq
# Need to remove the outgroup XELEGANS_MPEG75162_Mus and invariant sites to get variants
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=GADMA_raw_%J_stdout.txt
#SBATCH --error=GADMA_raw_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=GADMA_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
    --remove-indv XELEGANS_MPEG75162_Mus \
    --mac 1 \
    --recode --stdout | bgzip -c > GADMA.raw.vcf.gz


##
After filtering, kept 33 out of 34 Individuals
Outputting VCF file...
After filtering, kept 40761548 out of a possible 1115307525 Sites