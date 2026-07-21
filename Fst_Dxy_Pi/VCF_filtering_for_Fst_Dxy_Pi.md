# 4/20/2026
# VCF filtering of Xipho for Fst, Dxy, and Pi (I want invariant sites)

## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

# There is apparently no way to mask genotypes with RGQ < 20 in VCFtools, so I'm doing this in bcftools
# Any genotype with RGQ < 20 is marked as ./.
# Fst, Dxt, Pi filtering
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Fst_filt_%J_stdout.txt
#SBATCH --error=Fst_filt_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=Fst_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

bcftools filter \
  -e 'FORMAT/RGQ < 20' \
  -S . \
  -Oz \
  -o FINAL_XIPHO_raw.RGQ20.vcf.gz \
  FINAL_XIPHO_raw.vcf.gz

# create a filtered VCF containing only invariant sites
vcftools \
--gzvcf FINAL_XIPHO_raw.RGQ20.vcf.gz \
--remove-indv XELEGANS_MPEG75162_Mus \
--exclude-bed /scratch/amonc/xipho_revision/ref/repeat_regions.bed \
--max-maf 0 \
--minDP 5 \
--maxDP 50 \
--max-missing 0.50 \
--recode --stdout | bgzip -c > RGQ20.invariant.vcf.gz

# create a filtered VCF containing only variant sites
vcftools \
    --gzvcf FINAL_XIPHO_raw.RGQ20.vcf.gz \
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
    --max-missing 0.50 \
    --recode --stdout | bgzip -c > Fst_Dxy_Pi_variant.vcf.gz

# index both vcfs using tabix
tabix RGQ20.invariant.vcf.gz
tabix Fst_Dxy_Pi_variant.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
RGQ20.invariant.vcf.gz Fst_Dxy_Pi_variant.vcf.gz \
-O z -o Fst_Dxy_Pi_allsites.vcf.gz
tabix Fst_Dxy_Pi_allsites.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites.vcf.gz \
    --depth \
    --out Fst_Dxy_Pi_allsites_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites.vcf.gz \
    --missing-indv \
    --out Fst_Dxy_Pi_allsites_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Fst_Dxy_Pi_allsites_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```

### Final Fst_Dxy_Pi_allsites VCF

# Count SNPs job
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=count_snps_%J_stdout.txt
#SBATCH --error=count_snps_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=count_snps
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
# Get SNP count of final VCF (merged invariant + variant VCFs)
bcftools view -H Fst_Dxy_Pi_allsites.vcf.gz | wc -l > Fst_Dxy_Pi_snp_count.txt

# Count scaffolds
bcftools index -s Fst_Dxy_Pi_allsites.vcf.gz | wc -l