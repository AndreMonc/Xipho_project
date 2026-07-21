# 4/20/2026
# VCF filtering of Xipho for D statistics and ABBA-BABA (Fd across genome)

## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

#### Dstats filtering
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=Dstats_filt_%J_stdout.txt
#SBATCH --error=Dstats_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Dstats_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
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
    --maxDP 100 \
    --minGQ 20 \
    --max-missing 0.50 \
    --recode --stdout | bgzip -c > Dstats.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz \
    --depth \
    --out Dstats_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz \
    --missing-indv \
    --out Dstats_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Dstats_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Dstats_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 32015475 out of a possible 1115307525 Sites

### Final Dstats VCF

## Get index
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=tabix_dstats_%J_stdout.txt
#SBATCH --error=tabix_dstats_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=tabix_dstats
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

tabix Dstats.vcf.gz


### Get quick scaffold count
bcftools index -s Dstats.vcf.gz | wc -l

cp /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/