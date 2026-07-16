# 4/20/2026
# VCF filtering of Xipho for PopCluster

## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

## must remove outgroup

#### PopCluster filtering
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=popcluster_filt_%J_stdout.txt
#SBATCH --error=popcluster_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=popcluster_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

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
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --max-missing 0.75 \
    --thin 10000 \
    --recode --stdout | bgzip -c > PopCluster.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/PopCluster.vcf.gz \
    --depth \
    --out PopCluster_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/PopCluster.vcf.gz \
    --missing-indv \
    --out PopCluster_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/PopCluster.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/PopCluster_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/PopCluster_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 99634 out of a possible 1115307525 Sites

### Final PopCluster VCF

#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=tabix_popC_%J_stdout.txt
#SBATCH --error=tabix_popC_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=tabix_popC
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

tabix PopCluster.vcf.gz

### Get quick scaffold count
bcftools index -s PopCluster.vcf.gz | wc -l

cp /scratch/amonc/xipho_revision/vcf_filtering/PopCluster.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/