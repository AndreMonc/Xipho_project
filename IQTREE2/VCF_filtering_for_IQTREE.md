# 4/20/2026
# VCF filtering of Xipho for IQTREE

## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

#### IQ-TREE2 filtering
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=IQTREE_filt_%J_stdout.txt
#SBATCH --error=IQTREE_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=IQTREE_filt
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
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 100 \
    --minGQ 20 \
    --max-missing 0.75 \
    --thin 100 \
    --recode --stdout | bgzip -c > IQTREE.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz \
    --depth \
    --out IQTREE_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz \
    --missing-indv \
    --out IQTREE_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/IQTREE_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/IQTREE_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 4120529 out of a possible 1115307525 Sites

### Final IQTREE VCF
/scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz #also on ourdisk

# Index VCF
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=tabix_IQTREE_%J_stdout.txt
#SBATCH --error=tabix_IQTREE_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=tabix_IQTREE
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

tabix IQTREE.vcf.gz

# Count scaffolds (this command requires an index .tbi file)
bcftools index -s IQTREE.vcf.gz | wc -l # 372