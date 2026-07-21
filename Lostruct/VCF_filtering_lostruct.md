# 4/20/2026
# VCF filtering of Xipho for B allele stats and Q95, U20, and U50
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
#SBATCH --output=allele_stat_filt_%J_stdout.txt
#SBATCH --error=allele_stat_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=allele_stat_filt
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
    --max-missing 0.50 \
    --recode --stdout | bgzip -c > allele_stats.vcf.gz

# index VCF with tabix
tabix allele_stats.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz \
    --depth \
    --out allele_stats_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz \
    --missing-indv \
    --out allele_stats_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/allele_stats_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/allele_stats_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 27012938 out of a possible 1115307525 Sites

### Final allele_stats VCF


# Quick scaffold count (requires index); 414 scaffolds
bcftools index -s allele_stats.vcf.gz | wc -l


# Further filtering for lostruct in order to reduce snp count and speed
# processing time for inversion detection pipeline with lostruct
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=lostruct_filt_%J_stdout.txt
#SBATCH --error=lostruct_filt_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=lostruct_filt
#SBATCH --chdir=/scratch/amonc/xipho_revision/lostruct
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf allele_stats.vcf.gz \
    --max-missing 0.75 \
    --maf 0.05 \
    --recode --stdout | bgzip -c > lostruct.vcf.gz

# index VCF with tabix
tabix lostruct.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz \
    --depth \
    --out lostruct.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz \
    --missing-indv \
    --out lostruct.vcf.gz

cp /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision_old/vcfs/
cp /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision_old/vcfs/
cp /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision_old/vcfs/
cp /scratch/amonc/xipho_revision/lostruct/lostruct.vcf.gz.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision_old/vcfs/


# Now on sooner, so quick install of bcftools

wget https://github.com/samtools/htslib/releases/download/1.24/htslib-1.24.tar.bz2
tar -xvf htslib-1.24.tar.bz2 
module load GCC/14.3.0
module load bzip2/1.0.8-GCCcore-14.3.0
module load XZ/5.8.1-GCCcore-14.3.0 
./configure --prefix=/home/amonc/htslib-1.24
make
make install

liblzma

wget https://github.com/samtools/bcftools/releases/download/1.24/bcftools-1.24.tar.bz2

tar -xvf bcftools-1.24.tar.bz2
module load GCC/14.3.0
./configure --prefix=/home/amonc/bcftools-1.24
make
make install

# Quick scaffold count (requires index); 414 scaffolds
bcftools index -s lostruct.vcf.gz | wc -l