# 4/20/2026
# VCF filtering of Xipho for RAiSD (separate VCFs for Xingu and Belem populations)
## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

# create .txt file for Xingu individuals (no header line)
```
/scratch/amonc/xipho_revision/vcf_filtering/xin_indiv.txt
```
Xin_xsCOUFT0424_Mus
Xin_xsFTA012_Mus
Xin_xsGAPTO037_Mus
Xin_xsGAPTO271_Mus
Xin_xsGAPX047_Mus
Xin_xsMAYA066_Mus
Xin_xsMOP011_Mus
Xin_xsMRJ498_Mus
Xin_xsTP36025_Toe
Xin_xsTP36276_Toe
Xin_xsTP48649_Toe
Xin_xsTP81164_Toe
Xin_xsUHE455_Mus

# # create .txt file for Belem individuals (no header line)
```
/scratch/amonc/xipho_revision/vcf_filtering/bel_indiv.txt
```
Bel_xsFRC041_Mus
Bel_xsGUR156_Mus
Bel_xsLCA23_Mus
Bel_xsMLV157_Mus
Bel_xsRDP017_Mus
Bel_xsREBIO002_Mus
Bel_xsTP32151_Toe
Bel_xsTP37350_Toe
Bel_xsTP38597_Toe
Bel_xsTP51965_Toe

#### Xingu file
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Xingu_RAiSD_%J_stdout.txt
#SBATCH --error=Xingu_RAiSD_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Xingu_RAiSD
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
    --keep /scratch/amonc/xipho_revision/vcf_filtering/xin_indiv.txt \
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
    --max-missing 0.5 \
    --recode --stdout | bgzip -c > Xingu_RAiSD.vcf.gz

# index VCF with tabix
tabix Xingu_RAiSD.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD.vcf.gz \
    --depth \
    --out Xingu_RAiSD_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD.vcf.gz \
    --missing-indv \
    --out Xingu_RAiSD_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xingu_RAiSD_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 14183307 out of a possible 1115307525 Sites


### Final Xingu_RAiSD VCF


# Quick Xingu VCF scaffold count (requires index); 410 scaffolds
bcftools index -s Xingu_RAiSD.vcf.gz | wc -l

#### Belem file
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Belem_RAiSD_%J_stdout.txt
#SBATCH --error=Belem_RAiSD_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Belem_RAiSD
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
    --keep /scratch/amonc/xipho_revision/vcf_filtering/bel_indiv.txt \
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
    --max-missing 0.5 \
    --recode --stdout | bgzip -c > Belem_RAiSD.vcf.gz

# index VCF with tabix
tabix Belem_RAiSD.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf.gz \
    --depth \
    --out Belem_RAiSD_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf.gz \
    --missing-indv \
    --out Belem_RAiSD_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 10275610 out of a possible 1115307525 Sites

### Final Belem_RAiSD VCF


# Quick Xingu VCF scaffold count (requires index); 409 scaffolds
bcftools index -s Belem_RAiSD.vcf.gz | wc -l