# 5/20/2026
# VCF filtering of Xipho for selscan
# I need four VCFs: 1) Xingu and 2) Belem populations, and Tapajos + Xin and Tapajos + Bel. Then I can split the Tapajos + Xin
## Check individuals in raw VCF
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz | wc -l
```

# Individuals in Tapajos population
```
/scratch/amonc/xipho_revision/vcf_filtering/xipho_tapajos.txt
```
Tap_xsA08267_Mus
Tap_xsBR163-028_Mus
Tap_xsBR163-145_Mus
Tap_xsBR163-212_Mus
Tap_xsMPDS1217_Mus
Tap_xsMPDS1294_Mus
Tap_xsMSF111_Mus
Tap_xsPIME217_Mus
Tap_xsSER013_Mus
Tap_xsTM005_Mus

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

#### Tapajos file
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Tap_selscan_%J_stdout.txt
#SBATCH --error=Tap_selscan_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Tap_selscan
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

# create a filtered VCF
vcftools \
    --gzvcf FINAL_XIPHO_raw.vcf.gz \
    --keep /scratch/amonc/xipho_revision/vcf_filtering/xipho_tapajos.txt \
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
    --max-missing 1 \
    --recode --stdout | bgzip -c > Tap_selscan.vcf.gz

# index VCF with tabix
tabix Tap_selscan.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz \
    --depth \
    --out Tap_selscan_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz \
    --missing-indv \
    --out Tap_selscan_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 5492509 out of a possible 1115307525 Sites


# Xingu VCF for selscan
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Xin_selscan_%J_stdout.txt
#SBATCH --error=Xin_selscan_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Xin_selscan
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
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --max-missing 1 \
    --recode --stdout | bgzip -c > Xin_selscan.vcf.gz

# index VCF with tabix
tabix Xin_selscan.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz \
    --depth \
    --out Xin_selscan_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz \
    --missing-indv \
    --out Xin_selscan_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 1399203 out of a possible 1115307525 Sites


# Belem VCF for selscan
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=Bel_selscan_%J_stdout.txt
#SBATCH --error=Bel_selscan_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Bel_selscan
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
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --max-missing 1 \
    --recode --stdout | bgzip -c > Bel_selscan.vcf.gz

# index VCF with tabix
tabix Bel_selscan.vcf.gz

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz \
    --depth \
    --out Bel_selscan_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz \
    --missing-indv \
    --out Bel_selscan_missing

cp /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz.tbi /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp /scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
After filtering, kept 3550590 out of a possible 1115307525 Sites


# Quick VCF scaffold counts (requires index)
bcftools index -s Tap_selscan.vcf.gz | wc -l # 345 scaffolds
bcftools index -s Xin_selscan.vcf.gz | wc -l # 294 scaffolds
bcftools index -s Bel_selscan.vcf.gz | wc -l # 333 scaffolds

# Splitting VCFs 
# Ok, I have three VCF files filtered for running with selscan
# These files are the following:
/scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz
/scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz
/scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz

# I need to splice these files into separate VCFs for each scaffold

/scratch/amonc/xipho_revision/selscan

# Tapajos
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Tap_vcfsplit_%J_stdout.txt
#SBATCH --error=Tap_vcfsplit_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=Tap_vcfsplit
#SBATCH --chdir=/scratch/amonc/xipho_revision/selscan
#
#################################################
#!/bin/bash

# Input VCF
VCF=/scratch/amonc/xipho_revision/vcf_filtering/Tap_selscan.vcf.gz

# Output directory
OUTDIR=/scratch/amonc/xipho_revision/selscan/Tap_vcfs

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

# Xingu VCFs split
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Xin_vcfsplit_%J_stdout.txt
#SBATCH --error=Xin_vcfsplit_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=Xin_vcfsplit
#SBATCH --chdir=/scratch/amonc/xipho_revision/selscan
#
#################################################
#!/bin/bash

# Input VCF
VCF=/scratch/amonc/xipho_revision/vcf_filtering/Xin_selscan.vcf.gz

# Output directory
OUTDIR=/scratch/amonc/xipho_revision/selscan/Xin_vcfs

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

# Belem VCFs split
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Bel_vcfsplit_%J_stdout.txt
#SBATCH --error=Bel_vcfsplit_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=Bel_vcfsplit
#SBATCH --chdir=/scratch/amonc/xipho_revision/selscan
#
#################################################
#!/bin/bash

# Input VCF
VCF=/scratch/amonc/xipho_revision/vcf_filtering/Bel_selscan.vcf.gz

# Output directory
OUTDIR=/scratch/amonc/xipho_revision/selscan/Bel_vcfs

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