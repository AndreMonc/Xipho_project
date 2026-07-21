# 4/8/2026
# VCF filtering of Xipho for ReLERNN

## JUST the Tapajos population
#### Remove all non-Tapajos populations from dataset (34-->10 individuals).

## Check individuals in raw VCF, get list of Tap indivs
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz

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
```
## Add above individuals to a xipho_tapajos.txt file

# Download VCFtools v0.1.17
wget https://github.com/vcftools/vcftools/releases/download/v0.1.17/vcftools-0.1.17.tar.gz
tar -xf vcftools-0.1.17.tar.gz
./configure --prefix=/home/amonc
make
make install

# STEP 1: Initial filter to correct set of individuals
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Relernn_filt_step1_%J_stdout.txt
#SBATCH --error=Relernn_filt_step1_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Relernn_filt_step1
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz \
    --keep /scratch/amonc/xipho_revision/vcf_filtering/xipho_tapajos.txt \
    --recode --stdout | bgzip -c > xipho_relernn_tapajos.vcf.gz
```

## explore filter flags
```
bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' xipho_relernn_tapajos.vcf.gz  | head -10000
```
## number of ind: 10
```
bcftools query -l xipho_relernn_tapajos.vcf.gz | wc -l
```

#### Ok, everything checks out so far
# STEP 2: Bulk filtering for ReLERNN

# Do a fullfilter for ReLERNN
# Then identify scaffolds with <250 snps and then remove those along with sex chroms and small chroms (step 3)
#### minQ filter is redundant with the QUAL_filter GATK filter flag
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=ReLERNN_filt_pt1_%J_stdout.txt
#SBATCH --error=ReLERNN_filt_pt1_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=ReLERNN_filt_pt1
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf xipho_relernn_tapajos.vcf.gz \
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
    --recode --stdout | bgzip -c > RelERNN_biallelic_snps_pt1.vcf.gz
```
After filtering, kept 17596333 out of a possible 1115307525 Sites

## bcftools get uniq chromosome names
```
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_pt1.vcf.gz | uniq

```
## Check filter flags
```
bcftools query -f '%FILTER\n' RelERNN_biallelic_snps_pt1.vcf.gz
```

# Count SNPs in each chromosome
# This needs to be run on the Step 1 filtering VCF file for RelERNN (where the 250-snp threshold matters)
# Once I know those 250-snp scaffolds, I can drop those for use in final ReLERNN VCF and ARGweaver VCF

## Count # of snps per chromosome
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=snps_per_scaffold_%J_stdout.txt
#SBATCH --error=snps_per_scaffold_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=snps_per_scaffold
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_pt1.vcf.gz \
| awk '{count[$1]++} END {for (scaf in count) print scaf"\t"count[scaf]}' \
| sort -k1,1 \
> snps_per_scaffold.txt
```

# just other double checks on snp counts
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_pt1.vcf.gz | awk '{count[$1]++} END {for (c in count) print c, count[c]}'
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_pt1.vcf.gz | sort | uniq -c

## Scaffolds with less than 250 snps
/scratch/amonc/xipho_revision/vcf_filtering/scaff_ft250SNPs.bed


# Remove "bad" scaffolds: small (<110 kb), sex (Z and W), and few snps (< 250)
SEE ARGweaver file for merging of these bed files
Resulting merged bed file:
/scratch/amonc/xipho_revision/vcf_filtering/bad_scaffolds.bed

# STEP 3: Remove all "bad scaffolds" from ReLERNN before running

```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=exclude_bad_scaffs_%J_stdout.txt
#SBATCH --error=exclude_bad_scaffs_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=exclude_bad_scaffs
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf RelERNN_biallelic_snps_pt1.vcf.gz \
    --exclude-bed /scratch/amonc/xipho_revision/vcf_filtering/bad_scaffolds.bed \
    --recode --stdout | bgzip -c > RelERNN_biallelic_snps_tapajos.vcf.gz
```
After filtering, kept 17274416 out of a possible 17596333 Sites

### Final ReLERNN VCF
/scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz

# Save copy of VCF in ourdisk permanent storage
cp /scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/

## Check individuals again
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz
```

## Get final depth and missingness per individual
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=ReLERNN_depth_missing_%J_stdout.txt
#SBATCH --error=ReLERNN_depth_missing_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=ReLERNN_depth_missing
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz \
    --depth \
    --out ReLERNN_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz \
    --missing-indv \
    --out ReLERNN_missing
```
Outputs:


### Final ReLERNN Genome (bed) files
```
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_tapajos.vcf.gz | sort -u
bcftools query -f '%CHROM\n' RelERNN_biallelic_snps_tapajos.vcf.gz | awk '!seen[$0]++' | wc -l
```
So, I have 169 scaffolds in my relerrn VCF file. Now I will just make a .bed file for those scaffolds.
```
awk 'NR==FNR {len[$1]=$2; next} {print $1 "\t0\t" len[$1]}' scaff_len.txt relernn_scaffs.txt > ReLERNN_scaffolds.bed
```
Then I did a quick sort so scaffolds in numerical order 1....232 (169 scaffolds total)

```
/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed
```

# Get ReLERNN index file
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=ReLERNN_tabix_%J_stdout.txt
#SBATCH --error=ReLERNN_tabix_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=ReLERNN_tabix
#SBATCH --chdir=/ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs
#
#################################################
module load GCCcore/11.3.0

tabix ReLERNN_biallelic_snps_tapajos.vcf.gz