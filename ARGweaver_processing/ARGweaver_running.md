# ARGweaver begins again!
## Reviving this file on 8 April 2026, to run ARGweaver with new Xiphorhynchus reference VCF (to scaffolds)
## Will update what I need to

## 4/8/26
## Download ARGweaver
```
module load GCCcore/11.3.0 #I had a GCC issue on OSCER until I updated it; then installation worked
git clone https://github.com/CshlSiepelLab/ARGweaver.git 
cd ARGweaver
make
```

## Download SAMtools
```mkdir SAMtools
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=/home/amonc/ #adds to my main home bin folder
make
make install # This works when I use the prefix option after ./configure!! Yay!
```

## Download bedops
## /scratch/a_monc/postdoc/bedops/bedops
```mkdir bedops
cd bedops
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
```

## Download PHAST (PHAST 1.9.5); April 8, 2026
```
git clone https://github.com/CshlSiepelLab/phast.git
module load CMake/3.24.3-GCCcore-11.3.0
cmake -S . -B build -DCMAKE_EXE_LINKER_FLAGS="-lm"
cmake --build build
cmake --install build --prefix /home/amonc/
```

## Download htslib
```
mkdir htslib
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xf htslib-1.21.tar.bz2
cd htslib-1.21
./configure --prefix=/scratch/a_monc/postdoc/htslib-1.21
make 
make install
```

# Download bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
tar -zxvf bedtools-2.31.1.tar.gz
cd bedtools2
make

## My VCF file is the raw output VCF from snpArcher (no filtering done, just have flags for the GATK best practices sites)
## All sites VCF (invariant sites included); 1115307525 SNPs total

`/scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz`

## VCF Filtering
I want to remove the X. elegans outgroup individual, leaving 33 individuals
I should remove all scaffolds less than 110kb from the VCF. That way ARGweaver will not deal with them at all. These will go anyways during ARG processing.
I want to also remove sex chromosome scaffolds (Z and W)
I want to remove scaffolds with few snps (<250)--recommendation from ReLERNN. Do this once I have a filtered VCF ready for input to ReLERNN (to get more accurate sense on snp count before removing scaffolds)


# In summary, remove outgroup and "bad" scaffolds: small (<110 kb), sex (Z and W), and few snps (< 250)
Outgroup individual: XELEGANS_MPEG75162_Mus

Small scaffolds (from .fai file; 244 scaffolds: 
/scratch/amonc/xipho_revision/vcf_filtering/small_scaffolds_110kb.bed

Sex chromosome scaffolds (identified via RagTag; 33 scaffolds):
/scratch/amonc/xipho_revision/vcf_filtering/sex_chrom_scaffolds.bed

Few-snp scaffolds (put on pause until I run initial filtering for ReLERNN; 151 scaffolds):
/scratch/amonc/xipho_revision/vcf_filtering/scaff_ft250SNPs.bed


# Merge bed files of "bad" scaffolds

Manually merged with some excel work:
```
/scratch/amonc/xipho_revision/vcf_filtering/bad_scaffolds.bed
```
(Other double checking:
tail -n +2 small_scaffolds_110kb.bed > tmp1.bed
tail -n +2 sex_chrom_scaffolds.bed > tmp2.bed
tail -n +2 scaff_ft250SNPs.bed > tmp3.bed

cat tmp1.bed tmp2.bed tmp3.bed \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> merged.bed)

# For creating ARGweaver-ready VCF: exclude outgroup, bad scaffolds, and scaffolds that did not make it into final ReLERNN VCF

I know the scaffolds in the full all-sites VCF (477 scaffolds):
/scratch/amonc/xipho_revision/vcf_filtering/xipho_scaffs_all.bed

I also know the scaffolds retained in the ReLERNN VCF (169 scaffolds):
/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed

Now, I want to create a bed file of the scaffolds that are in xipho_scaffs_all.bed but not in ReLERNN_scaffolds.bed (should be 477 - 169 = 308 scaffolds):
awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' ReLERNN_scaffolds.bed xipho_scaffs_all.bed > non_ReLERNN_scaffs.bed

# Once I have final ReLERNN VCF scaffolds, I can prepare the ARGweaver VCF
See ReLERNN VCF filtering .md file for more info

# Filtering to get ARGweaver master VCF file (before chunking in 2 MB regions)
# With known ReLERNN scaffolds
(echo -e "chrom\tchromStart\tchromEnd"; cat ReLERNN_scaffolds.bed) > ReLERNN_scaffolds_with_header.bed

```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=ARG_vcf_master_%J_stdout.txt
#SBATCH --error=ARG_vcf_master_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=ARG_vcf_master
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz \
    --remove-indv XELEGANS_MPEG75162_Mus \
    --bed /scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds_with_header.bed \
    --recode --stdout | bgzip -c > ARGweaver_master.vcf.gz

tabix -p vcf ARGweaver_master.vcf.gz
```
## Check scaffolds output (to do on April 15, 2026)
tabix -l ARGweaver_master.vcf.gz | wc -l # 169 scaffolds
tabix -l ARGweaver_master.vcf.gz > ARGweaver_master_scaffs.txt

## VCF which I will then break up into 2MB chunks prior to running ARGweaver:
```
/scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_master.vcf.gz
```

# Save copy of VCF in ourdisk permanent storage
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=copy_ARG_vcf_%J_stdout.txt
#SBATCH --error=copy_ARG_vcf_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=copy_ARG_vcf
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
cp /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_master.vcf.gz /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
```
# Get depth and missingness
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=argweaver_stats_%J_stdout.txt
#SBATCH --error=argweaver_stats_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=argweaver_stats
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
module load GCCcore/11.3.0

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_master.vcf.gz \
    --depth \
    --out argweaver_depth

vcftools \
    --gzvcf /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_master.vcf.gz \
    --missing-indv \
    --out argweaver_missing

# send over to ourdisk folder
cp -r /scratch/amonc/xipho_revision/vcf_filtering/argweaver_depth.idepth /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/
cp -r /scratch/amonc/xipho_revision/vcf_filtering/argweaver_missing.imiss /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/vcfs/


## MASKS ########################################################################################
1. GATK poor sites mask--For ARGweaver
/scratch/amonc/xipho_revision/vcf_filtering/GATK_poor_sites.bed
2. Repeat mask--For ReLERNN and ARGweaver
/scratch/amonc/xipho_revision/ref/repeat_regions.bed
3. Site quality mask (qual <30 for variant sites)--For ARGweaver (dealing with genotype qualities directly in ARGweaver command)
/scratch/amonc/xipho_revision/vcf_filtering/lowQUAL_variant_regions.bed
4. Regions without recombination rate estimates (<50 snp for window)--Known after running ReLERNN
/scratch/amonc/xipho_revision/relernn/ReLERNN_ft50sites.bed

## 1. I want to create mask for poor GATK sites 
## What is the GATK-flagged poor-quality sites situation on snpArcher output?
## Below I list all filters in the raw vcf file, so "." is apparently equivalent to PASS
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=check_filters_%J_stdout.txt
#SBATCH --error=check_filters_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=check_filters
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
bcftools query -f '%FILTER\n' FINAL_XIPHO_raw.vcf.gz | sort -u > uniq.xipho.GATK.filters.txt
```
## Output:
.
FS_SOR_filter
FS_SOR_filter;MQ_filter
FS_SOR_filter;MQ_filter;RPRS_filter
FS_SOR_filter;RPRS_filter
MQ_filter
MQ_filter;RPRS_filter
QUAL_filter
RPRS_filter

## GATK poor quality mask
The command below creates a bed file for regions for sites with any GATK filter flag except "." (PASS) and "QUAL_filter" (which I am taking care of separately, since I only want to remove "QUAL_filter"-flagged sites that are variant).
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=GATK_bed_%J_stdout.txt
#SBATCH --error=GATK_bed_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=GATK_bed
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################
bcftools query \
  -i 'FILTER!="." && FILTER!~"QUAL_filter"' \
  -f '%CHROM\t%POS0\t%POS\n' \
  /scratch/amonc/xipho_revision/vcf_filtering/FINAL_XIPHO_raw.vcf.gz \
| sort -k1,1 -k2,2n \
| bedtools merge -i - > GATK_poor_sites.bed
```

# 2. Repeat region bed
SEPARATE .MD JUST FOR THIS PROCESS
/scratch/amonc/xipho_revision/ref/repeat_regions.bed

# 3. VCF minimum variant quality filter (using a bed-masking approach because I am not sure if --vcf-min-qual  <minQualScore> in ARGweaver would eliminate the invariant sites that don't have a qual score, which would be overly restrictive)
# This bed file only masks variant regions with QUAL <30 because QUAL for invariant sites is just too restrictive (confidence in a site being fully homozygous)
# This is equivalent to using the QUAL filter (from GATK) only for the variant sites
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=Qual_bed_%J_stdout.txt
#SBATCH --error=Qual_bed_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=Qual_bed
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#################################################
bcftools query \
  -i 'ALT!="." && (QUAL<30 || QUAL=".")' \
  -f '%CHROM\t%POS0\t%POS\n' \
  FINAL_XIPHO_raw.vcf.gz \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> lowQUAL_variant_regions.bed

---

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=Qual_bedv2_%J_stdout.txt
#SBATCH --error=Qual_bedv2_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=Qual_bedv2
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#################################################
bcftools query \
  -i 'ALT!="." && FILTER~"QUAL_filter"' \
  -f '%CHROM\t%POS0\t%POS\n' \
  FINAL_XIPHO_raw.vcf.gz \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> QUAL_filtered_variant_regions.bed

# Results for these two methods look the same
# Going with first file:
/scratch/amonc/xipho_revision/vcf_filtering/lowQUAL_variant_regions.bed

4. Regions without recombination rate estimates (<50 snp for window)--Known after running ReLERNN. (See ReLERNN notes for details)
/scratch/amonc/xipho_revision/relernn/ReLERNN_ft50sites.bed


######################################
1. GATK poor sites mask--For ARGweaver
/scratch/amonc/xipho_revision/vcf_filtering/GATK_poor_sites.bed
2. Repeat mask--For ReLERNN and ARGweaver
/scratch/amonc/xipho_revision/ref/repeat_regions.bed
3. Site quality mask (qual <30 for variant sites)--For ARGweaver (dealing with genotype qualities directly in ARGweaver command)
/scratch/amonc/xipho_revision/vcf_filtering/lowQUAL_variant_regions.bed
4. Regions without recombination rate estimates (<50 snp for window)--Known after running ReLERNN
/scratch/amonc/xipho_revision/relernn/ReLERNN_ft50sites.bed

# Merging bed files to create ARGweaver mask
cat \
  /scratch/amonc/xipho_revision/vcf_filtering/GATK_poor_sites.bed \
  /scratch/amonc/xipho_revision/ref/repeat_regions.bed \
  /scratch/amonc/xipho_revision/vcf_filtering/lowQUAL_variant_regions.bed \
  /scratch/amonc/xipho_revision/relernn/ReLERNN_ft50sites.bed \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> ARGweaver_mask.bed

###########################################################################

# Final ARGweaver mask
/scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed



###############################################################################
# ARGweaver commands

--vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \ # remove low and high depths and poor genotype qualities (for variant and invariant sites)
--maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \ # masks GATK poor sites, repeats, low qual variants, and sites with missing recombination rate estimates
--mask-cluster a,b \ # mask any region of length "b" that has "a" variant sites (possibly indicating alignment errors or mutational hotspots). Using a=2,b=5
--mask-Ns 33 \ # mask any regions with 33 missing lineages (out of 66 total in Xiphorhynchus dataset)
--mutrate 4.6e-9 \ # mutation rate after Ficedula study (Smeds et al., 2016), units of expected mutations per base per generation
--recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \ # estimate from ReLERNN
--popsize 391476 \ #Averaged from values in Thom GBE Table S16
--compress-seq 5 \ #site compression, compression results in a loss of resolution but speeds up processing time. Use conservatively.
--ntimes 20 \ #number of time points, default is 20
--maxtime 1e7 \ #best to have an old maxtime; same value as used in Sporophila paper (Hejase et al. 2020)
--delta 0.005 \ #default value of 0.01 (larger values give more resolution of recent time point sampling)--inspect discrete time values upon running. (Sporophila paper used 0.005)
--sample-step 10 \ # default value is 10
--iters 3000 \ # default is 1000, but I want to be sure of convergence


## 15 April 2026, split the ARGweaver VCF generated earlier into 2MB sliding windows with 100kb overlaps

```
/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed
bedtools makewindows -b /scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed -w 2000000 -s 1900000 > ARGweaver_windows.bed
```

# Here is the final windows file, that has the windows with recombination rates from ReLERNN:
/scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_windows.bed # 634 windows

# Create one VCF per bed file window. Each VCF will be run separately with ARGweaver.
```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=make_vcfs_%J_stdout.txt
#SBATCH --error=make_vcfs_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=make_vcfs
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#################################################

mkdir -p /scratch/amonc/xipho_revision/argweaver/vcfs

while IFS=$'\t' read -r chrom start end
do
    out="/scratch/amonc/xipho_revision/argweaver/vcfs/${chrom}_${start}_${end}.vcf.gz"

    bcftools view \
        -r "${chrom}:$((start+1))-${end}" \
        -Oz \
        -o "$out" \
        /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_master.vcf.gz

    tabix -p vcf "$out"
done < /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_windows.bed
```


# ARGweaver run ARRAY, started April 18
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARG_run
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%634
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARG_run_%A_%a.out
#SBATCH --error=logs/ARG_run_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

trap 'rmdir "${LOCKDIR}"' EXIT

FINAL_ARG="${OUT_NAME}.${TARGET_ITERS}.smc.gz"
LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"

# Skip completed jobs
if [[ -f "${FINAL_ARG}" ]]; then
    echo "Final output ${FINAL_ARG} already exists. Skipping."
    exit 0
fi

RESUME_FLAG=""
if [[ -f "${LOGFILE}" || -f "${STATSFILE}" ]]; then
    RESUME_FLAG="--resume"
    echo "Previous run detected; will resume."

    if [[ -f "${LOGFILE}" ]]; then
        LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)
        if [[ -n "${LAST_SAMPLE:-}" ]]; then
            echo "Last logged sample = ${LAST_SAMPLE}"
        fi
        if grep -q 'FINISH' "${LOGFILE}"; then
            echo "Log contains FINISH, but ${FINAL_ARG} was not found."
            echo "Proceeding cautiously with --resume."
        fi
    fi
else
    echo "No prior output found; starting fresh."
fi

/home/amonc/ARGweaver/bin/arg-sample \
  --vcf "${VCF_DIR}/${VCF_FILE}" \
  --region "${REGION}" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters "${TARGET_ITERS}" \
  ${RESUME_FLAG} \
  -o "${OUT_NAME}"

# Updated with a lock section

# ARGweaver run ARRAY, started April 18
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARG_run
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%634
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARG_run_%A_%a.out
#SBATCH --error=logs/ARG_run_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

LOCKDIR="${OUTPUT_DIR}/${OUT_NAME}.lock"

if ! mkdir "${LOCKDIR}" 2>/dev/null; then
    echo "Another job is already working on ${OUT_NAME}. Exiting."
    exit 0
fi

trap 'rmdir "${LOCKDIR}"' EXIT

FINAL_ARG="${OUT_NAME}.${TARGET_ITERS}.smc.gz"
LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"

# Skip completed jobs
if [[ -f "${FINAL_ARG}" ]]; then
    echo "Final output ${FINAL_ARG} already exists. Skipping."
    exit 0
fi

RESUME_FLAG=""
if [[ -f "${LOGFILE}" || -f "${STATSFILE}" ]]; then
    RESUME_FLAG="--resume"
    echo "Previous run detected; will resume."

    if [[ -f "${LOGFILE}" ]]; then
        LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)
        if [[ -n "${LAST_SAMPLE:-}" ]]; then
            echo "Last logged sample = ${LAST_SAMPLE}"
        fi
        if grep -q 'FINISH' "${LOGFILE}"; then
            echo "Log contains FINISH, but ${FINAL_ARG} was not found."
            echo "Proceeding cautiously with --resume."
        fi
    fi
else
    echo "No prior output found; starting fresh."
fi

/home/amonc/ARGweaver/bin/arg-sample \
  --vcf "${VCF_DIR}/${VCF_FILE}" \
  --region "${REGION}" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters "${TARGET_ITERS}" \
  ${RESUME_FLAG} \
  -o "${OUT_NAME}"


# To check on runs

# Started runs on April 18, 7PM (be sure to copy over ARGs to ourdisk by May 2)
squeue -u amonc -t RUNNING -h | wc -l # 435 on 4/24 12:40PM

squeue -u amonc -t PENDING -h | wc -l



ls *.0.smc.gz | wc -l

ls *.10.smc.gz | wc -l # 628 as of 4/21 2:05PM

ls *.100.smc.gz | wc -l # 628 as of 4/21 2:04PM

ls *.120.smc.gz | wc -l # 628 as of 4/21 2:07PM

ls *.130.smc.gz | wc -l # 627 as of 4/21 2:07PM

ls *.500.smc.gz | wc -l # 539 as of 4/21 2:03PM

ls *.1000.smc.gz | wc -l # 366 as of 4/21 2:02PM

ls *.2000.smc.gz | wc -l  # 110 as of 4/21 11:55AM, 309 as of 4/23 9:07AM, 493 as of 4/27 8:42AM

ls *.3000.smc.gz | wc -l # 81 as of 4/21 11:55AM, 128 as of 4/23 8:43AM, 183 as of 4/24 12:41PM, 366 by 4/27 8:41AM, 372 by 4/28 8:42AM



# copying over SMC files to ourdisk 

# Dry run
rsync -avn \
  /scratch/amonc/xipho_revision/argweaver/SMC_files/ \
  /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/SMC_files/

# Real run
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=rsync_args_%J_stdout.txt
#SBATCH --error=rsync_args_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=rsync_args
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver

set -euo pipefail

mkdir -p /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/SMC_files

rsync -av \
  /scratch/amonc/xipho_revision/argweaver/SMC_files/ \
  /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/SMC_files/


# Real run, back up VCFs
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=rsync_vcfs_%J_stdout.txt
#SBATCH --error=rsync_vcfs_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=rsync_vcfs
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver

set -euo pipefail

mkdir -p /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/SMC_files

rsync -av \
  /scratch/amonc/xipho_revision/argweaver/vcfs/ \
  /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/vcfs/


# checking details of jobs to see what worked
```
sacct -j 31116691 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ {print $2}' | \
sort | uniq -c
```

# result is 
377 COMPLETED
180 FAILED
77 RUNNING

# list of failed jobs
sacct -j 31116691 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ && $2=="FAILED" {print $1}' > failed_jobs.txt

# command to rework the failed jobs and try to save the runs
I have a file (/scratch/amonc/xipho_revision/argweaver/failed_jobs.txt) that lists the slurm job ID for each failed job in an argweaver job slurm array. Here is a peak inside this file (180 jobs failed):
```
31116691_17
31116691_24
31116691_27
31116691_28
31116691_29
31116691_30
31116691_33
31116691_34
```

I can use these failed job IDs to access the .out file for each of these jobs (in this directory: /scratch/amonc/xipho_revision/argweaver/logs), which have the following format (example below):
```
ARG_run_31116691_9.out
ARG_run_31116691_7.out
ARG_run_31064195_522.out
```

Note that I am specifically interested in the .out files for the 31116691 slurm run, which is the last one I ran.

Inside the .out file, is a short bit of information that allows us to find the OUT_NAME for the associated scaffold, which gives the basename for the associated ARGweaver output files. For example:

```
(base) [amonc@schooner3 logs]$ less ARG_run_31116691_486.out

Task 486
VCF_FILE=scaffold_69_1900000_3900000.vcf.gz
REGION=scaffold_69:1900001-3900000
OUT_NAME=scaffold69-486_out
Previous run detected; will resume.
Last logged sample = 1330
ARG_run_31116691_486.out (END)
```

To fix the failed runs we need to do a few operations on some of 5 argweaver output files for each genomic region. Here are the five outfiles (OUT_NAME and ITERS are placeholders for this example):

OUT_NAME.ITERS.sites.gz
OUT_NAME.ITERS.smc.gz
OUT_NAME.log
OUT_NAME.masked_regions.bed
OUT_NAME.stats


First, we need to place the numerically last iterations of the sites and smc files in a quarantine folder (/scratch/amonc/xipho_revision/argweaver/arg_output_quarantine). Thus, for example, if scaffold99-552_out.990.smc.gz and scaffold99-552_out.990.sites.gz were the last versions (990 is the highest iteration) of the smc and sites files for the scaffold99-552, then those files should go to quarantine.

Next, we need to wipe the contents of the log file from the ITER sample onwards. So for instance, a log file has the following:

```
sample 1429
resample_arg_regions: accept=0.420000
sample time:   9.4 m

prior:      -989444.122630
likelihood: -2676829.915597
joint:      -3666274.038228
nrecombs:   68616
noncompats: 6697
arglen:     12572254070781.550781
max memory: 2718.6 MB

sample 1430
resample_arg_regions: accept=0.455556
sample time:   9.7 m

prior:      -988571.752577
likelihood: -2676828.293660
joint:      -3665400.046237
nrecombs:   68514
noncompats: 6692
arglen:     12540316224119.605469
max memory: 2718.6 MB

sample 1431
resample_arg_regions: accept=0.377143
sample time:   9.2 m

prior:      -989309.448077
likelihood: -2676636.765692
joint:      -3665946.213769
nrecombs:   68581
noncompats: 6695
arglen:     12545724694926.533203
max memory: 2718.6 MB

sample 1432
resample_arg_regions: accept=0.445652
sample time:   9.7 m

prior:      -990493.752005
likelihood: -2676690.042089
joint:      -3667183.794094
nrecombs:   68721
noncompats: 6698
arglen:     12563629991984.548828
max memory: 2718.6 MB
```

So, if the last ITER of the .sites and .smc files is "1430" then the contents of the log file from the line "sample 1430" onwards should be wiped.


Next, the masked_regions.bed file can remain untouched.

Lastly, the .stats file contains content like the following:
```
resample	1428	-988981.830653	-957642.343788	-2677036.159680	-3666017.990332	68556	6710	12568873178510.578125
resample	1429	-989444.122630	-957946.037739	-2676829.915597	-3666274.038228	68616	6697	12572254070781.550781
resample	1430	-988571.752577	-957068.102142	-2676828.293660	-3665400.046237	68514	6692	12540316224119.605469
resample	1431	-989309.448077	-957826.577051	-2676636.765692	-3665946.213769	68581	6695	12545724694926.533203
resample	1432	-990493.752005	-958446.501373	-2676690.042089	-3667183.794094	68721	6698	12563629991984.548828
resample	1433	-990217.527264	-958558.956580	-2676656.867633	-3666874.394897	68672	6701	12567289435631.314453
resample	1434	-990817.705076	-959074.577331	-2676551.734473	-3667369.439550	68714	6710	12579097548349.958984
resample	1435	-990967.953659	-959121.612372	-2676601.559856	-3667569.513514	68734	6715	12582493311928.009766
resample	1436	-991151.255263	-959278.202703	-2676630.414466	-3667781.669728	68748	6716	12581090591140.894531
resample	1437	-990543.550468	-958854.175035	-2676814.680357	-3667358.230825	68695	6726	12572773825944.863281
resample	1438	-990640.876328	-959030.016741	-2676761.029640	-3667401.905968	68697	6721	12573574122866.720703
resample	1439	-990300.135315	-958322.409566	-2676728.686199	-3667028.821514	68679	6715	12563008250997.238281
```

So, if the last ITER (2nd column of the .stats file) of the .sites and .smc files is "1430", then the contents of the .stats file from line "sample 1430" onwards should be wiped. 







# Examples of errors show in different .err files
===== 31116691_486 =====
error: input ARG's sequence names do not match input sequences
*** Error in `/home/amonc/ARGweaver/bin/arg-sample': free(): invalid next size (normal): 0x00000000108d6c00 ***

===== 31116691_491 =====
error: status file is empty
error: resume failed.

===== 31116691_492 =====
error: status file is empty
error: resume failed.

===== 31116691_495 =====
error: input ARG's sequence names do not match input sequences
*** Error in `/home/amonc/ARGweaver/bin/arg-sample': free(): invalid next size (normal): 0x000000000890f2b0 ***

===== 31116691_508 =====
error: bad newick format (line 19727)
error: could not read ARG



# Bash script to fix argweaver fails
```
#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/argweaver
LOGDIR=${WORKDIR}/logs
FAILED=${WORKDIR}/failed_jobs.txt
OUTDIR=${WORKDIR}/SMC_files
QUAR=${WORKDIR}/arg_output_quarantine
ARRAY_ID=31116691

mkdir -p "$QUAR"

while read -r jobid; do
    [[ -z "$jobid" ]] && continue

    # Only process failed jobs from the final array run
    [[ "$jobid" == ${ARRAY_ID}_* ]] || continue

    task=${jobid#${ARRAY_ID}_}
    outlog="${LOGDIR}/ARG_run_${jobid}.out"

    if [[ ! -f "$outlog" ]]; then
        echo "MISSING .out: $outlog"
        continue
    fi

    OUT_NAME=$(awk -F= '/^OUT_NAME=/ {print $2; exit}' "$outlog")

    # Skip if this run was already processed
    if [[ -f "${OUTDIR}/${OUT_NAME}.log.pretruncate.bak" || \
          -f "${OUTDIR}/${OUT_NAME}.stats.pretruncate.bak" ]]; then
        echo "Skipping $OUT_NAME — backup already exists, likely already processed"
        continue
    fi

    if [[ -z "${OUT_NAME:-}" ]]; then
        echo "NO OUT_NAME FOUND: $outlog"
        continue
    fi

    echo
    echo "Processing $jobid  OUT_NAME=$OUT_NAME"

    # Find highest iteration for which both .smc.gz and .sites.gz exist
    last_iter=$(
        find "$OUTDIR" -maxdepth 1 -name "${OUT_NAME}.*.smc.gz" -printf "%f\n" |
        sed -E "s/^${OUT_NAME}\.([0-9]+)\.smc\.gz$/\1/" |
        sort -n |
        while read -r iter; do
            [[ -f "${OUTDIR}/${OUT_NAME}.${iter}.sites.gz" ]] && echo "$iter"
        done |
        tail -n 1
    )

    if [[ -z "${last_iter:-}" ]]; then
        echo "NO matched .smc.gz/.sites.gz pair found for $OUT_NAME"
        continue
    fi

    # Skip completed runs
    if [[ "$last_iter" -eq 3000 ]]; then
        echo "Skipping $OUT_NAME — completed run detected (iteration 3000)"
        continue
    fi

    echo "Last paired iteration = $last_iter"

    smc="${OUTDIR}/${OUT_NAME}.${last_iter}.smc.gz"
    sites="${OUTDIR}/${OUT_NAME}.${last_iter}.sites.gz"
    logfile="${OUTDIR}/${OUT_NAME}.log"
    statsfile="${OUTDIR}/${OUT_NAME}.stats"

    echo "Quarantine:"
    echo "  $smc"
    echo "  $sites"

    if [[ "${RUN:-0}" == "1" ]]; then
        mv -n "$smc" "$QUAR/"
        mv -n "$sites" "$QUAR/"
    fi

    # Truncate log from line: sample ITER
    if [[ -f "$logfile" ]]; then
        echo "Truncate log from: sample $last_iter"
        if [[ "${RUN:-0}" == "1" ]]; then
            cp -n "$logfile" "${logfile}.pretruncate.bak"
            awk -v iter="$last_iter" '
                $1=="sample" && $2==iter {exit}
                {print}
            ' "$logfile" > "${logfile}.tmp"
            mv "${logfile}.tmp" "$logfile"
        fi
    else
        echo "MISSING log: $logfile"
    fi

    # Truncate stats from row where column 2 == ITER
    if [[ -f "$statsfile" ]]; then
        echo "Truncate stats from iteration: $last_iter"
        if [[ "${RUN:-0}" == "1" ]]; then
            cp -n "$statsfile" "${statsfile}.pretruncate.bak"
            awk -v iter="$last_iter" '
                $2==iter {exit}
                {print}
            ' "$statsfile" > "${statsfile}.tmp"
            mv "${statsfile}.tmp" "$statsfile"
        fi
    else
        echo "MISSING stats: $statsfile"
    fi

done < "$FAILED"
```

# Saved as:
fix_failed_argweaver.sh
chmod +x fix_failed_argweaver.sh # make executable


# Run as sbatch
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --job-name=fix_ARG
#SBATCH --output=fix_ARG_%J.out
#SBATCH --error=fix_ARG_%J.err

RUN=1 /scratch/amonc/xipho_revision/argweaver/fix_failed_argweaver.sh

# Ok, now rerunning the ARGweaver script on morning of April 29, 20206
will be interesting to see how many jobs start after doing the cleanup of failed runs

# Again, checking on runs

# Started runs on April 18, 7PM (be sure to copy over ARGs to ourdisk by May 2)
squeue -u amonc -t RUNNING -h | wc -l # 435 on 4/24 12:40PM

squeue -u amonc -t PENDING -h | wc -l

ls *.10.smc.gz | wc -l # 628 as of 4/21 2:05PM

ls *.100.smc.gz | wc -l # 628 as of 4/21 2:04PM

ls *.120.smc.gz | wc -l # 628 as of 4/21 2:07PM

ls *.130.smc.gz | wc -l # 627 as of 4/21 2:07PM

ls *.500.smc.gz | wc -l # 539 as of 4/21 2:03PM

ls *.1000.smc.gz | wc -l # 366 as of 4/21 2:02PM

ls *.2000.smc.gz | wc -l  # 110 as of 4/21 11:55AM, 309 as of 4/23 9:07AM, 493 as of 4/27 8:42AM

ls *.3000.smc.gz | wc -l # 81 as of 4/21 11:55AM, 128 as of 4/23 8:43AM, 183 as of 4/24 12:41PM, 366 by 4/27 8:41AM, 372 by 4/28 8:42AM, 390 by 4/29 4:36



# checking details of jobs to see what worked
```
sacct -j 31140953 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ {print $2}' | \
sort | uniq -c
```

# result
    387 COMPLETED
    176 FAILED
     71 RUNNING


# list of failed jobs
sacct -j 31140953 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ && $2=="FAILED" {print $1}' > failed_jobs_31140953.txt


# Unable to figure out what the source of the name mismatch errors are, so I am aborting those runs and starting again for the failed runs
# command to rework the failed jobs and try to save the runs

I have a file (/scratch/amonc/xipho_revision/argweaver/failed_jobs_31140953.txt) that lists the slurm job ID for each failed job in an argweaver job slurm array. Here is a peak inside this file (176 jobs failed):
```
31140953_17
31140953_24
31140953_28
31140953_29
31140953_30
31140953_33
31140953_34
31140953_35
```

I can use these failed job IDs to access the .out file for each of these jobs (in this directory: /scratch/amonc/xipho_revision/argweaver/logs), which have the following filename format (example below):
```
ARG_run_31140953_17.out
ARG_run_31140953_24.out
ARG_run_31140953_28.out
```

Note that I am specifically interested in the .out files for the 31140953 slurm run, which is the last one I ran.

Inside the .out file, is a short bit of information that allows us to find the OUT_NAME for the associated scaffold, which gives us the basename for the associated ARGweaver output files. For example:

```
(base) [amonc@schooner1 logs]$ less ARG_run_31140953_34.out

Task 34
VCF_FILE=scaffold_2_19000000_21000000.vcf.gz
REGION=scaffold_2:19000001-21000000
OUT_NAME=scaffold2-34_out
Previous run detected; will resume.
Last logged sample = 2979
ARG_run_31140953_34.out (END)
```


Each OUT_NAME is the basename for five different outputs in ARGweaver. I was unable to fix some failed runs so I need to do move the argweaver output files for each failed genomic region so that ARGweaver can start fresh for these regions. All of the argweaver output files are in this folder: /scratch/amonc/xipho_revision/argweaver/SMC_files. Here are the five outfiles (OUT_NAME and ITERS are placeholders for this example):

OUT_NAME.ITERS.sites.gz
OUT_NAME.ITERS.smc.gz
OUT_NAME.log
OUT_NAME.masked_regions.bed
OUT_NAME.stats


For any failed runs that have not reached 3000 iterations, I want to place all of the associated argweaver output files (sites, smc, log, bed, and stats files--identified as failed by their corresponding basenames, though not included runs with 3000 iterations completed) in the quarantine folder (/scratch/amonc/xipho_revision/argweaver/arg_output_quarantine). This will allow argweaver to regenerate all these runs from scratch and to hopefully proceed without errors. This means that, for the sites and smc files, all the saved iterations (not just the last one) need to go to the quarantine folder.


# Bash command to clean out SMC_files folder of failed runs
#!/bin/bash
set -euo pipefail
shopt -s nullglob

FAILED_LIST="/scratch/amonc/xipho_revision/argweaver/failed_jobs_31140953.txt"
LOG_DIR="/scratch/amonc/xipho_revision/argweaver/logs"
SMC_DIR="/scratch/amonc/xipho_revision/argweaver/SMC_files"
QUAR_DIR="/scratch/amonc/xipho_revision/argweaver/arg_output_quarantine"

mkdir -p "$QUAR_DIR"

# Set to 1 after checking the dry run output
DO_MOVE=0

while read -r failed_id; do
    [[ -z "$failed_id" ]] && continue

    out_file="${LOG_DIR}/ARG_run_${failed_id}.out"

    if [[ ! -f "$out_file" ]]; then
        echo "WARNING: Missing out file: $out_file"
        continue
    fi

    OUT_NAME=$(grep -m 1 '^OUT_NAME=' "$out_file" | cut -d'=' -f2)

    if [[ -z "${OUT_NAME:-}" ]]; then
        echo "WARNING: Could not find OUT_NAME in $out_file"
        continue
    fi

    # Skip complete runs safely
    COMPLETE=0

    if [[ -f "${SMC_DIR}/${OUT_NAME}.3000.smc.gz" ]] || \
       [[ -f "${SMC_DIR}/${OUT_NAME}.3000.sites.gz" ]]; then
        COMPLETE=1
    elif [[ -f "${SMC_DIR}/${OUT_NAME}.log" ]]; then
        if grep -qE 'Last logged sample = 3000|sample[[:space:]]+3000' "${SMC_DIR}/${OUT_NAME}.log"; then
            COMPLETE=1
        fi
    fi

    if (( COMPLETE == 1 )); then
        echo "SKIP complete run: $OUT_NAME"
        continue
    fi

    files=()

    for pattern in \
        "${SMC_DIR}/${OUT_NAME}".*.sites.gz \
        "${SMC_DIR}/${OUT_NAME}".*.smc.gz \
        "${SMC_DIR}/${OUT_NAME}.log" \
        "${SMC_DIR}/${OUT_NAME}.masked_regions.bed" \
        "${SMC_DIR}/${OUT_NAME}.stats"
    do
        for f in $pattern; do
            [[ -e "$f" ]] && files+=("$f")
        done
    done

    if (( ${#files[@]} == 0 )); then
        echo "WARNING: No files found for $OUT_NAME"
        continue
    fi

    if (( ${#files[@]} == 0 )); then
        echo "WARNING: No files found for $OUT_NAME"
        continue
    fi

    echo "Quarantine candidate: $OUT_NAME"
    printf '  %s\n' "${files[@]}"

    if (( DO_MOVE == 1 )); then
        mv -v "${files[@]}" "$QUAR_DIR"/
    fi

done < "$FAILED_LIST"

echo "Done."


# Run as sbatch
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=arg_quarantine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/quarantine_%J.out
#SBATCH --error=logs/quarantine_%J.err
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver

set -euo pipefail

bash quarantine_failed_argweaver_outputs.sh


# Started runs on April 18, 7PM (be sure to copy over ARGs to ourdisk by May 2)
squeue -u amonc -t RUNNING -h | wc -l # 435 on 4/24 12:40PM

squeue -u amonc -t PENDING -h | wc -l

ls *.10.smc.gz | wc -l # 628 as of 4/21 2:05PM

ls *.100.smc.gz | wc -l # 628 as of 4/21 2:04PM

ls *.120.smc.gz | wc -l # 628 as of 4/21 2:07PM

ls *.130.smc.gz | wc -l # 627 as of 4/21 2:07PM

ls *.500.smc.gz | wc -l # 539 as of 4/21 2:03PM

ls *.1000.smc.gz | wc -l # 366 as of 4/21 2:02PM

ls *.2000.smc.gz | wc -l  # 110 as of 4/21 11:55AM, 309 as of 4/23 9:07AM, 493 as of 4/27 8:42AM

ls *.3000.smc.gz | wc -l # 81 as of 4/21 11:55AM, 128 as of 4/23 8:43AM, 183 as of 4/24 12:41PM, 366 by 4/27 8:41AM, 372 by 4/28 8:42AM, 390 by 4/29 4:36, 408 by 5/2 7:51AM, 414 by 5/4 12:55pm



# Checking to see what worked or failed
sacct -j 31291355 --format=JobID,State -n | awk '$1 ~ /_[0-9]+$/ {print $2}' | sort | uniq -c
    441 COMPLETED
      4 FAILED
    189 RUNNING


# discrepancy between 414 and 441, why?!

# Get list of tasks completed
sacct -j 31291355 --format=JobID,State -n \
| awk '$1 ~ /_[0-9]+$/ && $2=="COMPLETED" {split($1,a,"_"); print a[2]}' \
| sort -n > completed_441.txt

# Appears to be an issue with some files being locked! Interesting! That may have saved me from messing up the files.
# Clearing files with locks more than 12 hours old
find /scratch/amonc/xipho_revision/argweaver/SMC_files \
  -maxdepth 1 -name "*.lock" -type d -mmin +720 -exec rmdir {} \;


  # Checking to see what worked or failed
sacct -j 31366120 --format=JobID,State -n | awk '$1 ~ /_[0-9]+$/ {print $2}' | sort | uniq -c




# 7 May 2026; list of failed jobs
sacct -j 31366120 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ && $2=="FAILED" {print $1}' > failed_jobs_31366120.txt

# Failed window
scaffold_2_7600000_9600000.vcf.gz
scaffold2-28_out.log
scaffold2-28_out.stats
scaffold2-28_out.2170.sites.gz
scaffold2-28_out.2170.smc.gz
ARG_run_31366120_28.err

# Successful window on resume
scaffold2-31_out.1100.sites.gz
scaffold2-31_out.1100.smc.gz
scaffold_2_13300000_15300000.vcf.gz



# 7 May 2026
# I am having persistent issues with an error while running ARGweaver. Here is the error message:
```
(base) [amonc@schooner1 logs]$ less ARG_run_31366120_28.err

masked 14627 (3.7%) positions

model: 
  mu = 2.300000e-08
  rho = 7.500000e-08
  smc_prime = false
  ntimes = 20
  times = [0.000000,153.463573,424.682488,904.012522,1751.141055,3248.286447,5894.218246,10570.420787,18834.757085,33440.466280,59253.397090,104873.050894,185497.480121,327986.474363,579809.819649,1024861.716948,1811409.886515,3201490.519369,5658204.857537,10000000.000000]
  npop = 1
  popsizes = [391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0]

error: input ARG's sequence names do not match input sequences
```

I've attached the VCF file (scaffold_2_7600000_9600000.vcf.gz; which has the original ARG sequence names) and the ARGweaver-generated sites.gz file (scaffold2-28_out.2170.sites.gz), which I believe are read on resuming ARGweaver. First off, the VCF file does not have the names split into _1 or _2 haplotypes. I don't think that is the issue, because I have no problem resuming runs for other windows using the same script.


For comparison, I have attached a VCF (scaffold_2_13300000_15300000.vcf.gz)  and sites.gz file (scaffold2-31_out.1100.sites.gz) that worked no problem when resuming:
```
RESUME
arg-sample 1.0
start time: Tue Apr 21 08:37:04 2026
command: /home/amonc/ARGweaver/bin/arg-sample --vcf /scratch/amonc/xipho_revision/argweaver/vcfs/scaffold_2_13300000_15300000.vcf.gz --region scaffold_2:13300001-15300000 --vcf-genotype-filter DP<5;DP>50;GQ<20;RGQ<20 --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed --mask-cluster 2,5 --mask-Ns 33 --mutrate 4.6e-9 --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed --popsize 391476 --compress-seq 5 --ntimes 20 --maxtime 1e7 --delta 0.005 --sample-step 10 --iters 3000 --resume -o scaffold2-31_out
random seed: 1776778624
read input sites (chrom=scaffold_2, start=13300000, end=15300000, length=2000000, nseqs=66, nsites=355042)
Reading scaffold2-31_out.masked_regions.bed
Removed 0 sites overlapping mask (old=355042 , new=355042)
327292 sites are partially masked but otherwise invariant
masked 15147 (3.8%) positions

model: 
  mu = 2.300000e-08
  rho = 7.500000e-08
  smc_prime = false
  ntimes = 20
  times = [0.000000,153.463573,424.682488,904.012522,1751.141055,3248.286447,5894.218246,10570.420787,18834.757085,33440.466280,59253.397090,104873.050894,185497.480121,327986.474363,579809.819649,1024861.716948,1811409.886515,3201490.519369,5658204.857537,10000000.000000]
  npop = 1
  popsizes = [391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0,
              391476.0]

read input ARG (chrom=scaffold_2, start=13300000, end=15300000, nseqs=66)
max memory usage: 231.3 MB

Resample All Branches (3000 iterations)
--------------------------------------
sample 1101
resample_arg_regions: accept=0.854545
sample time:   4.8 m

prior:      -199168.229198
likelihood: -2969547.817035
joint:      -3168716.046233
nrecombs:   11334
noncompats: 3613
arglen:     8065901413178.759766
max memory: 374.1 MB

sample 1102
resample_arg_regions: accept=0.882927
sample time:   5.4 m

prior:      -199089.058332
likelihood: -2969684.253408
joint:      -3168773.311740
nrecombs:   11340
noncompats: 3622
arglen:     8048662045352.815430
max memory: 374.1 MB
```

Can you detect any reason why the first pair for VCF + sites files failed while the second pair succeeded?


# list of failed jobs
sacct -j 31366120 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ && $2=="FAILED" {print $1}' > failed_jobs_31366120.txt

# PLAN 
Ok, I think I have a plan. I need a script that first explores ARGweaver output (/scratch/amonc/xipho_revision/argweaver/SMC_files) to see if there are any completed genomic regions (with 3000 iterations). These regions will be identified by the presence of files ending with ".3000.smc.gz". These regions are then ignored.

For regions with fewer than 3000 iterations, I need to then inspect the first line of the .smc.gz files, which contain the sample names and is the source of the errors I've been encountering. If the name for a given individual is only a number (e.g., "64" or "51"), then that will certainly generate the name mismatch error. However, all alphanumeric names (like "Bel_xsTP32151_Toe_2") signify a correct smc.gz file. I need the script to then identify the largest iteration number of smc.gz file (for a given region) that has only alphanumeric names in the first line (underscores are fine).

Once I have identified the largest iteration of smc.gz file that has correct sample names (the last correct smc.gz file), I will need to do a few operations on some of the 5 argweaver output files for each genomic region. Here are the five outfiles (OUT_NAME and ITERS are placeholders for this example):

OUT_NAME.ITERS.sites.gz
OUT_NAME.ITERS.smc.gz
OUT_NAME.log
OUT_NAME.masked_regions.bed
OUT_NAME.stats


First, we need to place all sites and smc files with iterations higher than the last correct smc.gz file in a quarantine folder (/scratch/amonc/xipho_revision/argweaver/arg_output_quarantine). Thus, for example, if scaffold99-552_out.990.smc.gz were the last correct iteration, then scaffold99-552_out.1000.smc.gz and scaffold99-552_out.1010.smc.gz etc. should go to quarantine.

Next, we need to partially wipe the contents of the log file--specifically, we will wipe lines for ITERS above the last correct iteration. So for instance, a log file has the following:

```
sample 1429
resample_arg_regions: accept=0.420000
sample time:   9.4 m

prior:      -989444.122630
likelihood: -2676829.915597
joint:      -3666274.038228
nrecombs:   68616
noncompats: 6697
arglen:     12572254070781.550781
max memory: 2718.6 MB

sample 1430
resample_arg_regions: accept=0.455556
sample time:   9.7 m

prior:      -988571.752577
likelihood: -2676828.293660
joint:      -3665400.046237
nrecombs:   68514
noncompats: 6692
arglen:     12540316224119.605469
max memory: 2718.6 MB

sample 1431
resample_arg_regions: accept=0.377143
sample time:   9.2 m

prior:      -989309.448077
likelihood: -2676636.765692
joint:      -3665946.213769
nrecombs:   68581
noncompats: 6695
arglen:     12545724694926.533203
max memory: 2718.6 MB

sample 1432
resample_arg_regions: accept=0.445652
sample time:   9.7 m

prior:      -990493.752005
likelihood: -2676690.042089
joint:      -3667183.794094
nrecombs:   68721
noncompats: 6698
arglen:     12563629991984.548828
max memory: 2718.6 MB
```

So, if the last correct ITER of the .smc files is "1430", then the contents of the log file from the line "sample 1431" onwards should be wiped.


Next, the masked_regions.bed file can remain untouched.

Lastly, the .stats file contains content like the following:
```
resample	1428	-988981.830653	-957642.343788	-2677036.159680	-3666017.990332	68556	6710	12568873178510.578125
resample	1429	-989444.122630	-957946.037739	-2676829.915597	-3666274.038228	68616	6697	12572254070781.550781
resample	1430	-988571.752577	-957068.102142	-2676828.293660	-3665400.046237	68514	6692	12540316224119.605469
resample	1431	-989309.448077	-957826.577051	-2676636.765692	-3665946.213769	68581	6695	12545724694926.533203
resample	1432	-990493.752005	-958446.501373	-2676690.042089	-3667183.794094	68721	6698	12563629991984.548828
resample	1433	-990217.527264	-958558.956580	-2676656.867633	-3666874.394897	68672	6701	12567289435631.314453
resample	1434	-990817.705076	-959074.577331	-2676551.734473	-3667369.439550	68714	6710	12579097548349.958984
resample	1435	-990967.953659	-959121.612372	-2676601.559856	-3667569.513514	68734	6715	12582493311928.009766
resample	1436	-991151.255263	-959278.202703	-2676630.414466	-3667781.669728	68748	6716	12581090591140.894531
resample	1437	-990543.550468	-958854.175035	-2676814.680357	-3667358.230825	68695	6726	12572773825944.863281
resample	1438	-990640.876328	-959030.016741	-2676761.029640	-3667401.905968	68697	6721	12573574122866.720703
resample	1439	-990300.135315	-958322.409566	-2676728.686199	-3667028.821514	68679	6715	12563008250997.238281
```

So, if the last correct ITER of a given .smc.gz file is "1430", then the contents of the .stats file from line "resample 1431" onwards should be wiped. 

There may be some regions listed in my jobs file (see below; /scratch/amonc/xipho_revision/argweaver/jobs_file.txt) that have some ARGweaver output (like .bed or .log files) but no .smc.gz files. In those cases, these outputs can be moved to the quarantine folder (/scratch/amonc/xipho_revision/argweaver/arg_output_quarantine) prior to beginning a fresh run of this genomic region.

Ok, now this above functionality should be added to my ARG_run sbatch script below:

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARG_run
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%634
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARG_run_%A_%a.out
#SBATCH --error=logs/ARG_run_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

LOCKDIR="${OUTPUT_DIR}/${OUT_NAME}.lock"

if ! mkdir "${LOCKDIR}" 2>/dev/null; then
    echo "Another job is already working on ${OUT_NAME}. Exiting."
    exit 0
fi

trap 'rmdir "${LOCKDIR}"' EXIT

FINAL_ARG="${OUT_NAME}.${TARGET_ITERS}.smc.gz"
LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"

# Skip completed jobs
if [[ -f "${FINAL_ARG}" ]]; then
    echo "Final output ${FINAL_ARG} already exists. Skipping."
    exit 0
fi

RESUME_FLAG=""
if [[ -f "${LOGFILE}" || -f "${STATSFILE}" ]]; then
    RESUME_FLAG="--resume"
    echo "Previous run detected; will resume."

    if [[ -f "${LOGFILE}" ]]; then
        LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)
        if [[ -n "${LAST_SAMPLE:-}" ]]; then
            echo "Last logged sample = ${LAST_SAMPLE}"
        fi
        if grep -q 'FINISH' "${LOGFILE}"; then
            echo "Log contains FINISH, but ${FINAL_ARG} was not found."
            echo "Proceeding cautiously with --resume."
        fi
    fi
else
    echo "No prior output found; starting fresh."
fi

/home/amonc/ARGweaver/bin/arg-sample \
  --vcf "${VCF_DIR}/${VCF_FILE}" \
  --region "${REGION}" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters "${TARGET_ITERS}" \
  ${RESUME_FLAG} \
  -o "${OUT_NAME}"




# NEW AND MORE ROBUST SCRIPT. ARG_run_v3.sbatch
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARG_runv3
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%634
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARG_runv3_%A_%a.out
#SBATCH --error=logs/ARG_runv3_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
QUARANTINE_DIR=${WORKDIR}/arg_output_quarantine
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${QUARANTINE_DIR}"

# ----------------------------
# Helper function:
# Check whether the first line of an SMC file has valid sequence names.
# Valid names must contain at least one letter or underscore, and only:
# letters, numbers, and underscores.
# Pure numbers like 62, 63, 64, 65 are treated as corrupt.
# ----------------------------
smc_names_valid() {
    local smc="$1"
    local names_line
    local name

    names_line=$(zcat "$smc" 2>/dev/null | head -n 1 || true)

    if [[ ! "$names_line" =~ ^NAMES[[:space:]] ]]; then
        return 1
    fi

    for name in ${names_line#NAMES }; do

        # Purely numeric names are corrupt
        if [[ "$name" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        # Weird characters are suspicious/corrupt
        if [[ ! "$name" =~ ^[A-Za-z0-9_-]+$ ]]; then
            return 1
        fi
    done

    return 0
}

# ----------------------------
# Helper function:
# Extract iteration number from files like:
# OUT_NAME.1430.smc.gz
# OUT_NAME.1430.sites.gz
# ----------------------------
get_iter_from_arg_file() {
    local file="$1"
    basename "$file" | sed -E 's/.*\.([0-9]+)\.(smc|sites)\.gz/\1/'
}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

LOCKDIR="${OUTPUT_DIR}/${OUT_NAME}.lock"

if ! mkdir "${LOCKDIR}" 2>/dev/null; then
    echo "Another job is already working on ${OUT_NAME}. Exiting."
    exit 0
fi

trap 'rmdir "${LOCKDIR}"' EXIT

FINAL_ARG="${OUT_NAME}.${TARGET_ITERS}.smc.gz"
LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"

# ----------------------------
# Skip completed jobs only if the final SMC exists and has valid names.
# ----------------------------
if [[ -f "${FINAL_ARG}" ]]; then
    if smc_names_valid "${FINAL_ARG}"; then
        echo "Final output ${FINAL_ARG} already exists and appears valid. Skipping."
        exit 0
    else
        echo "Final output ${FINAL_ARG} exists but has invalid names. Will quarantine/recover."
    fi
fi

# ----------------------------
# Inspect existing SMC checkpoints.
# Find the largest valid iteration.
# ----------------------------
shopt -s nullglob

SMC_FILES=( ${OUT_NAME}.*.smc.gz )
LAST_VALID_ITER=""

if (( ${#SMC_FILES[@]} > 0 )); then

    echo "Inspecting existing SMC checkpoints for ${OUT_NAME}..."

    for smc in "${SMC_FILES[@]}"; do
        iter=$(get_iter_from_arg_file "$smc")

        if smc_names_valid "$smc"; then
            echo "VALID   ${smc}"

            if [[ -z "${LAST_VALID_ITER}" || "$iter" -gt "$LAST_VALID_ITER" ]]; then
                LAST_VALID_ITER="$iter"
            fi
        else
            echo "CORRUPT ${smc}"
        fi
    done

    echo "Last valid SMC iteration = ${LAST_VALID_ITER:-NONE}"

else
    echo "No SMC checkpoints found for ${OUT_NAME}."
fi

# ----------------------------
# If there are partial outputs but no valid SMC checkpoints,
# quarantine all outputs for this region and start fresh.
# ----------------------------
if [[ -z "${LAST_VALID_ITER}" ]]; then

    PARTIAL_OUTPUTS=( ${OUT_NAME}* )

    if (( ${#PARTIAL_OUTPUTS[@]} > 0 )); then
        echo "Partial outputs exist but no valid SMC checkpoint was found."
        echo "Quarantining all outputs for ${OUT_NAME}."

        for f in "${PARTIAL_OUTPUTS[@]}"; do
            echo "Moving $f to ${QUARANTINE_DIR}/"
            mv "$f" "${QUARANTINE_DIR}/"
        done
    fi

    RESUME_FLAG=""
    echo "Starting fresh."

else

    # ----------------------------
    # Quarantine all SMC/sites files above the last valid iteration.
    # ----------------------------
    echo "Quarantining SMC/sites files with iterations > ${LAST_VALID_ITER}."

    for f in ${OUT_NAME}.*.smc.gz ${OUT_NAME}.*.sites.gz; do
        [[ -e "$f" ]] || continue

        iter=$(get_iter_from_arg_file "$f")

        if [[ "$iter" -gt "$LAST_VALID_ITER" ]]; then
            echo "Moving $f to ${QUARANTINE_DIR}/"
            mv "$f" "${QUARANTINE_DIR}/"
        fi
    done

    # ----------------------------
    # Truncate log file after the last valid sample.
    # ----------------------------
    if [[ -f "${LOGFILE}" ]]; then
        echo "Truncating ${LOGFILE} after sample ${LAST_VALID_ITER}."

        TMP_LOG="${LOGFILE}.tmp"

        awk -v maxiter="${LAST_VALID_ITER}" '
            /^sample [0-9]+$/ {
                split($0,a," ")
                if (a[2] > maxiter) exit
            }
            { print }
        ' "${LOGFILE}" > "${TMP_LOG}"

        mv "${TMP_LOG}" "${LOGFILE}"
    fi

    # ----------------------------
    # Truncate stats file after the last valid iteration.
    # Iteration number is column 2.
    # ----------------------------
    if [[ -f "${STATSFILE}" ]]; then
        echo "Truncating ${STATSFILE} after iteration ${LAST_VALID_ITER}."

        TMP_STATS="${STATSFILE}.tmp"

        awk -v maxiter="${LAST_VALID_ITER}" '
            $2 > maxiter { exit }
            { print }
        ' "${STATSFILE}" > "${TMP_STATS}"

        mv "${TMP_STATS}" "${STATSFILE}"
    fi

    RESUME_FLAG="--resume"
    echo "Previous valid run detected; will resume from iteration ${LAST_VALID_ITER}."
fi

# ----------------------------
# Extra log diagnostics
# ----------------------------
if [[ -f "${LOGFILE}" ]]; then
    LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)

    if [[ -n "${LAST_SAMPLE:-}" ]]; then
        echo "Last logged sample after cleanup = ${LAST_SAMPLE}"
    fi

    if grep -q 'FINISH' "${LOGFILE}" && [[ ! -f "${FINAL_ARG}" ]]; then
        echo "Log contains FINISH, but ${FINAL_ARG} was not found."
        echo "Proceeding cautiously."
    fi
fi

# ----------------------------
# Run ARGweaver
# ----------------------------
/home/amonc/ARGweaver/bin/arg-sample \
  --vcf "${VCF_DIR}/${VCF_FILE}" \
  --region "${REGION}" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters "${TARGET_ITERS}" \
  ${RESUME_FLAG} \
  -o "${OUT_NAME}"

# Started runs on April 18, 7PM (be sure to copy over ARGs to ourdisk by May 2)
squeue -u amonc -t RUNNING -h | wc -l # 435 on 4/24 12:40PM

squeue -u amonc -t PENDING -h | wc -l

ls *.10.smc.gz | wc -l

ls *.10.smc.gz | wc -l # 628 as of 4/21 2:05PM

ls *.100.smc.gz | wc -l # 628 as of 4/21 2:04PM

ls *.120.smc.gz | wc -l # 628 as of 4/21 2:07PM

ls *.130.smc.gz | wc -l # 627 as of 4/21 2:07PM

ls *.500.smc.gz | wc -l # 539 as of 4/21 2:03PM

ls *.1000.smc.gz | wc -l # 366 as of 4/21 2:02PM

ls *.2000.smc.gz | wc -l  # 110 as of 4/21 11:55AM, 309 as of 4/23 9:07AM, 493 as of 4/27 8:42AM

ls *.3000.smc.gz | wc -l # 81 as of 4/21 11:55AM, 128 as of 4/23 8:43AM, 183 as of 4/24 12:41PM, 366 by 4/27 8:41AM, 372 by 4/28 8:42AM, 390 by 4/29 4:36, 408 by 5/2 7:51AM, 414 by 5/4 12:55pm; 455 by 5/8 8:07am



# checking details of jobs to see what worked
```
sacct -j 31404586 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ {print $2}' | \
sort | uniq -c
```
      3 FAILED
    631 RUNNING

# Ok, why are 631 still running?!

ls *-17_out.* -lht
scaffold2-31_out.1100.smc.gz

# Nice little script to find the first corrupted smc file for a given out_name
```
./find_first_corrupt_smc.sh scaffold37-357_out


./find_first_corrupt_smc.sh scaffold1-17_out
```

# New script 
I would like to update this script so that it instead finds the first corrupted smc.gz file (first smc.gz file with an incorrect sample name in the first line). Then the script needs to move this smc file and all other sites and smc files with a greater number of iterations to the arg quarantine folder. Similarly the log and stats files need to be updated so that they don't include information from the corrupted smc.gz iteration onwards.

Here is the script that needs updating (please provide a full script replacement):

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARG_runv3
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%634
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARG_runv3_%A_%a.out
#SBATCH --error=logs/ARG_runv3_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
QUARANTINE_DIR=${WORKDIR}/arg_output_quarantine
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${QUARANTINE_DIR}"

# ----------------------------
# Helper function:
# Check whether the first line of an SMC file has valid sequence names.
# Valid names must contain at least one letter or underscore, and only:
# letters, numbers, and underscores.
# Pure numbers like 62, 63, 64, 65 are treated as corrupt.
# ----------------------------
smc_names_valid() {
    local smc="$1"
    local names_line
    local name

    names_line=$(zcat "$smc" 2>/dev/null | head -n 1 || true)

    if [[ ! "$names_line" =~ ^NAMES[[:space:]] ]]; then
        return 1
    fi

    for name in ${names_line#NAMES }; do

        # Purely numeric names are corrupt
        if [[ "$name" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        # Weird characters are suspicious/corrupt
        if [[ ! "$name" =~ ^[A-Za-z0-9_-]+$ ]]; then
            return 1
        fi
    done

    return 0
}

# ----------------------------
# Helper function:
# Extract iteration number from files like:
# OUT_NAME.1430.smc.gz
# OUT_NAME.1430.sites.gz
# ----------------------------
get_iter_from_arg_file() {
    local file="$1"
    basename "$file" | sed -E 's/.*\.([0-9]+)\.(smc|sites)\.gz/\1/'
}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

LOCKDIR="${OUTPUT_DIR}/${OUT_NAME}.lock"

if ! mkdir "${LOCKDIR}" 2>/dev/null; then
    echo "Another job is already working on ${OUT_NAME}. Exiting."
    exit 0
fi

trap 'rmdir "${LOCKDIR}"' EXIT

FINAL_ARG="${OUT_NAME}.${TARGET_ITERS}.smc.gz"
LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"

# ----------------------------
# Skip completed jobs only if the final SMC exists and has valid names.
# ----------------------------
if [[ -f "${FINAL_ARG}" ]]; then
    if smc_names_valid "${FINAL_ARG}"; then
        echo "Final output ${FINAL_ARG} already exists and appears valid. Skipping."
        exit 0
    else
        echo "Final output ${FINAL_ARG} exists but has invalid names. Will quarantine/recover."
    fi
fi

# ----------------------------
# Inspect existing SMC checkpoints.
# Find the largest valid iteration.
# ----------------------------
shopt -s nullglob

SMC_FILES=( ${OUT_NAME}.*.smc.gz )
LAST_VALID_ITER=""

if (( ${#SMC_FILES[@]} > 0 )); then

    echo "Inspecting existing SMC checkpoints for ${OUT_NAME}..."

    for smc in "${SMC_FILES[@]}"; do
        iter=$(get_iter_from_arg_file "$smc")

        if smc_names_valid "$smc"; then
            echo "VALID   ${smc}"

            if [[ -z "${LAST_VALID_ITER}" || "$iter" -gt "$LAST_VALID_ITER" ]]; then
                LAST_VALID_ITER="$iter"
            fi
        else
            echo "CORRUPT ${smc}"
        fi
    done

    echo "Last valid SMC iteration = ${LAST_VALID_ITER:-NONE}"

else
    echo "No SMC checkpoints found for ${OUT_NAME}."
fi

# ----------------------------
# If there are partial outputs but no valid SMC checkpoints,
# quarantine all outputs for this region and start fresh.
# ----------------------------
if [[ -z "${LAST_VALID_ITER}" ]]; then

    PARTIAL_OUTPUTS=( ${OUT_NAME}* )

    if (( ${#PARTIAL_OUTPUTS[@]} > 0 )); then
        echo "Partial outputs exist but no valid SMC checkpoint was found."
        echo "Quarantining all outputs for ${OUT_NAME}."

        for f in "${PARTIAL_OUTPUTS[@]}"; do
            echo "Moving $f to ${QUARANTINE_DIR}/"
            mv "$f" "${QUARANTINE_DIR}/"
        done
    fi

    RESUME_FLAG=""
    echo "Starting fresh."

else

    # ----------------------------
    # Quarantine all SMC/sites files above the last valid iteration.
    # ----------------------------
    echo "Quarantining SMC/sites files with iterations > ${LAST_VALID_ITER}."

    for f in ${OUT_NAME}.*.smc.gz ${OUT_NAME}.*.sites.gz; do
        [[ -e "$f" ]] || continue

        iter=$(get_iter_from_arg_file "$f")

        if [[ "$iter" -gt "$LAST_VALID_ITER" ]]; then
            echo "Moving $f to ${QUARANTINE_DIR}/"
            mv "$f" "${QUARANTINE_DIR}/"
        fi
    done

    # ----------------------------
    # Truncate log file after the last valid sample.
    # ----------------------------
    if [[ -f "${LOGFILE}" ]]; then
        echo "Truncating ${LOGFILE} after sample ${LAST_VALID_ITER}."

        TMP_LOG="${LOGFILE}.tmp"

        awk -v maxiter="${LAST_VALID_ITER}" '
            /^sample [0-9]+$/ {
                split($0,a," ")
                if (a[2] > maxiter) exit
            }
            { print }
        ' "${LOGFILE}" > "${TMP_LOG}"

        mv "${TMP_LOG}" "${LOGFILE}"
    fi

    # ----------------------------
    # Truncate stats file after the last valid iteration.
    # Iteration number is column 2.
    # ----------------------------
    if [[ -f "${STATSFILE}" ]]; then
        echo "Truncating ${STATSFILE} after iteration ${LAST_VALID_ITER}."

        TMP_STATS="${STATSFILE}.tmp"

        awk -v maxiter="${LAST_VALID_ITER}" '
            $2 > maxiter { exit }
            { print }
        ' "${STATSFILE}" > "${TMP_STATS}"

        mv "${TMP_STATS}" "${STATSFILE}"
    fi

    RESUME_FLAG="--resume"
    echo "Previous valid run detected; will resume from iteration ${LAST_VALID_ITER}."
fi

# ----------------------------
# Extra log diagnostics
# ----------------------------
if [[ -f "${LOGFILE}" ]]; then
    LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)

    if [[ -n "${LAST_SAMPLE:-}" ]]; then
        echo "Last logged sample after cleanup = ${LAST_SAMPLE}"
    fi

    if grep -q 'FINISH' "${LOGFILE}" && [[ ! -f "${FINAL_ARG}" ]]; then
        echo "Log contains FINISH, but ${FINAL_ARG} was not found."
        echo "Proceeding cautiously."
    fi
fi

# ----------------------------
# Run ARGweaver
# ----------------------------
/home/amonc/ARGweaver/bin/arg-sample \
  --vcf "${VCF_DIR}/${VCF_FILE}" \
  --region "${REGION}" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters "${TARGET_ITERS}" \
  ${RESUME_FLAG} \
  -o "${OUT_NAME}"




---------------------------
# Beautiful new script
# Lock logic removed
# quarantining step improved--actually deleting now
# Verification of non-corrupted ITERS saved







#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARGsv11
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --array=1-634%300
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ARGsv11_%A_%a.out
#SBATCH --error=logs/ARGsv11_%A_%a.err

set -euo pipefail

module load GCCcore/11.3.0

JOBFILE=/scratch/amonc/xipho_revision/argweaver/jobs_file.txt
VCF_DIR=/scratch/amonc/xipho_revision/argweaver/vcfs
WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
VERIFY_DIR=${WORKDIR}/verified_iterations
TARGET_ITERS=3000

mkdir -p "${WORKDIR}/logs"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${VERIFY_DIR}"

smc_names_valid() {
    local smc="$1"
    local names_line
    local name

    names_line=$(zcat "$smc" 2>/dev/null | head -n 1 || true)

    if [[ ! "$names_line" =~ ^NAMES[[:space:]] ]]; then
        return 1
    fi

    for name in ${names_line#NAMES }; do
        if [[ "$name" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        if [[ ! "$name" =~ ^[A-Za-z0-9_-]+$ ]]; then
            return 1
        fi
    done

    return 0
}

get_iter_from_arg_file() {
    local file="$1"

    basename "$file" | sed -E 's/.*\.([0-9]+)\.(smc|sites)\.gz/\1/'
}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBFILE")

VCF_FILE=$(printf '%s\n' "$LINE" | cut -f1)
REGION=$(printf '%s\n' "$LINE" | cut -f2)
OUT_NAME=$(printf '%s\n' "$LINE" | cut -f3)

if [[ -z "${VCF_FILE}" || -z "${REGION}" || -z "${OUT_NAME}" ]]; then
    echo "Could not parse line ${SLURM_ARRAY_TASK_ID} from ${JOBFILE}"
    exit 1
fi

if [[ ! -f "${VCF_DIR}/${VCF_FILE}" ]]; then
    echo "Missing VCF: ${VCF_DIR}/${VCF_FILE}"
    exit 1
fi

cd "${OUTPUT_DIR}"

LOGFILE="${OUT_NAME}.log"
STATSFILE="${OUT_NAME}.stats"
VERIFY_FILE="${VERIFY_DIR}/${OUT_NAME}.verified"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VCF_FILE=${VCF_FILE}"
echo "REGION=${REGION}"
echo "OUT_NAME=${OUT_NAME}"
echo "VERIFY_FILE=${VERIFY_FILE}"

START_ITER=-1

if [[ -f "${VERIFY_FILE}" ]]; then
    START_ITER=$(cat "${VERIFY_FILE}")

    if [[ ! "${START_ITER}" =~ ^[0-9]+$ ]]; then
        echo "WARNING: ${VERIFY_FILE} does not contain a valid integer."
        echo "Ignoring verification checkpoint."
        START_ITER=-1
    fi
fi

echo "Previously verified through iteration: ${START_ITER}"

shopt -s nullglob

SMC_FILES=( ${OUT_NAME}.*.smc.gz )

FIRST_CORRUPT_ITER=""
LAST_VALID_ITER="${START_ITER}"
CHECKED_COUNT=0
SKIPPED_COUNT=0

if (( ${#SMC_FILES[@]} > 0 )); then

    echo "Inspecting SMC checkpoints for ${OUT_NAME}"

    while read -r smc; do
        iter=$(get_iter_from_arg_file "$smc")

        if [[ "$iter" -le "$START_ITER" ]]; then
            ((++SKIPPED_COUNT))
            continue
        fi

        ((++CHECKED_COUNT))

        if smc_names_valid "$smc"; then
            echo "VALID   ${smc}"

            LAST_VALID_ITER="$iter"
            echo "${iter}" > "${VERIFY_FILE}"

        else
            echo "CORRUPT ${smc}"

            FIRST_CORRUPT_ITER="$iter"
            break
        fi

    done < <(
        printf '%s\n' "${SMC_FILES[@]}" |
        sort -t. -k2,2n
    )

else
    echo "No SMC files found."
fi

echo
echo "Skipped previously verified files: ${SKIPPED_COUNT}"
echo "Checked files this run: ${CHECKED_COUNT}"
echo "LAST_VALID_ITER=${LAST_VALID_ITER:-NONE}"
echo "FIRST_CORRUPT_ITER=${FIRST_CORRUPT_ITER:-NONE}"
echo

if [[ -z "${FIRST_CORRUPT_ITER}" ]]; then

    if [[ -n "${LAST_VALID_ITER}" && "${LAST_VALID_ITER}" -ge 0 ]]; then
        echo "No corrupted SMC files detected among newly checked files."

        if [[ "${LAST_VALID_ITER}" -ge "${TARGET_ITERS}" ]]; then
            echo "Run already completed successfully."
            exit 0
        fi

        echo "Will resume from iteration ${LAST_VALID_ITER}"
        RESUME_FLAG="--resume"

    else
        echo "No valid checkpoints found."
        echo "Starting fresh."
        RESUME_FLAG=""
    fi

else

    echo "First corrupted iteration detected at ${FIRST_CORRUPT_ITER}"
    echo "Last valid iteration = ${LAST_VALID_ITER:-NONE}"

    echo
    echo "Deleting corrupted and downstream files..."

    DELETE_COUNT=0

    for (( iter=FIRST_CORRUPT_ITER; iter<=TARGET_ITERS; iter+=10 )); do

        smc="${OUT_NAME}.${iter}.smc.gz"
        sites="${OUT_NAME}.${iter}.sites.gz"

        if [[ -f "$smc" ]]; then
            echo "Deleting: $smc"
            rm -f -- "$smc"
            ((++DELETE_COUNT))
        fi

        if [[ -f "$sites" ]]; then
            echo "Deleting: $sites"
            rm -f -- "$sites"
            ((++DELETE_COUNT))
        fi

    done

    echo "Finished deleting files."
    echo "Deleted ${DELETE_COUNT} files."

    if [[ -f "${LOGFILE}" ]]; then
        echo
        echo "Truncating ${LOGFILE} before iteration ${FIRST_CORRUPT_ITER}"

        TMP_LOG="${LOGFILE}.tmp"

        awk -v baditer="${FIRST_CORRUPT_ITER}" '
            /^sample [0-9]+$/ {
                split($0,a," ")
                if (a[2] >= baditer)
                    exit
            }
            { print }
        ' "${LOGFILE}" > "${TMP_LOG}"

        mv "${TMP_LOG}" "${LOGFILE}"
    fi

    if [[ -f "${STATSFILE}" ]]; then
        echo
        echo "Truncating ${STATSFILE} before iteration ${FIRST_CORRUPT_ITER}"

        TMP_STATS="${STATSFILE}.tmp"

        awk -v baditer="${FIRST_CORRUPT_ITER}" '
            $2 >= baditer { exit }
            { print }
        ' "${STATSFILE}" > "${TMP_STATS}"

        mv "${TMP_STATS}" "${STATSFILE}"
    fi

    echo
    echo "Recovery complete."

    if [[ -n "${LAST_VALID_ITER}" && "${LAST_VALID_ITER}" -ge 0 ]]; then
        echo "${LAST_VALID_ITER}" > "${VERIFY_FILE}"
        echo "Verification checkpoint reset to ${LAST_VALID_ITER}"
        echo "Will resume from iteration ${LAST_VALID_ITER}"
        RESUME_FLAG="--resume"
    else
        rm -f "${VERIFY_FILE}"
        echo "No valid iterations remain."
        echo "Removed verification checkpoint."
        echo "Starting fresh."
        RESUME_FLAG=""
    fi
fi

# Ensure .stats does not extend beyond the last valid SMC checkpoint.
# This handles restored .stats files that contain extra iterations
# beyond the currently available/valid .smc.gz files.
if [[ -f "${STATSFILE}" && -n "${LAST_VALID_ITER}" && "${LAST_VALID_ITER}" -ge 0 ]]; then
    echo
    echo "Checking whether ${STATSFILE} extends beyond LAST_VALID_ITER=${LAST_VALID_ITER}"

    TMP_STATS="${STATSFILE}.tmp"

    awk -v lastvalid="${LAST_VALID_ITER}" '
        NR == 1 { print; next }

        $1 != "resample" { print; next }

        $2 <= lastvalid { print; next }

        $2 > lastvalid { next }
    ' "${STATSFILE}" > "${TMP_STATS}"

    mv "${TMP_STATS}" "${STATSFILE}"

    echo "Finished trimming ${STATSFILE} to resample iterations <= ${LAST_VALID_ITER}"
fi

if [[ -f "${LOGFILE}" ]]; then
    LAST_SAMPLE=$(grep -E '^sample [0-9]+' "${LOGFILE}" | tail -n 1 | awk '{print $2}' || true)

    if [[ -n "${LAST_SAMPLE:-}" ]]; then
        echo "Last logged sample after cleanup = ${LAST_SAMPLE}"
    fi
fi

/home/amonc/ARGweaver/bin/arg-sample \
    --vcf "${VCF_DIR}/${VCF_FILE}" \
    --region "${REGION}" \
    --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
    --maskmap /scratch/amonc/xipho_revision/vcf_filtering/ARGweaver_mask.bed \
    --mask-cluster 2,5 \
    --mask-Ns 33 \
    --mutrate 4.6e-9 \
    --recombmap /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed \
    --popsize 391476 \
    --compress-seq 5 \
    --ntimes 20 \
    --maxtime 1e7 \
    --delta 0.005 \
    --sample-step 10 \
    --iters "${TARGET_ITERS}" \
    ${RESUME_FLAG} \
    -o "${OUT_NAME}"




## Here is a managable file to identify files to delete, May 11
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=list_corrupt_arg
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/list_corrupt_arg_%J.out
#SBATCH --error=logs/list_corrupt_arg_%J.err

set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
JOBFILE=${WORKDIR}/jobs_file.txt
TARGET_ITERS=3000
OUT_LIST=${WORKDIR}/corrupt_and_downstream_arg_files.txt
SUMMARY=${WORKDIR}/corrupt_and_downstream_arg_summary.txt

mkdir -p "${WORKDIR}/logs"

: > "${OUT_LIST}"
: > "${SUMMARY}"

smc_names_valid() {
    local smc="$1"
    local names_line
    local name

    names_line=$(zcat "$smc" 2>/dev/null | head -n 1 || true)

    if [[ ! "$names_line" =~ ^NAMES[[:space:]] ]]; then
        return 1
    fi

    for name in ${names_line#NAMES }; do
        if [[ "$name" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        if [[ ! "$name" =~ ^[A-Za-z0-9_-]+$ ]]; then
            return 1
        fi
    done

    return 0
}

cd "${OUTPUT_DIR}"

while IFS=$'\t' read -r VCF_FILE REGION OUT_NAME; do

    [[ -z "${OUT_NAME:-}" ]] && continue

    FIRST_CORRUPT_ITER=""

    echo "Checking ${OUT_NAME}"

    for (( iter=0; iter<=TARGET_ITERS; iter+=10 )); do
        smc="${OUT_NAME}.${iter}.smc.gz"

        [[ ! -f "$smc" ]] && continue

        if ! smc_names_valid "$smc"; then
            FIRST_CORRUPT_ITER="$iter"
            echo "${OUT_NAME} first_corrupt_iter=${FIRST_CORRUPT_ITER}" >> "${SUMMARY}"
            break
        fi
    done

    [[ -z "${FIRST_CORRUPT_ITER}" ]] && continue

    for (( iter=FIRST_CORRUPT_ITER; iter<=TARGET_ITERS; iter+=10 )); do
        smc="${OUTPUT_DIR}/${OUT_NAME}.${iter}.smc.gz"
        sites="${OUTPUT_DIR}/${OUT_NAME}.${iter}.sites.gz"

        [[ -f "$smc" ]] && printf '%s\n' "$smc" >> "${OUT_LIST}"
        [[ -f "$sites" ]] && printf '%s\n' "$sites" >> "${OUT_LIST}"
    done

done < "${JOBFILE}"

echo "Done."
echo "File list: ${OUT_LIST}"
echo "Summary: ${SUMMARY}"
echo "Total files listed: $(wc -l < "${OUT_LIST}")"




























# Started runs on April 18, 7PM (be sure to copy over ARGs to ourdisk by May 2)

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=ARGcheck
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ARGcheck_%A_%a.out
#SBATCH --error=logs/ARGcheck_%A_%a.err

ls *.10.smc.gz | wc -l # 628 as of 4/21 2:05PM

ls *.100.smc.gz | wc -l # 628 as of 4/21 2:04PM

ls *.120.smc.gz | wc -l # 628 as of 4/21 2:07PM

ls *.130.smc.gz | wc -l # 627 as of 4/21 2:07PM

ls *.500.smc.gz | wc -l # 539 as of 4/21 2:03PM

ls *.1000.smc.gz | wc -l # 366 as of 4/21 2:02PM

ls *.2000.smc.gz | wc -l  # 110 as of 4/21 11:55AM, 309 as of 4/23 9:07AM, 493 as of 4/27 8:42AM

ls *.3000.smc.gz | wc -l # 81 as of 4/21 11:55AM, 128 as of 4/23 8:43AM, 183 as of 4/24 12:41PM, 366 by 4/27 8:41AM, 372 by 4/28 8:42AM, 390 by 4/29 4:36, 408 by 5/2 7:51AM, 414 by 5/4 12:55pm; 455 by 5/8 8:07am




squeue -u amonc -t RUNNING -h | wc -l # 435 on 4/24 12:40PM
squeue -u amonc -t PENDING -h | wc -l

# checking details of jobs to see what worked
```
sacct -j 31476748 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ {print $2}' | \
sort | uniq -c
```
491 COMPLETED
143 RUNNING

# checking details of jobs to see what worked
```
sacct -j 31488441 --format=JobID,State -n | \
awk '$1 ~ /_[0-9]+$/ {print $2}' | \
sort | uniq -c
```
59 COMPLETED
58 FAILED
50 RUNNING


# Jobs getting stuck at Quarantining step? Best to leave for 48 hrs I think
```

Task 17
VCF_FILE=scaffold_1_30400000_32400000.vcf.gz
REGION=scaffold_1:30400001-32400000
OUT_NAME=scaffold1-17_out
Inspecting SMC checkpoints for scaffold1-17_out
VALID   scaffold1-17_out.0.smc.gz
VALID   scaffold1-17_out.10.smc.gz
VALID   scaffold1-17_out.20.smc.gz
VALID   scaffold1-17_out.30.smc.gz
VALID   scaffold1-17_out.40.smc.gz
VALID   scaffold1-17_out.50.smc.gz
VALID   scaffold1-17_out.60.smc.gz
VALID   scaffold1-17_out.70.smc.gz
VALID   scaffold1-17_out.80.smc.gz
VALID   scaffold1-17_out.90.smc.gz
VALID   scaffold1-17_out.100.smc.gz
VALID   scaffold1-17_out.110.smc.gz
VALID   scaffold1-17_out.120.smc.gz
VALID   scaffold1-17_out.130.smc.gz
VALID   scaffold1-17_out.140.smc.gz

[ETC]

VALID   scaffold1-17_out.2400.smc.gz
VALID   scaffold1-17_out.2410.smc.gz
VALID   scaffold1-17_out.2420.smc.gz
CORRUPT scaffold1-17_out.2430.smc.gz

LAST_VALID_ITER=2420
FIRST_CORRUPT_ITER=2430

First corrupted iteration detected at 2430
Last valid iteration = 2420

Quarantining corrupted and downstream files...
```



I have a list of stat files in a .txt file (empty_stats_confirmed.txt) in my argweaver directory (/scratch/amonc/xipho_revision/argweaver; example lines below):
scaffold1-2_out.stats
scaffold1-4_out.stats
scaffold1-5_out.stats
scaffold1-6_out.stats
scaffold1-11_out.stats
scaffold2-24_out.stats
scaffold2-28_out.stats
scaffold2-29_out.stats
scaffold2-30_out.stats
scaffold2-31_out.stats

Now, I want to use that list to identify and copy over files from my ourdisk storage folder (/ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/SMC_files) to the SMC_files folder in scratch (/scratch/amonc/xipho_revision/argweaver/SMC_files)



while read -r f; do
    [[ ! -s "$f" ]] && echo "$f"
done < stats_file_list.txt > confirmed_empty_stats.txt


I have a list of all stats files (/scratch/amonc/xipho_revision/argweaver/stats_file_list.txt): 
scaffold1-1_out.stats
scaffold1-2_out.stats
scaffold1-3_out.stats
scaffold1-4_out.stats
scaffold1-5_out.stats
scaffold1-6_out.stats
scaffold1-7_out.stats
scaffold1-8_out.stats
scaffold1-9_out.stats
scaffold1-10_out.stats

And I now want to check if any of these are missing or empty in my SMC_files folder (/scratch/amonc/xipho_revision/argweaver/SMC_files)






I restored some of the empty .stat files with the hopes of properly truncating them this time around. However, I am wondering what might happen if an smc file is all good up to, say, 2300 iterations (none corrupted detected because corrupted ones had been deleted previously). Yet in this case, my newly restored .stat files might contain more iterations still, say 2906. I would like those extra iterations to still be removed, so that the stat file matches the current stage of the smc file. Here is my current script:


# 11 May 2026. I am having a really difficult time with the all-in-one script (eg, ARGsv11.sbatch). I will try to do this more piecemeal. 
# This is the new pipeline! Much faster to break up into a few scripts

1) A script to identify all the smc files at and above the corrupted iteration for each region.

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=list_corrupt_arg
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/list_corrupt_arg_%J.out
#SBATCH --error=logs/list_corrupt_arg_%J.err

set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
JOBFILE=${WORKDIR}/jobs_file.txt
TARGET_ITERS=3000
OUT_LIST=${WORKDIR}/corrupt_and_downstream_arg_files.txt
SUMMARY=${WORKDIR}/corrupt_and_downstream_arg_summary.txt

mkdir -p "${WORKDIR}/logs"

: > "${OUT_LIST}"
: > "${SUMMARY}"

smc_names_valid() {
    local smc="$1"
    local names_line
    local name

    names_line=$(zcat "$smc" 2>/dev/null | head -n 1 || true)

    if [[ ! "$names_line" =~ ^NAMES[[:space:]] ]]; then
        return 1
    fi

    for name in ${names_line#NAMES }; do
        if [[ "$name" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        if [[ ! "$name" =~ ^[A-Za-z0-9_-]+$ ]]; then
            return 1
        fi
    done

    return 0
}

cd "${OUTPUT_DIR}"

while IFS=$'\t' read -r VCF_FILE REGION OUT_NAME; do

    [[ -z "${OUT_NAME:-}" ]] && continue

    FIRST_CORRUPT_ITER=""

    echo "Checking ${OUT_NAME}"

    for (( iter=0; iter<=TARGET_ITERS; iter+=10 )); do
        smc="${OUT_NAME}.${iter}.smc.gz"

        [[ ! -f "$smc" ]] && continue

        if ! smc_names_valid "$smc"; then
            FIRST_CORRUPT_ITER="$iter"
            echo "${OUT_NAME} first_corrupt_iter=${FIRST_CORRUPT_ITER}" >> "${SUMMARY}"
            break
        fi
    done

    [[ -z "${FIRST_CORRUPT_ITER}" ]] && continue

    for (( iter=FIRST_CORRUPT_ITER; iter<=TARGET_ITERS; iter+=10 )); do
        smc="${OUTPUT_DIR}/${OUT_NAME}.${iter}.smc.gz"
        sites="${OUTPUT_DIR}/${OUT_NAME}.${iter}.sites.gz"

        [[ -f "$smc" ]] && printf '%s\n' "$smc" >> "${OUT_LIST}"
        [[ -f "$sites" ]] && printf '%s\n' "$sites" >> "${OUT_LIST}"
    done

done < "${JOBFILE}"

echo "Done."
echo "File list: ${OUT_LIST}"
echo "Summary: ${SUMMARY}"
echo "Total files listed: $(wc -l < "${OUT_LIST}")"





# Now, I have have a file (corrupt_and_downstream_arg_summary.txt) that lists the scaffold base name and first corrupt iteration for all incomplete regions (example lines here): 
scaffold1-10_out first_corrupt_iter=2550
scaffold1-14_out first_corrupt_iter=2250
scaffold1-16_out first_corrupt_iter=2040
scaffold1-17_out first_corrupt_iter=2430
scaffold1-19_out first_corrupt_iter=2420
scaffold1-20_out first_corrupt_iter=2350 


I would like to take this information and use it to truncate the .log and .stat files in the SMC_files folder (/scratch/amonc/xipho_revision/argweaver/SMC_files). This new script will clean out the iterations that are corrupted and all following iterations based on the info in the summary file (/scratch/amonc/xipho_revision/argweaver/corrupt_and_downstream_arg_summary.txt). That is all I want the script to do.



# Now I want a script to delete all the files listed in this file: /scratch/amonc/xipho_revision/argweaver/corrupt_and_downstream_arg_files.txt
Basically I just need a for loop to go row by row and delete the file listed. Here are a few example lines of the file contents:
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2550.smc.gz
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2550.sites.gz
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2560.smc.gz
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2560.sites.gz
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2570.smc.gz
/scratch/amonc/xipho_revision/argweaver/SMC_files/scaffold1-10_out.2570.sites.gz


# Script to delete smc and sites files that are corrupted or of higher iteration than a corrupted file
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=delete_corrupt_arg
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/delete_corrupt_arg_%J.out
#SBATCH --error=logs/delete_corrupt_arg_%J.err

set -euo pipefail

FILELIST=/scratch/amonc/xipho_revision/argweaver/corrupt_and_downstream_arg_files.txt

if [[ ! -f "${FILELIST}" ]]; then
    echo "Missing file list: ${FILELIST}"
    exit 1
fi

DELETE_COUNT=0
MISSING_COUNT=0

while read -r file; do

    [[ -z "${file:-}" ]] && continue

    if [[ -f "${file}" ]]; then
        echo "Deleting: ${file}"
        rm -f -- "${file}"
        ((++DELETE_COUNT))
    else
        echo "Missing: ${file}"
        ((++MISSING_COUNT))
    fi

done < "${FILELIST}"

echo
echo "Finished."
echo "Deleted files: ${DELETE_COUNT}"
echo "Already missing: ${MISSING_COUNT}"








# This is the truncating script:

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=truncate_arg_logs
#SBATCH --chdir=/scratch/amonc/xipho_revision/argweaver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/truncate_arg_logs_%J.out
#SBATCH --error=logs/truncate_arg_logs_%J.err

set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/argweaver
OUTPUT_DIR=${WORKDIR}/SMC_files
SUMMARY=${WORKDIR}/corrupt_and_downstream_arg_summary.txt

mkdir -p "${WORKDIR}/logs"

if [[ ! -f "${SUMMARY}" ]]; then
    echo "Missing summary file: ${SUMMARY}"
    exit 1
fi

cd "${OUTPUT_DIR}"

while read -r OUT_NAME ITER_FIELD; do

    [[ -z "${OUT_NAME:-}" ]] && continue

    FIRST_CORRUPT_ITER="${ITER_FIELD#first_corrupt_iter=}"

    if [[ ! "${FIRST_CORRUPT_ITER}" =~ ^[0-9]+$ ]]; then
        echo "Skipping malformed line: ${OUT_NAME} ${ITER_FIELD}"
        continue
    fi

    LAST_GOOD_ITER=$(( FIRST_CORRUPT_ITER - 10 ))

    if (( LAST_GOOD_ITER < 0 )); then
        echo "Skipping ${OUT_NAME}: LAST_GOOD_ITER would be negative"
        continue
    fi

    LOGFILE="${OUT_NAME}.log"
    STATSFILE="${OUT_NAME}.stats"

    echo
    echo "Processing ${OUT_NAME}"
    echo "First corrupt iteration: ${FIRST_CORRUPT_ITER}"
    echo "Last good checkpoint iteration to keep: ${LAST_GOOD_ITER}"

    if [[ -f "${LOGFILE}" ]]; then
        echo "Truncating ${LOGFILE}"

        TMP_LOG="${LOGFILE}.tmp"

        awk -v lastgood="${LAST_GOOD_ITER}" '
            /^sample [0-9]+$/ {
                split($0,a," ")
                if (a[2] > lastgood)
                    exit
            }
            { print }
        ' "${LOGFILE}" > "${TMP_LOG}"

        mv "${TMP_LOG}" "${LOGFILE}"
    else
        echo "No log file found: ${LOGFILE}"
    fi

    if [[ -f "${STATSFILE}" ]]; then
        echo "Truncating ${STATSFILE}"

        TMP_STATS="${STATSFILE}.tmp"

        awk -v lastgood="${LAST_GOOD_ITER}" '
            NR == 1 { print; next }

            $1 != "resample" { print; next }

            $2 <= lastgood { print; next }

            $2 > lastgood { next }
        ' "${STATSFILE}" > "${TMP_STATS}"

        mv "${TMP_STATS}" "${STATSFILE}"
    else
        echo "No stats file found: ${STATSFILE}"
    fi

done < "${SUMMARY}"

echo
echo "Done truncating logs and stats."


# Pipeline for running ARGs every 2 days. Once ARGs finish, run the list corrupt to clean out bad files.

1) ARG_run.sbatch
2) list_corrupt_args.sbatch
3) delete_corrupt_arg.sbatch
4) truncate_arg_logs.sbatch
5) repeat steps 1-4



# count valid 3000s = 335, May 11 at 8:57PM; 337, May 12 11:33am

# Ok, still getting a bunch of "corrupt" files. And, weirdly, when I list corrupt args I seem to have regressed to iterations well prior (like 100s of iterations) before good iterations. Like the resume is jumping back to far. I don't understand. ARGweaver will never finish at this rate. 

# Next steps,
Can I find a way to just fix the header line in the smc.gz files.
First,
1) Is there any evidence of truncation on the smc.gz with weird line one
--smaller file size?
--not reaching the end of the region?

# Does the indexing in smc files match order from VCFs?
scaffold1-4_out.2400.smc.gz

bcftools query -l scaffold_219_0_158132.vcf.gz


# There is no clear consistency in the indices for the smc.gz files





# Preparing files on LSU cluster
# Removing 337 vcfs for completed regions. Remaining are 297 vcfs.
#!/bin/bash

# Files
MAP_FILE="/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/index_vcfname.txt"
COMPLETED_FILE="/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/completed_indices.txt"
VCF_DIR="/scratch/a_monc/postdoc/xipho_revision/argweaver/vcfs"

# Loop through completed indices
while read -r idx; do

    # Skip empty lines
    [[ -z "$idx" ]] && continue

    # Get corresponding VCF filename from mapping file
    vcf=$(awk -F'\t' -v i="$idx" '$1 == i {print $2}' "$MAP_FILE")

    # If found, delete the VCF
    if [[ -n "$vcf" ]]; then
        target="${VCF_DIR}/${vcf}"

        if [[ -f "$target" ]]; then
            echo "Deleting: $target"
            rm "$target"
        else
            echo "Missing file: $target"
        fi
    else
        echo "No VCF found for index: $idx"
    fi

done < "$COMPLETED_FILE"







# Setup for LSU jobs


I want to create a bunch of different .qsub files to run on the HPC cluster. There is key job information stored in a jobs file (/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/jobs_file.filtered.txt). I will store .qsub scripts in a separate folder (/scratch/a_monc/postdoc/xipho_revision/argweaver/qsub_scripts). The output argweaver output will be stored in yet another folder (/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output).

My working directory for the .qsub creation script (this one) is: /scratch/a_monc/postdoc/xipho_revision/argweaver

The content of jobs_file.filtered.txt is three columns. Column 1 is the VCF name. Column 2 is the region name. Column 3 is the out_name. Here are a few lines as an example:
scaffold_1_1900000_3900000.vcf.gz       scaffold_1:1900001-3900000      scaffold1-2_out
scaffold_1_5700000_7700000.vcf.gz       scaffold_1:5700001-7700000      scaffold1-4_out
scaffold_1_7600000_9600000.vcf.gz       scaffold_1:7600001-9600000      scaffold1-5_out
scaffold_1_9500000_11500000.vcf.gz      scaffold_1:9500001-11500000     scaffold1-6_out
scaffold_1_17100000_19100000.vcf.gz     scaffold_1:17100001-19100000    scaffold1-10_out
scaffold_1_19000000_21000000.vcf.gz     scaffold_1:19000001-21000000    scaffold1-11_out
scaffold_1_24700000_26700000.vcf.gz     scaffold_1:24700001-26700000    scaffold1-14_out
scaffold_1_28500000_30500000.vcf.gz     scaffold_1:28500001-30500000    scaffold1-16_out

I will denote these components of the jobs_file.filtered.txt as [vcf_name], [region_name], and [out_name] in my directions below. So here is the basic composition of each qsub file:
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N [out_name]

cd /scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output

arg-sample \
  --vcf "[vcf_name]" \
  --region "[region_name]" \
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \
  --maskmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ARGweaver_mask.bed \
  --mask-cluster 2,5 \
  --mask-Ns 33 \
  --mutrate 4.6e-9 \
  --recombmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ReLERNN_clean_data.bed \
  --popsize 391476 \
  --compress-seq 5 \
  --ntimes 20 \
  --maxtime 1e7 \
  --delta 0.005 \
  --sample-step 10 \
  --iters 3000 \
  -o [out_name]
```

Note, I will need the full path for each [vcf_name] entry.

## Making a shell script for LSU
#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/a_monc/postdoc/xipho_revision/argweaver
JOBFILE=/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/jobs_file.filtered.txt
VCF_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/vcfs
QSUB_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/qsub_scripts
ARG_OUT=/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output

mkdir -p "$QSUB_DIR"
mkdir -p "$ARG_OUT"

cd "$WORKDIR"

while read -r vcf_name region_name out_name; do
  [[ -z "${vcf_name:-}" ]] && continue

  vcf_path="${VCF_DIR}/${vcf_name}"
  qsub_file="${QSUB_DIR}/${out_name}.qsub"

  cat > "$qsub_file" <<EOF
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N ${out_name}

cd ${ARG_OUT}

arg-sample \\
  --vcf ${vcf_path} \\
  --region ${region_name} \\
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \\
  --maskmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ARGweaver_mask.bed \\
  --mask-cluster 2,5 \\
  --mask-Ns 33 \\
  --mutrate 4.6e-9 \\
  --recombmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ReLERNN_clean_data.bed \\
  --popsize 391476 \\
  --compress-seq 5 \\
  --ntimes 20 \\
  --maxtime 1e7 \\
  --delta 0.005 \\
  --sample-step 10 \\
  --iters 3000 \\
  -o ${out_name}
EOF

  chmod +x "$qsub_file"
  echo "Created: $qsub_file"

done < "$JOBFILE"



# Submitting

for q in scaffold*_out.qsub; do
  echo "Submitting $q"
  qsub "$q"
done



# new arg attempt on OSCER
# just for the vcfs needing to be processed

I need to transfer the vcfs from the LSU hpc: /scratch/a_monc/postdoc/xipho_revision/argweaver/vcfs

To the OSCER hpc: /scratch/amonc/xipho_revision/new_arg_attempt/vcfs

I want to use tmux



## Making a shell script for OSCER
#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/new_arg_attempt
JOBFILE=/scratch/amonc/xipho_revision/new_arg_attempt/jobs_file.filtered.txt
VCF_DIR=/scratch/amonc/xipho_revision/new_arg_attempt/vcfs
SBATCH_DIR=/scratch/amonc/xipho_revision/new_arg_attempt/sbatch_scripts
ARG_OUT=/scratch/amonc/xipho_revision/new_arg_attempt/arg_output

mkdir -p "$SBATCH_DIR"
mkdir -p "$ARG_OUT"

cd "$WORKDIR"

while read -r vcf_name region_name out_name; do
  [[ -z "${vcf_name:-}" ]] && continue

  vcf_path="${VCF_DIR}/${vcf_name}"
  sbatch_file="${SBATCH_DIR}/${out_name}.sbatch"

  cat > "$sbatch_file" <<EOF
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=${out_name}_%J_stdout.txt
#SBATCH --error=${out_name}_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=${out_name}

module load GCCcore/11.3.0

cd ${ARG_OUT}

/home/amonc/ARGweaver/bin/arg-sample \\
  --vcf ${vcf_path} \\
  --region ${region_name} \\
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \\
  --maskmap /scratch/amonc/xipho_revision/new_arg_attempt/ARGweaver_mask.bed \\
  --mask-cluster 2,5 \\
  --mask-Ns 33 \\
  --mutrate 4.6e-9 \\
  --recombmap /scratch/amonc/xipho_revision/new_arg_attempt/ReLERNN_clean_data.bed \\
  --popsize 391476 \\
  --compress-seq 5 \\
  --ntimes 20 \\
  --maxtime 1e7 \\
  --delta 0.005 \\
  --sample-step 10 \\
  --iters 3000 \\
  -o ${out_name}
EOF

  chmod +x "$sbatch_file"
  echo "Created: $sbatch_file"

done < "$JOBFILE"


# Nice, ok, now submitting jobs.
for s in scaffold*_out.sbatch; do
  sbatch "$s"
done


for s in batch_*.qsub; do
  sbatch "$s"
done

# When I resubmit on OSCER, I will need to use --resume



# I need a script to verify each genomic region in the SMC_files output folder (/scratch/amonc/xipho_revision/argweaver/SMC_files). The regions are denoted by the basename (3rd column) in the jobs file (/scratch/amonc/xipho_revision/argweaver/jobs_file.txt). Here are a few lines of the jobs file (tab-delimited):
scaffold_1_0_2000000.vcf.gz	scaffold_1:1-2000000	scaffold1-1_out
scaffold_1_1900000_3900000.vcf.gz	scaffold_1:1900001-3900000	scaffold1-2_out
scaffold_1_3800000_5800000.vcf.gz	scaffold_1:3800001-5800000	scaffold1-3_out
scaffold_1_5700000_7700000.vcf.gz	scaffold_1:5700001-7700000	scaffold1-4_out
scaffold_1_7600000_9600000.vcf.gz	scaffold_1:7600001-9600000	scaffold1-5_out

I want to verify that:
1) the genomic region has the 0th to the 3000th iteration smc.gz files (e.g., scaffold1-1_out.0.smc.gz through scaffold1-1_out.3000.smc.gz at intervals of 10)
2) each associated .stat file has 3000 iterations (second column final line, should provide this information)

To make this script faster, I think it makes sense to first check the .stat files for 3000 iterations. If 3000 iterations are present for a region, then check the smc.gz files for that region across all iterations to make sure 0th to 3000th iterations are not corrupted. If these items are verified, then I want separate lists of the filenames for the qualifying smc.gz files, sites.gz files, .stats files, .log files, and .masked_regions.bed files from these regions. Thus, five lists in all--just for the complete regions. These scripts can be placed in a new "verified_regions" folder.

Make sense?

I anticipate running this verification script (as an .sbatch script) from the argweaver directory (/scratch/amonc/xipho_revision/argweaver).

In a separate script I want to use these five files to rsync the listed files to ourdisk folders (for permanent storage).

# I have not run the above yet. # left at 5:53pm on May 12, 2026.



## Making a shell script for OSCER -- resume version
#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/amonc/xipho_revision/new_arg_attempt
JOBFILE=/scratch/amonc/xipho_revision/new_arg_attempt/jobs_file.filtered.txt
VCF_DIR=/scratch/amonc/xipho_revision/new_arg_attempt/vcfs
SBATCH_DIR=/scratch/amonc/xipho_revision/new_arg_attempt/sbatch_scripts_resume
ARG_OUT=/scratch/amonc/xipho_revision/new_arg_attempt/arg_output

mkdir -p "$SBATCH_DIR"
mkdir -p "$ARG_OUT"

cd "$WORKDIR"

while read -r vcf_name region_name out_name; do
  [[ -z "${vcf_name:-}" ]] && continue

  vcf_path="${VCF_DIR}/${vcf_name}"
  sbatch_file="${SBATCH_DIR}/${out_name}.sbatch"

  cat > "$sbatch_file" <<EOF
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=${out_name}_%J_stdout.txt
#SBATCH --error=${out_name}_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=${out_name}

module load GCCcore/11.3.0

cd ${ARG_OUT}

/home/amonc/ARGweaver/bin/arg-sample \\
  --vcf ${vcf_path} \\
  --region ${region_name} \\
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \\
  --maskmap /scratch/amonc/xipho_revision/new_arg_attempt/ARGweaver_mask.bed \\
  --mask-cluster 2,5 \\
  --mask-Ns 33 \\
  --mutrate 4.6e-9 \\
  --recombmap /scratch/amonc/xipho_revision/new_arg_attempt/ReLERNN_clean_data.bed \\
  --popsize 391476 \\
  --compress-seq 5 \\
  --ntimes 20 \\
  --maxtime 1e7 \\
  --delta 0.005 \\
  --sample-step 10 \\
  --iters 3000 \\
  --resume \\
  -o ${out_name}
EOF

  chmod +x "$sbatch_file"
  echo "Created: $sbatch_file"

done < "$JOBFILE"


find /scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output -maxdepth 1 -type f -name "*.stats" -printf "%f\n" > stats_file_list_already_started.txt
find /scratch/amonc/xipho_revision/argweaver/SMC_files \
  -maxdepth 1 -type f -name "*.stats" -printf "%f\n" \
  > stats_file_list_already_started.txt


  34 jobs running on lsu hpc
  29 are actually running

  Do I have any pending jobs on the LSU hpc?
qstat -u a_monc | grep ' Q ' | wc -l
qstat -u a_monc | grep ' R ' | wc -l

I have a couple of files that will help me with a new script. First, a jobs file--here are a few lines of the jobs file (tab-delimited; /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/jobs_file.filtered.txt):
scaffold_1_0_2000000.vcf.gz	scaffold_1:1-2000000	scaffold1-1_out
scaffold_1_1900000_3900000.vcf.gz	scaffold_1:1900001-3900000	scaffold1-2_out
scaffold_1_3800000_5800000.vcf.gz	scaffold_1:3800001-5800000	scaffold1-3_out
scaffold_1_5700000_7700000.vcf.gz	scaffold_1:5700001-7700000	scaffold1-4_out
scaffold_1_7600000_9600000.vcf.gz	scaffold_1:7600001-9600000	scaffold1-5_out

The third column of the jobs file has base names for different genomic regions. 

I have another file (/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/stats_file_list_already_started.txt), with names of .stats file for genomic regions I have already started running. Here are a few lines from this file:
scaffold10-141_out.stats
scaffold102-557_out.stats
scaffold10-142_out.stats
scaffold10-140_out.stats
scaffold10-137_out.stats
scaffold10-143_out.stats

I want to generate a new list that has only the regions found in the jobs file but not in the stats_file_list_already_started.txt file. Also, I want this list to be formatted as the basename.qsub. (Next I will need a command to submit all of the .qsub jobs)

# Nice, ok, now submitting jobs.
for s in scaffold9-*_out.qsub; do
  qsub "$s"
done


for s in scaffold*_out.sbatch; do
  sbatch "$s"
done

for s in batch*.qsub; do
  qsub "$s"
done

hpc_thomlab02


# Next up, run the rsync of verified regions, which should give me 337 clean regions in ourdisk




I am getting corrupted ARGs on OSCER but not on the LSU HPC. I am running the scripts the same way pretty much. The script specifications are different though.

#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N ${out_name}


vs.

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=${out_name}_%J_stdout.txt
#SBATCH --error=${out_name}_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=${out_name}

What is the memory allocation to the lsu jobs?

How can I check the actual memory usage of an lsu job?



# Completed jobs on LSU HPC
May 17, 2026, 20 files at 3000 iters



I want to run a parallel job on Mike. I think this will allow me to pack more simultaneous argweaver jobs on the cluster. Right now, I'm maxing out at 64 simulatanous jobs, which are separate .qsub scripts. I think I can do better.

#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/a_monc/postdoc/xipho_revision/argweaver
JOBFILE=/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/jobs_file.filtered.txt
VCF_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/vcfs
QSUB_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/qsub_scripts
ARG_OUT=/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output

mkdir -p "$QSUB_DIR"
mkdir -p "$ARG_OUT"

cd "$WORKDIR"

while read -r vcf_name region_name out_name; do
  [[ -z "${vcf_name:-}" ]] && continue

  vcf_path="${VCF_DIR}/${vcf_name}"
  qsub_file="${QSUB_DIR}/${out_name}.qsub"

  cat > "$qsub_file" <<EOF
#!/bin/bash
#PBS -A hpc_thomlab02
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N ${out_name}

cd ${ARG_OUT}

arg-sample \\
  --vcf ${vcf_path} \\
  --region ${region_name} \\
  --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \\
  --maskmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ARGweaver_mask.bed \\
  --mask-cluster 2,5 \\
  --mask-Ns 33 \\
  --mutrate 4.6e-9 \\
  --recombmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ReLERNN_clean_data.bed \\
  --popsize 391476 \\
  --compress-seq 5 \\
  --ntimes 20 \\
  --maxtime 1e7 \\
  --delta 0.005 \\
  --sample-step 10 \\
  --iters 3000 \\
  -o ${out_name}
EOF

  chmod +x "$qsub_file"
  echo "Created: $qsub_file"

done < "$JOBFILE"


This is what I have available on the cluster for parallel (module avail command):
parallel-netcdf/1.12.2/intel-2021.5.0                     
parallel-netcdf/1.12.2/intel-2021.5.0-intel-mpi-2021.5.1  
parallel/20210922/intel-2021.5.0  


I want to keep the same basic setup for my runs. Also, I just want to add a simple check for a 0 iteration smc file (which, if present, triggers the addition of the --resume flag to teh arg-sample command)



# New Mike script




#!/bin/bash
set -euo pipefail

WORKDIR=/scratch/a_monc/postdoc/xipho_revision/argweaver
JOBFILE=/ddnA/work/a_monc/postdoc/xipho_revision/argweaver/jobs_file.filtered.txt
VCF_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/vcfs
QSUB_DIR=/scratch/a_monc/postdoc/xipho_revision/argweaver/qsub_scripts_parallel
ARG_OUT=/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output

PPN=64
PARALLEL_JOBS=16
REGIONS_PER_QSUB=16

mkdir -p "$QSUB_DIR" "$ARG_OUT"
cd "$WORKDIR"

batch=1
count=0
cmdfile=""

while read -r vcf_name region_name out_name; do
  [[ -z "${vcf_name:-}" ]] && continue

  if (( count % REGIONS_PER_QSUB == 0 )); then
    batch_pad=$(printf "%04d" "$batch")
    cmdfile="${QSUB_DIR}/batch_${batch_pad}.commands.tsv"
    qsub_file="${QSUB_DIR}/batch_${batch_pad}.qsub"

    : > "$cmdfile"

    cat > "$qsub_file" <<EOF
#!/bin/bash
#PBS -A hpc_thomlab02
#PBS -l nodes=1:ppn=${PPN}
#PBS -l walltime=48:00:00
#PBS -q workq
#PBS -N arg_batch_${batch_pad}
#PBS -o ${QSUB_DIR}/batch_${batch_pad}.out
#PBS -e ${QSUB_DIR}/batch_${batch_pad}.err

set -euo pipefail

module load parallel/20210922/intel-2021.5.0

cd ${WORKDIR}

parallel --colsep '\\t' --jobs ${PARALLEL_JOBS} --joblog ${QSUB_DIR}/batch_${batch_pad}.parallel.log '
  vcf_path={1}
  region_name={2}
  out_name={3}
  outdir={4}

  mkdir -p "\$outdir"
  cd "\$outdir"

  resume_flag=""
  if [[ -s "\${out_name}.0.smc.gz" ]]; then
    resume_flag="--resume"
  fi

  echo "Starting \$out_name in \$outdir with resume_flag=\$resume_flag"

  arg-sample \\
    --vcf "\$vcf_path" \\
    --region "\$region_name" \\
    --vcf-genotype-filter "DP<5;DP>50;GQ<20;RGQ<20" \\
    --maskmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ARGweaver_mask.bed \\
    --mask-cluster 2,5 \\
    --mask-Ns 33 \\
    --mutrate 4.6e-9 \\
    --recombmap /ddnA/work/a_monc/postdoc/xipho_revision/argweaver/ReLERNN_clean_data.bed \\
    --popsize 391476 \\
    --compress-seq 5 \\
    --ntimes 20 \\
    --maxtime 1e7 \\
    --delta 0.005 \\
    --sample-step 10 \\
    --iters 3000 \\
    \$resume_flag \\
    -o "\$out_name"
' :::: ${cmdfile}
EOF

    chmod +x "$qsub_file"
    echo "Created qsub: $qsub_file"

    ((++batch))
  fi

  vcf_path="${VCF_DIR}/${vcf_name}"
  outdir="${ARG_OUT}/${out_name}"

  printf "%s\t%s\t%s\t%s\n" \
    "$vcf_path" "$region_name" "$out_name" "$outdir" >> "$cmdfile"

  ((++count))

done < "$JOBFILE"

echo "Created $((batch - 1)) parallel qsub scripts in:"
echo "$QSUB_DIR"



# Test run to reorganize folders
#!/bin/bash
set -euo pipefail

ARG_OUT=/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output

cd "$ARG_OUT"

for file in *; do
    [[ -f "$file" ]] || continue

    # Extract region basename
    # Keeps:
    # scaffold1-1_out
    # from:
    # scaffold1-1_out.0.smc.gz
    # scaffold1-1_out.stats
    # scaffold1-1_out.log
    # scaffold1-1_out.sites.gz
    # scaffold1-1_out.masked_regions.bed
    # etc.

    base=$(echo "$file" | sed -E 's/^(.+_out)\..*$/\1/')

    # Skip files that do not match expected ARGweaver naming
    [[ "$base" == *_out ]] || continue

    mkdir -p "$base"

    echo mv "$file" "$base/"
done




# Some manual attempts
/scratch/a_monc/postdoc/xipho_revision/arg_manual_attempts/

cp -r /scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output/scaffold2-33_out /scratch/a_monc/postdoc/xipho_revision/arg_manual_attempts/
cp scaffold2-33_out.qsub /scratch/a_monc/postdoc/xipho_revision/arg_manual_attempts/scaffold2-33_out

I need a quick command to delete all smc and sites files in a folder with iter = 2920 or higher.

find . -maxdepth 1 -type f \( -name "*.smc.gz" -o -name "*.sites.gz" \) \
| awk -F. '
{
    iter=$(NF-2)
    if (iter >= 1000)
        print
}' \
| xargs -r rm -f



# Verification saving script
#!/bin/bash
#SBATCH --job-name=verify_arg_regions
#SBATCH --account=hpc_thomlab02
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=verify_arg_regions_%j.out
#SBATCH --error=verify_arg_regions_%j.err
#SBATCH --chdir=/scratch/a_monc/postdoc/xipho_revision/argweaver

set -euo pipefail
shopt -s nullglob

BASE_DIR="/scratch/a_monc/postdoc/xipho_revision/argweaver"
ARG_OUT="${BASE_DIR}/arg_output"

VERIFY_DIR="${BASE_DIR}/verified_output"
SMC_DEST="${VERIFY_DIR}/verified_smcs"
SITES_DEST="${VERIFY_DIR}/verified_sites"
STATS_DEST="${VERIFY_DIR}/verified_stats"
LOGS_DEST="${VERIFY_DIR}/verified_logs"
BEDS_DEST="${VERIFY_DIR}/verified_beds"

COMPLETED_LIST="${VERIFY_DIR}/completed_uncorrupted_regions.txt"
FAILED_LIST="${VERIFY_DIR}/failed_or_corrupted_regions.txt"

mkdir -p "$SMC_DEST" "$SITES_DEST" "$STATS_DEST" "$LOGS_DEST" "$BEDS_DEST"

: > "$COMPLETED_LIST"
: > "$FAILED_LIST"

is_corrupt_smc () {
    local smc="$1"

    local names_line
    names_line=$(gzip -cd "$smc" 2>/dev/null | sed -n '1p') || return 0

    [[ -z "$names_line" ]] && return 0

    # Remove leading NAMES if present, then test each name.
    # Corrupt = any name that is solely numeric.
    echo "$names_line" | awk '
    {
        start = 1
        if ($1 == "NAMES") start = 2
        for (i = start; i <= NF; i++) {
            if ($i ~ /^[0-9]+$/) {
                exit 1
            }
        }
        exit 0
    }'
}

for region_dir in "$ARG_OUT"/*; do
    [[ -d "$region_dir" ]] || continue

    basename_region=$(basename "$region_dir")
    echo "Checking $basename_region"

    good=1
    reason=""

    for iter in $(seq 0 3000); do
        smc="${region_dir}/${basename_region}.${iter}.smc.gz"

        if [[ ! -s "$smc" ]]; then
            good=0
            reason="missing_smc_iter_${iter}"
            break
        fi

        if ! is_corrupt_smc "$smc"; then
            good=0
            reason="corrupt_smc_iter_${iter}"
            break
        fi
    done

    if [[ "$good" -eq 1 ]]; then
        echo "$basename_region" >> "$COMPLETED_LIST"

        mkdir -p "${SMC_DEST}/${basename_region}"
        mkdir -p "${SITES_DEST}/${basename_region}"

        cp -p "${region_dir}"/*.smc.gz "${SMC_DEST}/${basename_region}/"

        sites_files=( "${region_dir}"/*.sites.gz )
        if (( ${#sites_files[@]} > 0 )); then
            cp -p "${sites_files[@]}" "${SITES_DEST}/${basename_region}/"
        fi

        stats_files=( "${region_dir}"/*.stats )
        if (( ${#stats_files[@]} > 0 )); then
            cp -p "${stats_files[@]}" "$STATS_DEST/"
        fi

        log_files=( "${region_dir}"/*.log )
        if (( ${#log_files[@]} > 0 )); then
            cp -p "${log_files[@]}" "$LOGS_DEST/"
        fi

        bed_files=( "${region_dir}"/*.bed )
        if (( ${#bed_files[@]} > 0 )); then
            cp -p "${bed_files[@]}" "$BEDS_DEST/"
        fi

        echo "  COMPLETE: copied files"

    else
        echo -e "${basename_region}\t${reason}" >> "$FAILED_LIST"
        echo "  FAILED: $reason"
    fi
done

echo
echo "Done."
echo "Completed regions:"
wc -l "$COMPLETED_LIST"

echo "Failed/corrupt regions:"
wc -l "$FAILED_LIST"


# Nice, ok, now submitting jobs.
for s in scaffold*_out.qsub; do
  qsub "$s"
done


for s in scaffold*_out.sbatch; do
  sbatch "$s"
done

for s in batch*.qsub; do
  qsub "$s"
done

hpc_thomlab02


# Transfer of files from Mike to OSCER
SRC=/scratch/a_monc/postdoc/xipho_revision/argweaver/arg_output
REMOTE=amonc@schooner.oscer.ou.edu
DEST=/ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/argweaver_2000iters
LOGDIR=$SRC/transfer_logs_2000

ssh $REMOTE "mkdir -p $DEST/verified_smcs $DEST/verified_sites $DEST/verified_stats $DEST/verified_logs $DEST/verified_beds"

# Get masked percentages from .log files. Save below as bash file.

#!/bin/bash
set -euo pipefail

LOG_DIR="/ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/argweaver_2000iters/verified_logs"
OUT_FILE="/ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/argweaver/masked_percentages.txt"

tmpfile=$(mktemp)

for log in "$LOG_DIR"/*.log; do
    [[ -e "$log" ]] || continue

    logfile=$(basename "$log")

    awk -v file="$logfile" '
        /^masked [0-9]+ \([0-9.]+%\) positions/ {
            positions = $2
            percentage = $3

            gsub(/[()%]/, "", percentage)

            print positions "\t" percentage "\t" file
            exit
        }
    ' "$log" >> "$tmpfile"
done

{
    echo -e "masked_positions\tmasked_percentage\tlog_file"
    sort -k2,2nr "$tmpfile"
} > "$OUT_FILE"

rm -f "$tmpfile"

echo "Saved sorted results to: $OUT_FILE"