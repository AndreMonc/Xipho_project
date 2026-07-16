# 4/22/26
# Implementation of haplotype-based selection scan on Xiphorhynchus data
# Xiphorhynchus revision

# Download and compile selscan v3.0.0
module load GCCcore/11.3.0
git clone https://github.com/szpiech/selscan.git
cd selscan && git checkout main
cd src && make

# 
I want a script that runs selscan for three subprograms (--nsl --ihs --xpnsl). For --xpnsl, my reference population is Tapajos and my alternate population is Xingu. Since each of these populations has separate VCFs (1 for each scaffold), it is critical that the VCF for the same scaffold is run for each population.

All my vcfs for Xingu (1 per scaffold), are here: /scratch/amonc/xipho_revision/selscan/Xin_vcfs
All my vcfs for Tapajos (1 per scaffold), are here: /scratch/amonc/xipho_revision/selscan/Tap_vcfs

I don't expect the number of scaffolds to be the same for each population, so it won't be possible to run --xpnsl for all scaffolds (not all are paired between populations). 

Additionally, I want to run the "--unphased" flag, since my VCFs are unphased. The "--out" flag can be "Xin" to indicate that I am focused on Xingu populations across these three statistics. I also want to normalize all the output files. I think I'd like to do a separate normalization first (i.e. "selscan norm [--ihs|--nsl|--xpnsl]" for all output files) and then windowing with 50kb windows (i.e., "--norm-files *.out.norm --bp-win --winsize 50000")

(I've attached the selscan manual and a recent usage guide paper)

# Xingu nSL and iHS runs
```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=selscan_Xin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/scratch/amonc/xipho_revision/selscan
#SBATCH --output=logs/selscan_Xin_%j.out
#SBATCH --error=logs/selscan_Xin_%j.err

module load GCCcore/11.3.0

set -eu
shopt -s nullglob

THREADS=${SLURM_CPUS_PER_TASK:-8}

XIN_DIR=/scratch/amonc/xipho_revision/selscan/Xin_vcfs

OUTDIR=/scratch/amonc/xipho_revision/selscan/Xin_selscan_results
RAW=${OUTDIR}/raw
NORM=${OUTDIR}/norm
WINDOWS=${OUTDIR}/windows

mkdir -p "$RAW" "$NORM" "$WINDOWS" logs

# Function: get the scaffold name from a one-scaffold VCF
get_scaffold () {
    local vcf="$1"
    bcftools view -H "$vcf" | awk 'NF {print $1; exit}'
}

echo "Running selscan for Xingu nSL and iHS..."

for xin_vcf in "$XIN_DIR"/*.vcf.gz; do
    scaff=$(get_scaffold "$xin_vcf")
    safe_scaff=$(echo "$scaff" | sed 's/[^A-Za-z0-9_.-]/_/g')

    echo "Xingu scaffold: $scaff"

    # nSL
    selscan \
        --nsl \
        --vcf "$xin_vcf" \
        --unphased \
        --threads "$THREADS" \
        --out "${RAW}/Xin.${safe_scaff}"

    # iHS
    selscan \
        --ihs \
        --vcf "$xin_vcf" \
        --pmap \
        --unphased \
        --threads "$THREADS" \
        --out "${RAW}/Xin.${safe_scaff}"
done

echo "Normalizing all output files together by statistic..."

cp "$RAW"/*.{ihs,nsl}.out "$NORM"/ 2>/dev/null || true

cd "$NORM"

# Normalize raw outputs
ihs_files=( *.ihs.out )
nsl_files=( *.nsl.out )

if (( ${#ihs_files[@]} > 0 )); then
    selscan norm --ihs --files "${ihs_files[@]}" --bins 100
fi

if (( ${#nsl_files[@]} > 0 )); then
    selscan norm --nsl --files "${nsl_files[@]}" --bins 100
fi

echo "Creating 10 kb windows from normalized files..."

ihs_norm=( *.ihs.out.100bins.norm )
nsl_norm=( *.nsl.out.100bins.norm )

if (( ${#ihs_norm[@]} > 0 )); then
    selscan norm --ihs --norm-files "${ihs_norm[@]}" --bp-win --winsize 10000
fi

if (( ${#nsl_norm[@]} > 0 )); then
    selscan norm --nsl --norm-files "${nsl_norm[@]}" --bp-win --winsize 10000
fi

mv *.windows "$WINDOWS"/ 2>/dev/null || true

echo "Done."
echo "Raw outputs:        $RAW"
echo "Normalized outputs: $NORM"
echo "Window outputs:     $WINDOWS"
```



# Belem nSL and iHS runs
```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=selscan_Bel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/scratch/amonc/xipho_revision/selscan
#SBATCH --output=logs/selscan_Bel_%j.out
#SBATCH --error=logs/selscan_Bel_%j.err

module load GCCcore/11.3.0

set -eu
shopt -s nullglob

THREADS=${SLURM_CPUS_PER_TASK:-8}

Bel_DIR=/scratch/amonc/xipho_revision/selscan/Bel_vcfs
TAP_DIR=/scratch/amonc/xipho_revision/selscan/Tap_vcfs

OUTDIR=/scratch/amonc/xipho_revision/selscan/Bel_selscan_results
RAW=${OUTDIR}/raw
NORM=${OUTDIR}/norm
WINDOWS=${OUTDIR}/windows
LOGS=${OUTDIR}/logs

mkdir -p "$RAW" "$NORM" "$WINDOWS" "$LOGS" logs

# Function: get the scaffold name from a one-scaffold VCF
get_scaffold () {
    local vcf="$1"
    bcftools view -H "$vcf" | awk 'NF {print $1; exit}'
}

declare -A TAP_BY_SCAFFOLD

echo "Indexing Tapajos VCFs by scaffold..."
for tap_vcf in "$TAP_DIR"/*.vcf.gz; do
    scaff=$(get_scaffold "$tap_vcf")
    TAP_BY_SCAFFOLD["$scaff"]="$tap_vcf"
done

echo "Running selscan for Belem nSL, Belem iHS, and paired Bel-vs-Tap XPnSL..."

for Bel_vcf in "$Bel_DIR"/*.vcf.gz; do
    scaff=$(get_scaffold "$Bel_vcf")
    safe_scaff=$(echo "$scaff" | sed 's/[^A-Za-z0-9_.-]/_/g')

    echo "Belem scaffold: $scaff"

    # nSL: Belem only
    selscan \
        --nsl \
        --vcf "$Bel_vcf" \
        --unphased \
        --threads "$THREADS" \
        --out "${RAW}/Bel.${safe_scaff}"

    # iHS: Belem only, using physical map because no genetic map available
    selscan \
        --ihs \
        --vcf "$Bel_vcf" \
        --pmap \
        --unphased \
        --threads "$THREADS" \
        --out "${RAW}/Bel.${safe_scaff}"

    # XPnSL: only if Tapajos has the same scaffold
    if [[ -n "${TAP_BY_SCAFFOLD[$scaff]+set}" ]]; then
        tap_vcf="${TAP_BY_SCAFFOLD[$scaff]}"

        selscan \
            --xpnsl \
            --vcf-ref "$tap_vcf" \
            --vcf "$Bel_vcf" \
            --unphased \
            --threads "$THREADS" \
            --out "${RAW}/Bel.${safe_scaff}"
    else
        echo "Skipping XPnSL for $scaff: no matching Tapajos VCF found." \
            | tee -a "${LOGS}/unpaired_scaffolds.txt"
    fi
done

echo "Normalizing all output files together by statistic..."

cp "$RAW"/*.{ihs,nsl,xpnsl}.out "$NORM"/ 2>/dev/null || true

cd "$NORM"

# Normalize raw outputs
ihs_files=( *.ihs.out )
nsl_files=( *.nsl.out )
xpnsl_files=( *.xpnsl.out )

if (( ${#ihs_files[@]} > 0 )); then
    selscan norm --ihs --files "${ihs_files[@]}" --bins 100
fi

if (( ${#nsl_files[@]} > 0 )); then
    selscan norm --nsl --files "${nsl_files[@]}" --bins 100
fi

if (( ${#xpnsl_files[@]} > 0 )); then
    selscan norm --xpnsl --files "${xpnsl_files[@]}"
fi

echo "Creating 10 kb windows from normalized files..."

ihs_norm=( *.ihs.out.norm )
nsl_norm=( *.nsl.out.norm )
xpnsl_norm=( *.xpnsl.out.norm )

if (( ${#ihs_norm[@]} > 0 )); then
    selscan norm --ihs --norm-files "${ihs_norm[@]}" --bp-win --winsize 10000
fi

if (( ${#nsl_norm[@]} > 0 )); then
    selscan norm --nsl --norm-files "${nsl_norm[@]}" --bp-win --winsize 10000
fi

if (( ${#xpnsl_norm[@]} > 0 )); then
    selscan norm --xpnsl --norm-files "${xpnsl_norm[@]}" --bp-win --winsize 10000
fi

mv *.windows "$WINDOWS"/ 2>/dev/null || true

echo "Done."
echo "Raw outputs:        $RAW"
echo "Normalized outputs: $NORM"
echo "Window outputs:     $WINDOWS"
```





-------------
# Notes: 
# XP-nSL

# iHH12 (can't take unphased data, so not using)

# Notes on selscan usage

I can use Tapajos as reference pop and Bel or Xin as the object population for XP-nSL

genotype depth of 5-50x (Genomic Evidence for Elevational Segregation and Adaptive Introgression in Prunellidae Radiation; Mol. Ecology 2025)

XP-nSL, iHS and nSL all run; Used (--pmap) to calculate IHS (Genomic Footprints of Hybridisation in North Atlantic Eels (Anguilla anguilla and A. rostrata; Molecular ecology 2025); sliding windows of 10kb, min. of 10 snps; overlapping regions of iHS and nSL merged to create list of candidate regions under selection

Potential Adaptive Introgression From Dogs in Iberian Grey Wolves (Canis lupus) -- Molecular Ecology
"Second, we estimated XP-nSL to identify sites under selection but also present in several different haplotypes, thus ignoring long stretches of homozygosity (Szpiech et al. 2021). We normalised the results with the norm function of the selscan package and only reported sites that were found to be putatively under selection after normalisation (‘crit’ value = 1, corresponding to Z- score > |2|; Voight et al. 2006)."

# Things to keep in mind
Run separately for each chromosome (or scaffold)
No missing data
Normalize data
--run iHS with --pmap
--unphased flag for unphased data (default is false)
selscan default is 100kb! Very rough resolution, run on 50kb windows (lowest recommended in 2026 guide paper for dense snp datasets)

# Several options with different pros and cons

# nSL (Ferrer-Admetlla et al. 2014)
A Ferrer-Admetlla et al. (2014) On detecting incomplete soft or hard selective sweeps
	using haplotype structure. Molecular Biology and Evolution 31: 1275-1291.

- single population statistic
- less power with fewer than < 50 individuals (see Szpiech 2024)
- performs well with unphased data
- does not require genetic map information

To calculate nSL:
./selscan --nsl --vcf <vcf> --out <outfile>

# Step 1, run nsl
selscan --nsl --vcf mydata.vcf.gz --out genome
# Step 2, run normalization
selscan norm --nsl --files chr*.nsl.out --bins 100 # normalization step

# XP-nSL (Szpiech et al. 2021)
ZA Szpiech et al. (2021) Application of a novel haplotype-based scan for local adaptation 
	to study high-altitude adaptation in rhesus macaques. Evolution Letters 
	doi: https://doi.org/10.1002/evl3.232

- More power than nSL, working well with as few as 10 individuals (see Szpiech 2024)
To calculate XP-nSL:
./selscan --xpnsl --vcf <vcf> --vcf-ref <vcf> --out <outfile>

# Step 1, run xpnsl
selscan --xpnsl --vcf popA.vcf.gz --vcf-ref popB.vcf.gz
# Step 2, run normalization
selscan norm --xpnsl --files *.xpnsl.out


# iHS (Voight et al., 2006)
BF Voight et al. (2006) A map of recent positive selection in the human 
	genome. PLoS Biology 4: e72.
- single population statistic
- requires large sample size >100 individuals (see Szpiech 2024)

# iHH12 also looks good--for soft sweeps; no unphased version


# Pages 14-15 of Selection scans and downstream analysis with selscan (2026) is really helpful for thinking about value of each statistic.


# Also, paper on using with unphased data (Szpiech 2024): 
ZA Szpiech (2024) selscan 2.0: scanning for sweeps in unphased data. Bioinformatics, 40(1), btae006.
	doi: https://doi.org/10.1093/bioinformatics/btae006






