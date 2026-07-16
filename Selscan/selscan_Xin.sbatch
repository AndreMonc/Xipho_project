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
