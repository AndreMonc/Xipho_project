## For snpArcher run, I need to download all X. spixii samples onto OU cluster

# First, Download SRA toolkit
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-alma_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz
# I added bin folder to my bash profile path
which fastq-dump

# Running a slurm array to download all 33 individuals from NCBI
#!/bin/bash
#SBATCH --account=general
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --job-name=sra_dl
#SBATCH --output=logs/sra_dl_%A_%a.out
#SBATCH --error=logs/sra_dl_%A_%a.err
#SBATCH --array=1-33

set -euo pipefail

# -------- user settings --------
ACCESSION_LIST="spixii_sra_list.txt"
BASE_DIR="/scratch/amonc/sra_project"
SRA_DIR="${BASE_DIR}/sra_cache"
FASTQ_DIR="${BASE_DIR}/fastq"
TMP_DIR="${BASE_DIR}/tmp/${SLURM_ARRAY_TASK_ID}"
THREADS="${SLURM_CPUS_PER_TASK}"
# --------------------------------

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR"

# Load tools if available
# module purge
# module load pigz
module load SRA-Toolkit/3.0.3-gompi-2022a
# Uncomment those if OSCER provides them

# If using your own sratoolkit install:
#export PATH=/scratch/amonc/sratoolkit.3.3.0-centos_linux64/bin:$PATH

# Get accession for this array task
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSION_LIST")

echo "[$(date)] Starting accession: $ACC"
echo "Threads: $THREADS"

# Download the .sra file
prefetch "$ACC" \
    --output-directory "$SRA_DIR"

# Convert to FASTQ
# Using --split-files for paired-end output
# Using --skip-technical to avoid technical reads
fasterq-dump "$SRA_DIR/$ACC" \
    --split-files \
    --skip-technical \
    --threads "$THREADS" \
    --temp "$TMP_DIR" \
    --outdir "$FASTQ_DIR"

# Compress
if [[ -f "$FASTQ_DIR/${ACC}_1.fastq" ]]; then
    pigz -p "$THREADS" "$FASTQ_DIR/${ACC}_1.fastq"
fi

if [[ -f "$FASTQ_DIR/${ACC}_2.fastq" ]]; then
    pigz -p "$THREADS" "$FASTQ_DIR/${ACC}_2.fastq"
fi

if [[ -f "$FASTQ_DIR/${ACC}.fastq" ]]; then
    pigz -p "$THREADS" "$FASTQ_DIR/${ACC}.fastq"
fi

# Remove downloaded .sra after successful conversion
rm -rf "$SRA_DIR/$ACC"

# Clean temp
rm -rf "$TMP_DIR"

echo "[$(date)] Finished accession: $ACC"

# Make executable and then run
chmod +x submit_sra_array.sbatch 
sbatch submit_sra_array.sbatch

# Get read counts
#!/bin/bash
#SBATCH --account=general
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --job-name=fastq_stats
#SBATCH --output=logs/fastq_stats_%A_%a.out
#SBATCH --error=logs/fastq_stats_%A_%a.err

cd /scratch/amonc/sra_project/fastq

seqkit stats *.fastq.gz -T > read_counts.tsv


# Download missing fastq file

#!/bin/bash
#SBATCH --account=general
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --job-name=last_download
#SBATCH --output=logs/last_download_%A_%a.out
#SBATCH --error=logs/last_download_%A_%a.err

module load SRA-Toolkit/3.0.3-gompi-2022a

# 1) Download the .sra file
prefetch SRR27456428 --max-size 50G

# 2) Move into the downloaded SRA directory
cd SRR27456428

# 3) Convert to paired fastq files
fasterq-dump SRR27456428 \
  --split-files \
  --skip-technical \
  --threads 4 \
  --temp /scratch/amonc/xipho_revision/xspixii_download/tmp

# 4) Compress the fastq files to fastq.gz
pigz -p 4 SRR27456428_1.fastq SRR27456428_2.fastq

# executable
chmod +x last_download.sbatch 


# Check last downloaded files
seqkit stats SRR27456428*.fastq.gz -T > read_counts_last.tsv