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
