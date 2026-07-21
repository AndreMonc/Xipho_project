# PIXY for Coari project with Gymnopithys and Phlegopsis
# conda 26.1.1
conda create --name pixy
conda activate pixy

conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib

version 2.0.0.beta14


# pixy .sbatch file, run on OU HPC
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=pixy
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=pixy_%j.out
#SBATCH --error=pixy_%j.err
#SBATCH --chdir=/scratch/amonc/xipho_revision/fst_dxy_pi

source ~/.bashrc
conda activate pixy

set -eo pipefail

OUTDIR=/scratch/amonc/xipho_revision/fst_dxy_pi/pixy
mkdir -p "$OUTDIR"

test -w "$OUTDIR" || {
  echo "ERROR: pixy output folder is not writable: $OUTDIR"
  exit 1
}

pixy --stats pi fst dxy \
  --vcf /scratch/amonc/xipho_revision/fst_dxy_pi/Fst_Dxy_Pi_allsites.vcf.gz \
  --populations /scratch/amonc/xipho_revision/fst_dxy_pi/xiph_pops.txt \
  --n_cores 4 \
  --bed_file /scratch/amonc/xipho_revision/fst_dxy_pi/windows.bed \
  --output_folder "$OUTDIR" \
  --output_prefix pixy_output


# Part of log file below:
2026-06-06 04:27:43,966 pixy.__main__:main:821 [INFO]: [pixy] NOTE: The following chromosomes/scaffolds did not have sufficient data to estimate FST: scaffold_234, scaffold_211, scaffold_248, scaffold_265, scaffold_252, scaffold_418, scaffold_156, scaffold_301, scaffold_426, scaffold_464, scaffold_389, scaffold_394, scaffold_385, scaffold_308, scaffold_383, scaffold_257, scaffold_451, scaffold_442, scaffold_154, scaffold_444, scaffold_311, scaffold_423, scaffold_242, scaffold_330, scaffold_217, scaffold_206, scaffold_190, scaffold_447, scaffold_341, scaffold_236, scaffold_298, scaffold_428, scaffold_281, scaffold_323, scaffold_244, scaffold_407, scaffold_427, scaffold_388, scaffold_197, scaffold_347, scaffold_334, scaffold_473
2026-06-06 04:27:43,998 pixy.__main__:main:976 [INFO]: [pixy] All calculations complete at 04:27:43 on 2026-06-06
2026-06-06 04:27:43,998 pixy.__main__:main:981 [INFO]: [pixy] Time elapsed: 05:50:39
2026-06-06 04:27:43,998 pixy.__main__:main:982 [INFO]: [pixy] Output files written to /scratch/amonc/xipho_revision/fst_dxy_pi/pixy
2026-06-06 04:27:43,998 pixy.__main__:main:991 [INFO]: [pixy] If you use pixy in your research, please cite the following paper: Korunes, KL and K Samuk. pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Mol Ecol Resour. 2021 Jan 16. doi: 10.1111/1755-0998.13326.
