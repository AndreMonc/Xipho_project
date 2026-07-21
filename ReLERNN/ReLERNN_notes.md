# ReLERNN, for recombination rates
#### Xiphorhyncus manuscript revision
#### 12 April 2026

## create conda environment with older python version
```
conda create --name ReLERNN python==3.8
conda activate ReLERNN
```

## Installed on 12 April 2026
Installing while in ReLERNN environment
```
wget https://github.com/kr-colab/ReLERNN/archive/refs/tags/v1.0.0.tar.gz
tar -xzvf v1.0.0.tar.gz
cd ReLERNN
pip install .
pip install h5py
pip install tensorflow
```

## Run ReLERNN example file
```
cd examples
./example_pipeline.sh
```

# Final ReLERNN VCF
```
/scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf.gz
gunzip -c RelERNN_biallelic_snps_tapajos.vcf.gz > RelERNN_biallelic_snps_tapajos.vcf

```

## Created genome (.bed) file
```
/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed
```

# Running on OU HPC: schooner
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=relernn_%J_stdout.txt
#SBATCH --error=relernn_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=relernn
#SBATCH --chdir=/scratch/amonc/xipho_revision/relernn
#
#################################################

# Initialize conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activate your environment
conda activate ReLERNN

SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
MU="4.6e-9"
GENTIME="3"
DIR="/scratch/amonc/xipho_revision/relernn"
VCF="/scratch/amonc/xipho_revision/vcf_filtering/RelERNN_biallelic_snps_tapajos.vcf"
GENOME="/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed"

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedGenTime ${GENTIME} \
    --assumedMu ${MU} \
    --nTrain 13000 \
    --nVali 2000 \
    --nTest 100

# Train network
${TRAIN} \
    --projectDir ${DIR} \

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR}
```


12 April 2026
Sent off relernn job
Will get results back and need to create a mask for regions without recombination rate (for input to ARGweaver)

14 April 2026
### ReLERNN output is not complete (due to default --minSites 50)
#### Thus, I think I need to mask the ARGweaver input VCF for the ReLERNN windows with <50 sites (basically, those missing from the ReLERNN output but within the range of the ReLERNN genome file)
#### So, I need to take the bedfile difference between this file (the 169 ReLERNN scaffolds):
```
/scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed
```
## and this file of the good >= 50-site ReLERNN windows:
```
awk 'BEGIN{OFS="\t"} {$4=""; sub(/\t\t/, "\t"); print}' ReLERNN_data.bed > ReLERNN_clean_data.bed # Getting rid of nsites column
/scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed
```
## bedtools command to identify regions with < 50 sites
```
bedtools subtract -a /scratch/amonc/xipho_revision/vcf_filtering/ReLERNN_scaffolds.bed -b /scratch/amonc/xipho_revision/relernn/ReLERNN_clean_data.bed > ReLERNN_ft50sites.bed
```
## Sweet, finally have the bed file for regions with fewer than 50 sites (No recomb data estimated)
## This will serve as part of the ARGweaver mask
```
/scratch/amonc/xipho_revision/relernn/ReLERNN_ft50sites.bed
```

