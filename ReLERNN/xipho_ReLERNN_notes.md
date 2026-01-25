# ReLERNN
#### Okay, revisiting these directions on 22 November 2024, more than a year after running!
##### Let's see if we can't get this up and running efficiently!

## created conda environment with older python. Ok starting here on 11/27/24 to install on SuperMike3 (smic is down over Thanksgiving break)
I think I have my order mixed up. I think correct should be create the conda environment, then install ReLERNN within that.
```
conda create --name ReLERNN python==3.8

conda activate ReLERNN
```

## Reinstalled on 17 August 2023
Installing while in ReLERNN environment
```
cd /scratch/a_monc/postdoc/xipho_project/ReLERNN/software
git clone https://github.com/kr-colab/ReLERNN.git
cd ReLERNN
pip install .
```
## After installing ReLERNN with many python versions, 3.8 seemed to work with no issues!!
pip install for hp5y and tensorflow separately
So, 
```
pip install h5py
pip install tensorflow
```
## Run ReLERNN example file
```
cd examples
./example_pipeline.sh
```

#### 2/3/22 and basically from scratch on 5/11/23
#### 5/11/23 coming back to this in order to get recombination rates per chromosome for xipho and phleg
#### Coming back on 22 November 24 to get recombination rates for just xipho (new reference used)

## setup of conda environment (ignoring this 11/22/24 and any CUDA/GPU stuff)
```
conda create --name tf python=3.9
conda activate tf
```
### Just going to use the ReLERNN env with python 3.8
### Need to filter my raw VCF file from snpArcher

*Only use Tapajos individuals*
*filter into a Z version* and an autosomes version

# Final ReLERNN VCFs
## Autosomes
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf
```
## Z chromosome
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf
```

## Created genome (bed) files. For autosome genome file only including scaffolds with > 250 SNPs (omits 7 scaffolds)
## See xipho_VCF_filering_new_for_ReLERNN.md file for more details on VCF filtering and genome file creation
## Autosomes
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.bed
```
## Z chromosome
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.bed
```

## Ok, setting up the all chromes minus Z run.
## Running on SuperMike
```
/scratch/a_monc/postdoc/ReLERNN/ReLERNN
```
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q workq
#PBS -N ReLERNN_xiphoTapajos_autosomes

cd /scratch/a_monc/postdoc/xipho_project/ReLERNN/autosomes

source activate ReLERNN

SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
MU="6.9e-9"
GENTIME="3"
DIR="/scratch/a_monc/postdoc/xipho_project/ReLERNN/autosomes/"
VCF="/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf"
GENOME="/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.bed"

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
    --projectDir ${DIR} \
    --minSites 200
```


# Z only run--Tapajos pops
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=48:00:00
#PBS -q workq
#PBS -N ReLERNN_xiphoTapajos_Zonly

cd /scratch/a_monc/postdoc/xipho_project/ReLERNN/Z_chrom

source activate ReLERNN

SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
MU="6.9e-9"
GENTIME="3"
DIR="/scratch/a_monc/postdoc/xipho_project/ReLERNN/Z_chrom/"
VCF="/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf"
GENOME="/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.bed"

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



## Making 10-kb window file from the 49 scaffolds for which I estimated ReLERNN recombination rates
### Manually created the two-column genome .txt file with scaffold and scaffold length
```
bedtools makewindows -g /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/ReLERNN_49scaffolds_genome_file.txt -w 10000 > 10kbwindows_for_recomb_rates.bed
```

#### Perfect
#### Now run my window_adjust.py program to get the weighted averages for the ReLERNN windows
#### Running locally on cluster, supermike3, 28 Nov 2024
#### Took 18 min to run the below qsub file
#### Input for ARGweaver

```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=6:00:00
#PBS -q workq
#PBS -N window_adjust_all

cd /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN

python window_adjust_allwindows.py --dataFile ReLERNN_Tapajos_auto_and_Z.txt --windowFile 10kbwindows_for_recomb_rates.bed
```

## Now running again with the 200 site cutoff (which will help me know what regions have most reliable recomb rate estimates)
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=1:00:00
#PBS -q workq
#PBS -N window_adjust_200siteCutoff

cd /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN

python window_adjust_200siteCutoff.py --dataFile ReLERNN_Tapajos_auto_and_Z.txt --windowFile 10kbwindows_for_recomb_rates.bed
```

### Ok, so I see now (11/29/2024) that the ReLERNN output is not complete (due to default --minSites 50)
#### Thus, I think I need to mask the ARGweaver input VCF for the ReLERNN windows with <50 sites (basically, those missing from the ReLERNN output but within the range of the ReLERNN genome file)
#### So, I need to take the bedfile difference between this file (the 49 scaffolds used in ReLERNN/ARGweaver):
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_bed_ReLERNN_scaffolds.bed
```
## and this file of the good >= 50-site ReLERNN windows:
```
/ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/ReLERNN_Tapajos_auto_and_Z.bed
```
## bedtools command
```
bedtools subtract -a /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_bed_ReLERNN_scaffolds.bed -b /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/ReLERNN_Tapajos_auto_and_Z.bed > ReLERNN_min50site_mask.bed
```
## Sweet, finally have the ReLERNN Min50Sites mask 
```
/ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/ReLERNN_min50site_mask.bed
```
#### Now, I need to merge this with the main ARGweaver mask
#### Then, I need to create the official recombmap bed file

