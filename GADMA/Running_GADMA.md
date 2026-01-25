# Goal:
#### I need to run different demographic models to test whether there was a bottleneck in the Belem population
# 19 January 2026

## GADMA esto no funcionó, ver mas abajo que si funciono
#crear el ambiente
conda create -n gadma_env python=3.10 
conda activate gadma_env 
conda config --add channels bioconda 
conda config --add channels conda-forge

#instalar GADMA
conda install -c bioconda gadma
pip install momi

#verificar la instalacion
gadma --test

si funciona bien deberia salir un mensjae asi:
--Finish pipeline--

--Test passed correctly--
Thank you for using GADMA!


#### no pude instalarlo, genera muchas fallas y al final parece que es mejor cambiar la version de python 
# Below works:
## 1. Crear ambiente con Python 3.9
```
conda create -n gadma_env python=3.9 -y
```
## 2. Activar el ambiente
```
conda activate gadma_env
```
## 3. Instalar dependencias básicas desde conda
```
conda install -c conda-forge -c bioconda numpy=1.24.4 h5py=3.10.0 scikit-allel=1.3.7 scipy matplotlib pandas -y
```

## 4. Instalar GADMA y librerías específicas desde pip
```
pip install gadma==2.0.3 moments==1.0b0 momi
```

## 5. Probar la instalación
```
gadma --test
```

Funcionoooooo

## Successfully installed 
attrs-25.4.0 autograd-1.8.0 cython-3.2.4 dadi-2.4.3 demes-0.2.3 demesdraw-0.4.1 gadma-2.0.3 jsonschema-4.25.1 jsonschema-specifications-2025.9.1 moments-1.0b0 moments-popgen-1.3.1 momi-2.1.21 msprime-1.3.4 newick-1.11.0 nlopt-2.7.1 pandas-2.2.2 pysam-0.23.3 referencing-0.36.2 rpds-py-0.27.1 ruamel.yaml-0.16.12 tskit-0.6.4

## Most importantly, I installed moments, the default demographic inference tool for GADMA
moments-1.0b0


## Pop file
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xiph_pops.txt`



## Running GADMA--one option, but not with full .yaml
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N GADMA_run_test

cd /scratch/a_monc/postdoc/xipho_project/GADMA


gadma --input vcf_file.vcf,popmap_file -o run_test
```

## Calculation of effective sequence length
L = (X - Y) / X * Nseq

Nseq, Total length of sequence (length of the reference genome, xipho_elegans_ragtagRef_no_W.fa) = 1120117357
X, Total number of snps received from this data (my raw SNP count from snpArcher) = 46675869
Y, SNPs filtered out (X - SNPs in GADMA VCF file), (46675869-100913) = 46574956

## Final Sequence calculation
L = (46675869 - 46574956) / 46675869 * 1120117357
L = 100913 / 46675869 * 1120117357
L = 2421688

## Full run of GADMA with .yaml
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N GADMA_run1

cd /scratch/a_monc/postdoc/xipho_project/GADMA

source activate gadma_env

gadma -p param_file_xipho.yaml -o /scratch/a_monc/postdoc/xipho_project/GADMA/run1
```

## Run 2 with fixing of divergence times for both splits
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q checkpt
#PBS -N GADMA_run2

cd /scratch/a_monc/postdoc/xipho_project/GADMA

source activate gadma_env

gadma -p param_file_xipho_run2.yaml -o /scratch/a_monc/postdoc/xipho_project/GADMA/run2
```
## Run 3 with better projection (from easySFS)
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N GADMA_run3

cd /scratch/a_monc/postdoc/xipho_project/GADMA

source activate gadma_env

gadma -p param_file_xipho_run3.yaml -o /scratch/a_monc/postdoc/xipho_project/GADMA/run3
```

## Run 4 created range 1 SD above and below the two split times
#### doing only 1 repeat to improve time
#### droppping wall time to 4 hrs
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N GADMA_run4

cd /scratch/a_monc/postdoc/xipho_project/GADMA

source activate gadma_env

gadma -p param_file_xipho_run4.yaml -o /scratch/a_monc/postdoc/xipho_project/GADMA/run4
```

## Calculation of effective sequence length for stricter VCF (used in Run 5, below)
L = (X - Y) / X * Nseq

Nseq, Total length of sequence (length of the reference genome, xipho_elegans_ragtagRef_no_W.fa) = 1120117357
X, Total number of snps received from this data (my raw SNP count from snpArcher) = 46675869
Y, SNPs filtered out (X - SNPs in GADMA VCF file), (46675869-88930) = 46586939

## Final Sequence calculation
L = (46675869 - 46586939) / 46675869 * 1120117357
L = 88930 / 46675869 * 1120117357
L = 2134123

# Run 5, using a new VCF file (stricter filtering; VCF Dataset 5)
#### Updated projection appropriate for new VCF file
#### Updated sequence length, given fewer SNPs
#### 20 repeats (in yaml)
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=24:00:00
#PBS -q bigmem
#PBS -N GADMA_run5

cd /scratch/a_monc/postdoc/xipho_project/GADMA

source activate gadma_env

gadma -p param_file_xipho_run5.yaml -o /scratch/a_monc/postdoc/xipho_project/GADMA/run5
```
