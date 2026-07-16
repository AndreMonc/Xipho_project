# Goal:
## Get allele stats for both the ABBA and BABA configurations of the populations
## Following pipeline found here: https://github.com/AndreMonc/allele_stats

# Additional notes here

## Will use this VCF file:
`/scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf.gz`

# Uncompress VCF file
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=allele_stats_unzip_%J_stdout.txt
#SBATCH --error=allele_stats_unzip_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=allele_stats_unzip
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################

gunzip -c allele_stats.vcf.gz > allele_stats.vcf
```


## Explanation for which populations are A and B. Figure S1.
First, I want Belem to be population A (fixed for ancestral allele) and Tapajos to be population C (sister to A+B and fixed for alternative allele)
Then, I allow Xingu (popB) to vary in genotype at these sites I expect Xingu individuals to show more alternate alleles nearer to Tapajos, and that there 
will be more overall sites that match this ABBA pattern than match the BABA pattern.


## check line number for header
```
grep -n "#CHROM" /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf
```

= line 527

#### (using 526 for skipRows flag)

## Create the windows over which to calculate allele stats
#### Ok, first, I want to use the .fai for the reference (new Xiphorhynchus elegans pacbio genome)
```
cp /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/snpArcher_wd/results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.fai /scratch/amonc/xipho_revision/allele_stats
```
/scratch/amonc/xipho_revision/allele_stats/xiph_elegans_ref.fa.fna.fai

## Next, create a genome file map
`awk -v OFS='\t' {'print $1,$2'} /scratch/amonc/xipho_revision/allele_stats/xiph_elegans_ref.fa.fna.fai > genome_file.txt`
## Then create the windows with bedtools (non-overlapping 10kb windows)
`bedtools makewindows -g genome_file.txt -w 10000 > windows.bed`

## Get the python file
`wget https://github.com/AndreMonc/allele_stats/blob/main/allele_stats.py`

# Create python environment
conda create -n allele_stats_env \
    python=3.11.4 \
    pandas=2.1.4

conda activate allele_stats_env

pip install numpy==1.26.1


## allele stats Run1 (Belem as popA--fixed for ancestral allele)

```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=175G
#SBATCH --output=allele_stats_Bel_popA_%J_stdout.txt
#SBATCH --error=allele_stats_Bel_popA_%J_stderr.txt
#SBATCH --time=01:00:00
#SBATCH --job-name=allele_stats_Bel_popA
#SBATCH --chdir=/scratch/amonc/xipho_revision/allele_stats/bel_popA
#
#################################################
source ~/miniforge3/etc/profile.d/conda.sh
conda activate allele_stats_env

python allele_stats.py --vcfFile /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf --skipRows 526 \
--windowFile windows.bed --popKey popKey.txt --popA belem \
--popB xingu --popC tapajos
```

## Flip analysis. Just to compare.
## Explanation for which populations are A and B.
For this analysis I want Xingu to be population A (fixed for ancestral allele) and Tapajos to be population C (sister to A+B and fixed for alternative allele)
Then, I allow Belem to vary in genotype at these sites (I expect Belem individuals to show few alternative alleles and no spatial pattern--no gene flow with Tapajos)

## allele stats Run2 (Belem as popB--allowed to vary in genotypes)
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=175G
#SBATCH --output=allele_stats_Xin_popA_%J_stdout.txt
#SBATCH --error=allele_stats_Xin_popA_%J_stderr.txt
#SBATCH --time=01:00:00
#SBATCH --job-name=allele_stats_Xin_popA
#SBATCH --chdir=/scratch/amonc/xipho_revision/allele_stats/xin_popA
#
#################################################
conda activate allele_stats_env

python allele_stats.py --vcfFile /scratch/amonc/xipho_revision/vcf_filtering/allele_stats.vcf --skipRows 526 \
--windowFile windows.bed --popKey popKey.txt --popA xingu \
--popB belem --popC tapajos
```