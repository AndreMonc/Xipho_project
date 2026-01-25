# Goal:
## Get allele stats for both the ABBA and BABA configurations of the populations
## Following pipeline found here: https://github.com/AndreMonc/allele_stats

# Additional notes here

## Will use my Fst VCF file:
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf`

#### Get SNP count using grep (27114000)
`grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf|wc -l`

## Explanation for which populations are A and B. Figure S1.
First, I want Belem to be population A (fixed for ancestral allele) and Tapajos to be population C (sister to A+B and fixed for alternative allele)
Then, I allow Xingu (popB) to vary in genotype at these sites I expect Xingu individuals to show more alternate alleles nearer to Tapajos, and that there 
will be more overall sites that match this ABBA pattern than match the BABA pattern.


## check line number for header
`grep -n "#CHROM" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf`

= line 322

#### (using 321 for skipRows flag)

## Create the windows over which to calculate allele stats
#### Ok, first, I want to use the .fai for the pseudochromosome reference (new Xiphorhynchus elegans pacbio genome aligned to chromosome-level Chiroxiphia genome)
`/ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa.fai`
## Next, create a genome file map
`awk -v OFS='\t' {'print $1,$2'} /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa.fai > genome_file.txt`
## Then create the windows with bedtoos (non-overlapping 10kb windows)
`bedtools makewindows -g genome_file.txt -w 10000 > windows.bed`

## Get the python file
`wget https://github.com/AndreMonc/allele_stats/blob/main/allele_stats.py`


## allele stats Run1 (Belem as popA--fixed for ancestral allele)
```#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q bigmem
#PBS -N allele_stats_Bel_popA

cd /scratch/a_monc/postdoc/xipho_project/allele_stats

python allele_stats.py --vcfFile /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf --skipRows 321 \
--windowFile windows.bed --popKey popKey.txt --popA belem \
--popB xingu --popC tapajos
```

###
Ok, I got an error
IndexError: list index out of range

And I think it may be due to the presence of some genotypes without two alleles??

## Well, I ran the following command to get all the unique genotypes from the vcf file (Fst_vcf_final.recode.vcf):
`/scratch/a_monc/postdoc/bcftools-1.3/bcftools view /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf | /scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f'[%GT\n]' | sort -u`
## Result:
./.
0/0
0|0
0/1
0|1
1/1
1|1
### So, there are pipes that I need to account for. I made updates to three sections of the code to deal with pipes. Fixed now.


## Flip analysis. Just to compare.
## Explanation for which populations are A and B. Figure S2.
For this analysis I want Xingu to be population A (fixed for ancestral allele) and Tapajos to be population C (sister to A+B and fixed for alternative allele)
Then, I allow Belem to vary in genotype at these sites (I expect Belem individuals to show few alternative alleles and no spatial pattern--no gene flow with Tapajos)

## allele stats Run2 (Belem as popB--allowed to vary in genotypes)
```#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q bigmem
#PBS -N allele_stats_Bel_popB

cd /scratch/a_monc/postdoc/xipho_project/allele_stats

python allele_stats.py --vcfFile /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf --skipRows 321 \
--windowFile windows.bed --popKey popKey.txt --popA xingu \
--popB belem --popC tapajos
```