# 12/5/24 (created the out.idepth and out.imiss on 11/12/24)
## Alrighty!!! Filtering on the new and improved VCF dataset from snpArcher (with new Xipho reference)

## First I want to check the depth and missingness across individuals within the dataset
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=24:00:00
#PBS -q checkpt
#PBS -N xipho_depth_missingness

cd /scratch/a_monc/postdoc/xipho_project/vcftools_output

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --depth
    --out FINAL_pseudo_XIPHO_raw.vcf

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --missing-indv
    --out FINAL_pseudo_XIPHO_raw.vcf
```

## bcftools get uniq chromosome names; 272 scaffolds, looks good
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz | uniq
```
## number of ind = 33
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz | wc -l
```
## Ok, everything checks out so far
## Filtering similarly to the ReLERNN filtering (although note that I am filtering flags explicitly--the mask used in ReLERNN is specific to the ARGweaver vcf, which doesn't include scaffolds <110 kb in length)
#### I am removing the GapRepeatIntervals with a bed mask (one component of ARGweaver mask)
#### I am removing the GATK flagged sites with specific VCFtools flags
#### I am not masking out the the regions that ReLERNN could not get recomb rates for (fewer than 50 sites). Relernn estimates based on just Tapajos VCF, so I feel it is too restrictive for a general Fst scan.

```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N Fst_VCF_filtering_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --exclude-bed /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_GapRepeatIntervals_final.bed \
    --remove-filtered FS_SOR_filter \
    --remove-filtered MQ_filter \
    --remove-filtered RPRS_filter \
    --remove-indels \
    --minQ 30 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.5 \
    --recode \
    --out Fst_vcf_final
# After filtering, kept 27114000 out of a possible 46675869 Sites
# SNP archer does output SNPs (via BWA and GATK) for soft-masked regions in the reference
```

## Conducting a test to see if the GapRepeat bed makes any difference (if snpArcher includes snps from soft-masked regions, then it will make a difference)
#### If snpArcher does not include snps from soft-masked regions, then it won't make a difference
#### So note, below I am not including the --exclude-bed flag
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N Fst_VCF_filtering_test_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --remove-filtered FS_SOR_filter \
    --remove-filtered MQ_filter \
    --remove-filtered RPRS_filter \
    --remove-indels \
    --minQ 30 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.5 \
    --recode \
    --out Fst_vcf_final_test
```
#### After filtering, kept 30136268 out of a possible 46675869 Sites


## Ok, lets assess final VCF file for Fst
## bcftools get uniq chromosome names; 186 scaffolds, so definitely lost some (272 to start with), but makes sense
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf | uniq
```
## SNP count: 27114000
```
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf|wc -l
```
## explore filter flags
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf | head -10000
```
## Get final depth and missingness per individual
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset3

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf \
    --missing-indv

```

## Create appropriate python environment for running simon martin scripts
## https://github.com/simonhmartin/genomics_general/issues/124
```
conda create --name genomics_general python=3.9
conda activate genomics_general
pip install numpy==1.26.0
```

## Following Simon Martin genomics general page
## convert my VCF file to geno file
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=2:00:00
#PBS -q checkpt
#PBS -N Fst_VCF_2_geno

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

python /scratch/a_monc/postdoc/genomics_general/VCF_processing/parseVCF.py -i /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf | bgzip > Fst_VCF_final.geno.gz
```

## Check # of individuals in VCF: 33, good
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf | wc -l
```
## Check individuals in VCF
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_vcf_final.recode.vcf
```
## run popgenWindows.py to get window data
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=2:00:00
#PBS -q checkpt
#PBS -N Fst_stats

source activate genomics_general

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

python /scratch/a_monc/postdoc/genomics_general/popgenWindows.py -w 10000 -g Fst_VCF_final.geno.gz -o Fst_VCF_final.csv.gz -f phased -T 4 -p Bel -p Tap -p Xin --popsFile xiph_pops.txt --writeFailedWindows
```
#### Got the warning messages, but no np.NaN errors, so this looks good. Also, got 110166 Fst windows.
#### Good to go!

## Final Fst file:
```
/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/Fst_VCF_final.csv.gz
```