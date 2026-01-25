# Create VCF dataset for GADMA analysis

```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N GADMA_VCF_filtering_xipho

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
    --max-missing 0.75 \
    --recode \
    --out GADMA_vcf_pt1
```
After filtering, kept 25141922 out of a possible 46675869 Sites


## Thin down to one snp per 10kb
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N GADMA_VCF_filt_final_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_pt1.recode.vcf \
    --thin 10000 \
    --recode \
    --out GADMA_vcf_final
```

## After filtering, kept 100913 out of a possible 25141922 Sites


## Ok, lets assess final VCF file for GADMA analysis
## bcftools get uniq chromosome names
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf | uniq
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf | uniq | wc -l
```

## SNP count: 100913
```
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf|wc -l
```

## explore filter flags
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf | head -10000
```

## Get final depth and missingness per individual
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N GADMA_VCF_depth_missing

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf \
    --missing-indv
```


## Trying a more strict dataset for demographic modeling (only 5% missingness)

## Create STRICT VCF dataset for GADMA analysis
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N GADMA_VCF_STRICT_xipho

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
    --max-missing 0.95 \
    --recode \
    --out GADMA_STRICT_vcf_pt1
```
#### Lenient (original) VCF pt 1: After filtering, kept 25141922 out of a possible 46675869 Sites
#### Strict VCF pt 1: After filtering, kept 2415232 out of a possible 46675869 Sites


## Thin down to one snp per 10kb
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N GADMA_STRICT_VC_final_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_pt1.recode.vcf \
    --thin 10000 \
    --recode \
    --out GADMA_STRICT_vcf_final
```
#### Lenient (original) final: After filtering, kept 100913 out of a possible 25141922 Sites
#### Strict (current filtering) final: After filtering, kept 88930 out of a possible 2415232 Sites

## Final VCF strict version
```
/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf
```
## Get scaffold count
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf | uniq | wc -l
```

####
## Get depth and missingnes
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N GADMA_STRICT_VCF_depth_missing

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf \
    --missing-indv
```