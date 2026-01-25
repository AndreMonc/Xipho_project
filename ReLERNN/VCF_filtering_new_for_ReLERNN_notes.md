# 11/22/24
# VCF filtering of Xipho for ReLERNN

## JUST the Tapajos population
#### Remove all non-Tapajos populations from dataset 33-->X individuals.

## Check individuals in raw VCF, get list of Tap indivs
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz

Tap_xsA08267_Mus
Tap_xsBR163-028_Mus
Tap_xsBR163-145_Mus
Tap_xsBR163-212_Mus
Tap_xsMPDS1217_Mus
Tap_xsMPDS1294_Mus
Tap_xsMSF111_Mus
Tap_xsPIME217_Mus
Tap_xsSER013_Mus
Tap_xsTM005_Mus
```
## Add above individuals to a xipho_tapajos.txt file

## Using the ARGweaver vcf file, which has the small scaffolds (<110 kb>) removed. That way the recombmap will correspond fully with the ARG estimate regions
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=24:00:00
#PBS -q checkpt
#PBS -N ReLERNN_Tapajos

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver.recode.vcf \
    --keep /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/xipho_tapajos.txt \
    --recode \
    --out ReLERNN_Tapajos

# bcftools get uniq chromosome names; 84 scaffolds, looks good
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver.recode.vcf | uniq
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos.recode.vcf | uniq | wc -l
```

## explore filter flags
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' ReLERNN_Tapajos.recode.vcf  | head -10000
```
## number of ind
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l ReLERNN_Tapajos.recode.vcf | wc -l
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l ReLERNN_Tapajos.recode.vcf
```
#### Ok, everything checks out so far
## Filter to a autosome dataset and a Z chrom dataset without bad regions
#### This time around adding a minQ filter too!
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N ReLERNN_Tapajos_autosomes_pt1

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos.recode.vcf \
    --exclude-bed /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_Argweaver_final_mask.bed \
    --not-chr Chromosome_Z_RagTag \
    --remove-indels \
    --minQ 30 \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --recode \
    --out ReLERNN_Tapajos_autosomes_hiqual_pt1
```
kept 16888455 out of a possible 40728547 Sites


## Second part (site filtering based on missingness)
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N ReLERNN_Tapajos_autosomes_pt2

# filtering sites by missing data; separating just to make sure of order of operations (so that max-missing filter takes into account the genotype filters)
vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_pt1.recode.vcf \
    --max-missing 0.5 \
    --recode \
    --out ReLERNN_Tapajos_autosomes_hiqual_final
# After filtering, kept 16751683 out of a possible 16888455 Sites
```

## bcftools get uniq chromosome names; only 55 scaffolds, interesting. SNP count matches that of VCFtools, good.
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos_autosomes_hiqual_pt1.recode.vcf | uniq
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf | uniq
grep -v "^#" ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf|wc -l
```
## Check filter flags
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%FILTER\n' ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf | uniq > uniq.poorGATK.test.txt
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%FILTER\n' xipho_poorGATK.recode.vcf | uniq > uniq.poorGATK.txt
sort -u uniq.poorGATK.test.txt > truly.uniq.poorGATK.test.txt
```

## Ok, filter to Z chrom dataset
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N ReLERNN_Tapajos_Z_pt1

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos.recode.vcf \
    --exclude-bed /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_Argweaver_final_mask.bed \
    --chr Chromosome_Z_RagTag \
    --remove-indels \
    --minQ 30 \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --maxDP 50 \
    --minGQ 20 \
    --recode \
    --out ReLERNN_Tapajos_Z_hiqual_pt1
# After filtering, kept 315591 out of a possible 40728547 Sites
```

## bcftools get uniq chromosome names--only Z. Good. SNP count with grep matches VCFtools count.
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos_Z_hiqual_pt1.recode.vcf | uniq
grep -v "^#" ReLERNN_Tapajos_Z_hiqual_pt1.recode.vcf|wc -l
```

## Filtering out sites based on missingness
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N ReLERNN_Tapajos_Z_pt2

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_pt1.recode.vcf \
    --max-missing 0.5 \
    --recode \
    --out ReLERNN_Tapajos_Z_hiqual_final
# After filtering, kept 241074 out of a possible 315591 Sites
grep -v "^#" ReLERNN_Tapajos_Z_hiqual_final.recode.vcf|wc -l
```

## bcftools get uniq chromosome names. 55 chromosome in autosome file and 1 chromosome in Z file.
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf | uniq
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' ReLERNN_Tapajos_Z_hiqual_final.recode.vcf | uniq
```
## explore filter flags
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf  | head -10000
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' ReLERNN_Tapajos_Z_hiqual_final.recode.vcf  | head -10000
```
## number of ind = 10
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf | wc -l
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l ReLERNN_Tapajos_Z_hiqual_final.recode.vcf | wc -l
```

### Final ReLERNN VCFs
### Autosomes
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf
```
### Z chromosome
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf
```
## Get final depth and missingness per individual--autosome VCF
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset5

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf \
    --missing-indv
```

## Get final depth and missingness per individual--Z chrom VCF
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset6

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf \
    --missing-indv
```

## Actual I realized I should check for number of SNPs in each VCF chromosome
## Count # of snps per chromosome in VCF autosome file
```
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.recode.vcf | cut -f 1 | sort | uniq -c
```
## Scaffolds with less than 250 snps (will not include in genome file)
```
scaffold_195
scaffold_120
scaffold_160
scaffold_200
scaffold_210
scaffold_212
scaffold_136
```
## Count # of snps per chromosome in VCF Z chrom file (241,074 SNPs > 250 SNPs)
```
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.recode.vcf | cut -f 1 | sort | uniq -c
```

### Final ReLERNN Genome (bed) files
## Autosomes, only 48 autosomes, so that more than 250 SNPs per scaffold
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_autosomes_hiqual_final.bed
```
## Z chromosome
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ReLERNN_Tapajos_Z_hiqual_final.bed
```

