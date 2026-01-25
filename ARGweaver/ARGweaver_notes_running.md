# ARGweaver begins again!
## Reviving this file on 17 November 2024, to run ARGweaver with new Xiphorhynchus reference VCF
## Will update what I need to

## 4/27/23
## Redownloaded on 6 September 2023
```git clone https://github.com/CshlSiepelLab/ARGweaver.git 
cd ARGweaver
make
```

## Download SAMtools
```mkdir SAMtools
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=/scratch/a_monc/postdoc/SAMtools/
make
make install # This works when I use the prefix option after ./configure!! Yay!
```

## Add SAMtools to bash profile
`/scratch/a_monc/postdoc/SAMtools/bin`

## Download bedops
## /scratch/a_monc/postdoc/bedops/bedops
```mkdir bedops
cd bedops
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
```
Add bedops to path or will not work

## Download PHAST
## Leaving this from my April installation, not sure if I actually use it or not
```mkdir PHAST
wget http://www.netlib.org/clapack/clapack.tgz # Didn't work so had to scp to cluster
tar -xvzf clapack.tgz
cd CLAPACK-3.2.1
cp make.inc.example make.inc && make f2clib && make blaslib && make lib
```
```
git clone https://github.com/CshlSiepelLab/phast.git
cd src
make CLAPACKPATH=/scratch/a_monc/postdoc/PHAST/CLAPACK-3.2.1
```

## Download htslib (But looks like I have bgzip and tabix through conda as default)
## So not actually using installed material below from htslib
```
mkdir htslib
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xf htslib-1.21.tar.bz2
cd htslib-1.21
./configure --prefix=/scratch/a_monc/postdoc/htslib
make 
make install
```

## "ARGweaver assumes that any site which does not appear in the input file is invariant", so I will need to create a bunch of masks
## My VCF file is the raw output VCF from snpArcher (no filtering done, just have flags for the GATK best practices sites):

`/ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz`

## ARGWEAVER analysis
## So, what am I working with in terms of individuals and chromosomes?
```bcftools query -l FINAL_pseudo_XIPHO_raw.vcf.gz #33 individuals
bcftools query -f '%CHROM\n' FINAL_pseudo_XIPHO_raw.vcf.gz | uniq #34 chromosomes (Chiroxiphia as refererence) + 238 xipho scaffolds = 272 unique "blocks" in CHROM field
```

## VCF Filtering
I should just remove all contigs less than 110kb from the VCF. That way ARGweaver will not deal with them at all. These will go anyways during ARG processing.
Add a --mac 1 filter to get rid of all monomorphic sites (shouldn't be any I don't think).
Don't use remove-filtered-all, I will mask GATK-flagged poor-quality sites.
Checked missingness across individuals on raw vcf; result was no outliers--all <6% missingness. Checked mean depth across individuals on raw vcf; no major outliers, lowest was 7.9, highest was 23.2

## Raw input is 46,675,869 SNPs
```#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=24:00:00
#PBS -q checkpt
#PBS -N FiltVCF_for_ARGweaver_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --out xipho_ARGweaver \
    --mac 1 \
    --exclude-bed small_scaffold_110kb_ARGfilt.bed \
    --recode
```

## Check scaffolds output
`/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' xipho_ARGweaver.recode.vcf | uniq | wc -l`
## Output 40,728,547 SNPs and 84 scaffolds--just what I expected, sweet


## Now create VCF with only the 49 scaffolds (autosome + Z chromosome) for which I have recombination data (see ReLERNN for more details)
## (In xipho_ReLERNN_new.md file I have the beds separated into two: one each for autosomes and Z chromosome)
## First create combined bed file from the two ReLERNN bedfiles (manually paste the Z chrom at end of file)
## So here is the 49-scaffold bed file:
`/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_bed_ReLERNN_scaffolds.bed`

# What I actually need is a bed file with the scaffolds to exclude (84 scaffs - 49 ReLERNN scaffs = 35 scaffolds to exclude)
# Here is this "exclude" bed file, with scaffolds to be removed from ARGweaver VCF:
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/nonReLERNN_scaffs_35.bed`

```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=20
#PBS -l walltime=3:00:00
#PBS -q checkpt
#PBS -N final_ARG_vcf_filter

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver.recode.vcf \
    --out xipho_ARGweaver_final \
    --exclude-bed /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/nonReLERNN_scaffs_35.bed \
    --recode
```
## After filtering, kept 40558259 out of a possible 40728547 Sites

## Check output
```bcftools query -l xipho_ARGweaver_final.recode.vcf | wc -l #33 individuals
bcftools query -f '%CHROM\n' xipho_ARGweaver_final.recode.vcf | uniq #49 chromosomes, all looks good
bcftools query -f '%CHROM\n' xipho_ARGweaver_final.recode.vcf | wc -l
```

## VCF which I will then break up into 2MB chunks prior to running ARGweaver:
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf`

## Get final depth and missingness per individual
```#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset4

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf \
    --missing-indv
```

## What is the GATK-flagged poor-quality sites situation on snpArcher output?
## Below I list all filters in the raw vcf file, so I guess "." is apparently equivalent to PASS
```/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%FILTER\n' FINAL_pseudo_XIPHO_raw.vcf.gz | uniq > uniq.xipho.GATK.filters.txt
sort -u uniq.xipho.GATK.filters.txt > truly.uniq.xipho.GATK.filters.txt
```
## Output:
.
FS_SOR_filter
FS_SOR_filter;MQ_filter
FS_SOR_filter;MQ_filter;RPRS_filter
FS_SOR_filter;RPRS_filter
MQ_filter
MQ_filter;RPRS_filter
RPRS_filter

## So, here I want to create a VCF with just the poor quality GATK sites, so can then create a mask for these sites
Weirdly, this VCF commend doesn't provide the expected output with remove-filtered ".", so using keep-filtered for poor quality sites.
This keeps 2,187,454 SNPs out of 40,728,547

```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N create_poorGATK_VCF

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver.recode.vcf \
    --out xipho_poorGATK \
    --keep-filtered FS_SOR_filter  \
    --keep-filtered MQ_filter  \
    --keep-filtered RPRS_filter  \
    --recode
```

## Double-check flags on poorGATK output VCF
```/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%FILTER\n' xipho_poorGATK.recode.vcf | uniq > uniq.poorGATK.txt
sort -u uniq.poorGATK.txt > truly.uniq.poorGATK.txt
```

## Convert this GATK filters VCF to bed file (run interactively on cluster). Then use that bed file as a mask for ARGweaver.
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N vcf2bed_GATK

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

/scratch/a_monc/postdoc/bedops/bin/vcf2bed --max-mem=8G --sort-tmpdir=${PWD} < xipho_poorGATK.recode.vcf > xipho_poorGATK.bed
```

## What about a 50-bp buffer on both ends of the GATK poor sites? Sounds reasonable to me.
`https://bedtools.readthedocs.io/en/latest/content/tools/slop.html`

## old command 
`bedtools slop -b 50 -i xipho_poorGATK.bed -g GATKpoorSites_genomefile.txt > xipho_poorGATK_withSlop.bed`
## new command 
`bedtools slop -b 50 -i xipho_poorGATK.bed -g /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_genome_file.txt > xipho_poorGATK_withSlop.bed`

###########################################################################
## ARGweaver Masks --> masked sites/regions will be considered "unknown" (rather than invariant)
## individual depth masks (see bed_depth_maps_from_bams folder on my laptop for details on how prepared)
`/ddnA/work/a_monc/postdoc/xipho_project/Depth_from_BAMs_500bp_windows/ind_mask_file.txt`
## GATK poor-quality sites
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_poorGATK_withSlop.bed`
## Mask for N/n sites and Repeats (see Get_repeat_mask_bedfile.md) in the reference genomes
`/ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_GapRepeatIntervals_final.bed`
## High missing data (Used in my first ARGweaver estimate) but I don't need this mask! Masking genotypes is sufficient 
## Indels--ARGweaver will just ignore these, so don't need to mask I don't think (I masked in first ARGweaver estimate)
## low coverage, high coverage, SNP quality
`--vcf-geno-type-filter ~DP<5;DP>50;GQ<20~`
## min base quality
`--vcf-min-qual 30`

## Merge the GATK poor-quality sites and GapRepeat beds
## Bedops v2.4.41 to merge Xipho bed files for Argweaver mask
```
bedops --merge /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_poorGATK_withSlop.bed /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_GapRepeatIntervals_final.bed > xipho_Argweaver_final_mask.bed
```

## Not-Final ARGweaver mask (yes, confusing given the file name!). Updated the mask one for time after running ReLERNN (see below)
`/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_Argweaver_final_mask.bed`

## Actually, I realized I needed to also merge in to this bed mask the regions that ReLERNN could not get recomb rates for (fewer than 50 sites)
```
bedops --merge /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_Argweaver_final_mask.bed /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/ReLERNN_min50site_mask.bed > xipho_Argweaver_truly_final_mask.bed
```

## Truly final ARGweaver mask
`/ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/xipho_Argweaver_truly_final_mask.bed`


###############################################################################

## Other ARGweaver input
## Using mutation rate of 2.3e-9 mutations/site/year after Ficedula study (Smeds et al., 2016; Direct estimate of the rate of germline mutation in a bird)
## Using generation time of 3 years, which is the average of estimated values for the genus Xiphorhynchus in Appendix S4 (Bird et al. 2020; Generation lengths of the world’s birds and their implications for extinction risk)
## Need priors for three items: mutation rate, recombination rate, and coalescence (based on population sizes)--RELERNN comes in here
`--mutrate 6.9e-9` #units of expected mutations per base per generation (using Ficedula rate of 2.3e-9 mutations/site/year x generation time of 3 years)
`--recombrate <rate>` #units of the probability of a recombination between two neighboring bases per generation.

## Using ReLERNN to make the recombmap
## Need to filter VCFs for ReLERNN (autosomes and Z)
## I am providing the raw relernn results (combined for autosomes and Z, and not including windows with <50 sites. I also masked these missing low-site windows with the ARGweaver truly final mask)
`--recombmap`

`/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_recomb_map_ReLERNN.bed`

`--popsize`

Ne for Xiphorhynchus: 391476 #Averaged from values in Thom GBE Table S16

## Time discretization
`--ntimes 20` #number of time points, default is 20.
`--maxtime 1e7` #the maximum time point in the model, in units of generations, default is 200,000 generations. Use a time that is greater than divergence within group
`--delta 0.005` #default value of 0.01 (larger values give more resolution of recent time point sampling)--inspect discrete time values upon running. (Sporophila paper used 0.005)

## Sampling frequency
`--sample-step 20` #default of 10

## Site compression
To do: calculate the frequency of non-singleton variants

`--compress-seq 5` #site compression, compression results in a loss of resolution but speeds up processing time. Use conservatively.

## Number of MCMC iterations
```
--resume
--iters 2000 #to add iterations if wanted
```

## 21 November 2024, split the ARGweaver VCF generated earlier into 2MB sliding windows with 100kb overlaps
## Manually created the following 49-scaffold bed file (the ReLERNN scaffolds, necessary for recombmap)
## Next step is to make windows off that bed file (2MB sliding windows with 100kb overlaps):
```
/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_bed_ReLERNN_scaffolds.bed
bedtools makewindows -b /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_bed_ReLERNN_scaffolds.bed -w 2000000 -s 1900000 > ARGweaver_final_windows.bed
```
# Here is the final windows file, that has the windows with recombination rates from ReLERNN:
`/ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_final_windows.bed` #604 windows, solid. Nice improvement from the previous set.

## Make a big series of bed files, each with a single window for a VCF file
```
counter=1; cat /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_final_windows.bed | while read LINE; 
do 
    name=`echo $LINE | awk '{print $1}'` 
    echo "$LINE" > "/scratch/a_monc/postdoc/xipho_project/ARGweaver/bed2MBwindows/$name-$counter.bed"
    ((counter++)); done

ls -lR *.bed | wc -l
ls 
```

## My ARGweaver VCF file needs to be gzipped.
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=2:00:00
#PBS -q checkpt
#PBS -N bgzip_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

/scratch/a_monc/postdoc/htslib-1.21/bgzip -c /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf > xipho_ARGweaver_final.recode.vcf.gz
```

## Index my ARGweaver VCF
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=2:00:00
#PBS -q checkpt
#PBS -N tabix_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

/scratch/a_monc/postdoc/htslib-1.21/tabix -p vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf.gz
```

## Loop over the bed window files in order to make all the vcf window files.
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N bed_to_vcfWindows_xipho

for i in /scratch/a_monc/postdoc/xipho_project/ARGweaver/bed2MBwindows/*.bed
do 
  bcftools view -R $i /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xipho_ARGweaver_final.recode.vcf.gz > $i.vcf
done 
```


## Move VCF window files to their respective folder (26 Nov 2024)
```
mv *.vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows
ls -lR *.vcf | wc -l
```

## rename files .bed.vcf to just .vcf
`find . -name "*.bed.vcf" -exec sh -c 'mv "$1" "${1%.bed.vcf}.vcf"' _ {} \;`


## Once I have all VCF windows, gzip them all
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N bgzip_xiph_windows

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows

for i in /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows/*.vcf
do 
  /scratch/a_monc/postdoc/htslib-1.21/bgzip -c $i > $i.gz  
done 
```

## then index all the VCF windows
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N tabix_xipho_windows

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows

for i in /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows/*.vcf.gz
do 
  /scratch/a_monc/postdoc/htslib-1.21/tabix -p vcf $i
done
```

## Made all my ARGweaver jobs using my custom ARGjobs_xipho.py file
## Output 604 qsub files (without resume flag)
## Started running ARGweaver on December 4
## 83 jobs finished by Dec 9

## This command sends off the first 83 jobs (all Chromosome 1)
```
for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome1-*.qsub; do
    qsub $file
done
```

## Can I submit more than 86 on smic and some stay in queue?? Trying to add Chromosome 10 jobs
```
for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome10-*.qsub; do
    qsub $file
done
```

## Now submitting more jobs on supermike3
```
for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome11-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome12-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome13-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome14-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome15-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome16-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome17-*.qsub; do
    qsub $file
done
```


## After submitting these first jobs 175 jobs, I created the "resume" versions of the ARGweaver qsubs
```
for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome19-*.qsub; do
    qsub $file
done


for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome2-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome3-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome4-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume/Chromosome9-*.qsub; do
    qsub $file
done

for file in /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/qsubs_without_resume_argweaver3/ChromosomeZ-*.qsub; do
    qsub $file
done
```

## Example qsub file (same setup for all windows). If a job did not finish, I resubmitted the job and included the "--resume" flag to the arg-sample command.
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N Chromosome1-1

cd /project/gthom/a_monc/xipho_project/ARG_tutorial

/scratch/a_monc/postdoc/ARGweaver/bin/arg-sample --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows/Chromosome_1_RagTag-1.vcf.gz \
--region Chromosome_1_RagTag:1-2000000 \
--vcf-min-qual 30 \
--vcf-genotype-filter "DP<5;DP>50;GQ<20" \
--ind-maskmap /ddnA/work/a_monc/postdoc/xipho_project/Depth_from_BAMs_500bp_windows/ind_mask_file.txt \
--maskmap /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/xipho_Argweaver_truly_final_mask.bed \
--mask-cluster 2,5 \
--mutrate 6.9e-9 \
--recombmap /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_recomb_map_ReLERNN.bed \
--popsize 391476 \
--compress-seq 5 \
--ntimes 20 \
--maxtime 1e7 \
--delta 0.005 \
--sample-step 50 \
--iters 2000 \
-o Chromosome1-1_out
```

## Ok in total 592 windows finished out of 604 windows
## 12 windows had errors in ARGweaver, all error messages copied to my "keeping track of submissions file"
## So 98.0% windows finished (592/604)

## Now moving forward (March 20, 2025) with processing ARG smc files for the 592 windows. I suspect a few more windows may get dropped during processing. We'll see