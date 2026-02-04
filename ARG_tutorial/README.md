# ARGweaver setup tutorial



# Goal
Provide a step-by-step example of what goes into setting up and running ARGweaver. I draw heavily on insights/recommendations from Hubisz and Siepel (2020):

Hubisz, M., & Siepel, A. (2020). Inference of ancestral recombination graphs using ARGweaver. Pp. 231-266 in: Dutheil, J. Y. (ed). Statistical population genomics. Methods in Molecular Biology, 2090. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-0199-0_10

# Notes

- To begin this tutorial, I recommend downloading the ARG_tutorial folder from GitHub and using this folder as your working directory. I've prepared the various bash commands assuming this to be the case.

- I first created a small VCF file to accompany this mini tutorial by filtering down a massive VCF to just two small scaffolds and thinning to a maximum of 1 snp per 10 bp. Here are the two scaffolds in the VCF, along with their lengths:

`Chromosome_33_RagTag,length=2809104`

`scaffold_241,length=100119`

In an actual study, users may only want to estimate the Ancestral Recombination Graph for a small genomic region (still, it would be good to include a buffer on each side of that region, say 50 kb) or across the whole genome, as in our Xiphorhynchus study. Also, you would not want to conduct any site thinning--I just did that to get the VCF small enough to upload onto Github. 

One of the key things to remember--and which requires some careful attention--is that ARGweaver assumes that any sites not in the VCF are invariant. However, there are often many poor quality sites that we don't want to include in our ARGweaver analysis. To address this issue, ARGweaver conveniently allows us to mask poor quality sites/regions so that they are considered unknown rather than invariant. Using an all-sites VCF would also be a great option, as this would explicitly provide information for the invariant sites. Either way, there will still be some subsets of the genome, such as repetitive regions, that are generally helpful to mask.

# Programs to download

- bcftools
- htslib
- vcftools
- bedtools
- bedops
- miniconda3
- fgbio
- Reseqtools
- ARGweaver

Note, you will need to add each of these programs to your PATH or use full paths
to the excecutables to run. Adding the programs to your PATH may look something like the following, depending on your setup:
```
echo 'export PATH="/where/to/install/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile
```
My installation examples below assume a 64-bit Linux system with an x86_64 architecture.

# Installation instructions

- bcftools:
```
wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
tar -xf bcftools-1.22.tar.bz2
cd bcftools-1.22
./configure --prefix=/where/to/install
make 
make install
```
- htslib
```
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xf htslib-1.21.tar.bz2
cd htslib-1.21
./configure --prefix=/where/to/install
make 
make install
```
- vcftools 
```
wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
tar -xf vcftools-0.1.16.tar.gz
cd vcftools-0.1.16
./configure --prefix=/where/to/install
make
make install
```

- bedtools 
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
tar -zxvf bedtools-2.31.1.tar.gz
cd bedtools2
make
```

- bedops
```
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
```

- miniconda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
- fgbio
```
conda install fgbio
```
- Reseqtools
```
wget https://github.com/hewm2008/Reseqtools/archive/refs/tags/v0.25.tar.gz
tar -xvzf v0.25.tar.gz
cd Reseqtools-0.25
tar -xvzf iTools_Code20180520.tar.gz

# note that iTools excecutable found in .../Reseqtools-0.25/iTools_Code/bin
```
- ARGweaver
```
git clone https://github.com/CshlSiepelLab/ARGweaver.git 
cd ARGweaver
make
```
# Preparation of VCF file input
## Step 1. Take "raw" tutorial VCF file and double check what we're working with

count number of individuals (should = 33):

`bcftools query -l xipho_ARGweaver_tutorial.vcf|wc -l`

#count number of snps (should = 89,353):

`grep -v "^#" xipho_ARGweaver_tutorial.vcf|wc -l`

view scaffolds (should = Chromosome_33_RagTag and scaffold_241):

`bcftools query -f '%CHROM\n' xipho_ARGweaver_tutorial.vcf | uniq`

check out length of scaffolds (just peak inside VCF with `less` and look at header, which will also show other scaffold lengths carried over from the original VCF header)

`less xipho_ARGweaver_tutorial.vcf`


## Step 2a. Remove scaffolds below 110 kb in size
I'm choosing 110 kb as a lower limit on scaffold size because I will eventually look at ARG stats over 10-kb windows and I will drop 50-kb on the ends of each scaffold. Thus, 110 kb is the minimum scaffold size that I can use. We'll use a bed file to remove the one scaffold under 110 kb, since this way can easy scale up when we have more scaffolds to remove.

```
vcftools \
    --vcf xipho_ARGweaver_tutorial.vcf \
    --out xipho_ARGweaver_tutorial_above110kb \
    --exclude-bed scaffolds_under_110kb.bed \
    --recode
```

After filtering, kept 89041 out of a possible 89353 Sites

## Step 2b. Optional--removal of other scaffolds based on other criteria (such as a lack of recombination data)
Later, when setting up the command to run the ARGweaver function `arg-sample`, you will need to either provide a single recombination rate or a recombination map with recombination rates across the genome (as a bed file) to ARGweaver. If you estimate a recombination map with, say, ReLERNN (Adrion et al. 2020) it is likely that your map will not contain recombination data for certain smallish scaffolds. In that case you can follow the filtering method of Step 2a again and remove scaffolds with no recombination map data. For the purposes of this tutorial, however, we have recombination data for Chromosome_33_RagTag (`recomb_map.bed`), so we'll keep moving forward.

## Step 3. Zip and index the VCF file

```
bgzip -c xipho_ARGweaver_tutorial_above110kb.recode.vcf > xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz
tabix -p vcf xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz
```

## Step 4. Create bed file for all scaffolds in the VCF file 
(Provided in this tutorial `final_scaffolds.bed`). To create this, first do a final check of what scaffolds are included in the VCF:
```
bcftools query -f '%CHROM\n' xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz | uniq
```

Then, you can do some quick manual work in excel, grabbing the lengths of the appropriate chromosomes/scaffolds from the VCF header, and creating the final scaffolds bed file.

## Step 5. Create the bed map for slicing up the VCF
For our Xiphorhynchus paper we created 2 Mb sliding windows with 100kb overlap, which was quite manageable to process. However, since we only have one small chromosome at this point--long enough for only one 2 MB window--I'm going to instead use 500 kb sliding windows with 100 kb overlap, just to generate a few more windows.

```
bedtools makewindows -b final_scaffolds.bed -w 500000 -s 400000 > sliding_windows.bed
```

## Step 6. Create individual bed files with coordinates for each VCF window we want to create

```
mkdir -p bed_windows

counter=1
cat sliding_windows.bed | while read LINE; do
    name=$(echo $LINE | awk '{print $1}')
    echo "$LINE" > bed_windows/"${name}-${counter}.bed"
    ((counter++))
done
```

This creates a new folder `bed_windows` and outputs 8 bed files, each containing the coordinates for one sliding window.

## Step 7. Slice VCF into chunks for parallel processing
```
mkdir -p vcf_windows

for i in ./bed_windows/*.bed; do
    base=$(basename "$i" .bed)
    bcftools view -R "$i" xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz -Oz \
        -o vcf_windows/"${base}.vcf.gz"
done
```

## Step 8. Index each VCF file
```
for i in ./vcf_windows/*.vcf.gz; do
    tabix -p vcf "$i"
done
```

Ok, excellent!! We've now got our VCF files chunked and ready to go! Now, before we can set up our ARGweaver analysis, we have a few more steps involving masking of low quality sites across the genome.


# Masking 
There are two ways to mask in ARGweaver: individual-based masks and global masks that affect all individuals. For this tutorial, we will focus on the global masks, since this should be sufficient for most projects.

Our basic strategy will be to create bed files to map poor-quality regions for which we can't be confident in the genotype calls. Then we will combine these bed files to create the final mask map, which will tell ARGweaver to treat these regions as unknown.

## Step 9. Create poor-quality sites bed file based on filter flags

First, we need to assess what filter flags we have in our VCF file (these flags are often added during the GATK pipeline and are generally filtered out of datasets before analyses, since they are consider poor quality according to GATK best practices):

```
bcftools query -f '%FILTER\n' xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz | sort -u > uniq.xipho.GATK.filters.txt
```

This gives us the following output:
```
.
FS_SOR_filter
FS_SOR_filter;MQ_filter
FS_SOR_filter;MQ_filter;RPRS_filter
FS_SOR_filter;RPRS_filter
MQ_filter
MQ_filter;RPRS_filter
RPRS_filter
```

The "." is equivalent to a PASS. The other filter flags boil down to the following:
```
FS_SOR_filter
MQ_filter
RPRS_filter 
```

We can then filter our original VCF to just these sites and create a bed file from that VCF (`vcf2bed` is a function within bedops):
```
vcftools \
    --gzvcf xipho_ARGweaver_tutorial_above110kb.recode.vcf.gz \
    --out xipho_poorGATK \
    --keep-filtered FS_SOR_filter  \
    --keep-filtered MQ_filter  \
    --keep-filtered RPRS_filter  \
    --recode
```

After filtering, kept 44414 out of a possible 89041 Sites

```
vcf2bed --max-mem=8G --sort-tmpdir=${PWD} < xipho_poorGATK.recode.vcf > xipho_poorGATK.bed
```

## Step 10. Create a bed file for the N sites and repetitive sites in our reference genome

N sites and repetitive sites in the reference genome could mislead our estimation of ARGs--better to treat these sites as unknown regions. Our reference for the Xiphorhynchus paper is on Dryad (`xipho_elegans_ragtagRef_no_W.fa`), but it is too big to have here on GitHub. So, I will just show this step and include the bed file output here in the tutorial folder.

First, we'll turn soft masked regions (repetive regions identified by Repeat Masker) into hard-masked regions: 
```
fgbio HardMaskFasta --input=xipho_elegans_ragtagRef_no_W.fa --output=xipho_elegans_ragtagRef_no_W_hardMasked.fa
```

Second, we'll get a bed file for all the hard-masked regions of the reference genome:
```
iTools Fatools findN -InPut xipho_elegans_ragtagRef_no_W_hardMasked.fa -OutPut xipho_elegans_ragtagRef_no_W_hardMaskedIntervals.bed
```

Third, we can optionally merge abutting and nearby regions (separated by 100 bp or less) in the bed file:
```
bedtools merge -d 100 -i xipho_elegans_ragtagRef_no_W_hardMaskedIntervals.bed > xipho_merged_GapRepeatIntervals.bed
```

## Step 11. Merge poor-quality and gap/repeat bed files for final bed file mask
```
bedops --merge xipho_poorGATK.bed xipho_merged_GapRepeatIntervals.bed > xipho_Argweaver_mask.bed
```
Great! Now we have a full map of the poor quality GATK-flagged sites and N + repetive sites across the genome.

## Step 12. Let's estimate some ARGs!
We use the ARGweaver function `arg-sample` to actually carry out estimation of the Ancestral Recombination Graph and get output smc.gz files. Here is a full `arg-sample` command that incorporates the various files that we filtered and created during this tutorial. The command below takes the first VCF file window that we created: Chromosome_33_RagTag-1.vcf.gz.

```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -q single
#PBS -N Chromosome_33_RagTag-1

cd /project/gthom/a_monc/xipho_project/ARG_tutorial

./ARGweaver/bin/arg-sample --vcf ./vcf_windows/Chromosome_33_RagTag-1.vcf.gz \
--region Chromosome_33_RagTag:1-500000 \
--vcf-min-qual 30 \
--vcf-genotype-filter "DP<5;DP>50;GQ<20" \
--maskmap xipho_Argweaver_mask.bed \
--mask-cluster 2,5 \
--mutrate 4.6e-9 \
--recombmap recomb_map.bed \
--popsize 391476 \
--compress-seq 5 \
--ntimes 20 \
--maxtime 1e7 \
--delta 0.005 \
--sample-step 50 \
--iters 2000 \
-o Chromosome_33_RagTag-1_out
```

# Congratulations!! 
You have now set up the estimation of an Ancestral Recombination Graph with ARGweaver!
You can use custom scripts to prepare the jobs across all the VCF windows for which you want to estimate ARGs. Ultimately, the exact setup of these job scripts will vary based on your HPC requirements, paths, etc.

For each `arg-sample` run, you will need to update the `--vcf`, `--region`, and `-o` flags to run different VCF windows. If running on an HPC (highly recommended!), you will want to update the job name too for each run. In this tutorial folder, I have included an `ARGweaver_jobinfo` file that I used for the Xiphorhynchus paper, that organizes the key info that differed across each `arg-sample` run. Then, I used a custom python script `ARGjobs_xipho.py` (also in the tutorial folder) to prepare all the jobs. If a job terminates before the desired number of MCMC iterations, you can simply add the flag `--resume` to carry on from where the ARG estimation left off.

Once the `arg-sample` runs finish, you will have a set of smc.gz files that you can plug into the nicely prepared [ARG analysis pipeline](https://github.com/CshlSiepelLab/bird_capuchino_analysis/tree/master/ARG_analysis) used for Sporophila seedeaters (Hejase et al. 2020).





