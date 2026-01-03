# 1/31/25 (created the out.idepth and out.imiss on 11/12/24)
# Alrighty!!! Filtering on the new and improved VCF dataset from snpArcher (with new Xipho reference)
# Goal create vcf dataset for cluster analyses (sNMF or PopCluster)

# First I want to check the depth and missingness across individuals within the dataset
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


# bcftools get uniq chromosome names; 272 scaffolds, looks good
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz | uniq

# number of ind = 33
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz | wc -l

# Ok, everything checks out so far
# Filtering similarly to the ReLERNN filtering (although note that I am filtering flags explicitly--the mask used in ReLERNN is specific to the ARGweaver vcf, which doesn't include scaffolds <110 kb in length)
# I am removing out the GapRepeatIntervals with a bed mask (one component of ARGweaver mask)
# I am removing the GATK flagged sites with specific VCFtools flags
# Masking out the 
# I am not masking out the the regions that ReLERNN could not get recomb rates for (fewer than 50 sites). Relernn estimates based on just Tapajos VCF, so I feel it is too restrictive for a general Fst scan.

#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N Cluster_VCF_filtering_xipho

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
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.75 \
    --recode \
    --out cluster_vcf_pt1
# After filtering, kept 10218793 out of a possible 46675869 Sites
# SNP archer does output SNPs (via BWA and GATK) for soft-masked regions in the reference


# Thin down to one snp per 10kb
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N Cluster_VCF_filt_final_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_pt1.recode.vcf \
    --thin 10000 \
    --recode \
    --out cluster_vcf_final
# After filtering, kept 99298 out of a possible 10218793 Sites


# Ok, lets assess final VCF file for cluster analysis
# bcftools get uniq chromosome names; 169 scaffolds, so definitely lost some but that makes sense (272 scaffolds to start with)
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf | uniq
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf | uniq | wc -l

# SNP count: 99298
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf|wc -l

# explore filter flags
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf | head -10000


# Get final depth and missingness per individual
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset2

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/cluster_vcf_final.recode.vcf \
    --missing-indv

# Downloading popcluster (first had to enter my information on webpage and agree to terms and conditions)
wget https://cms.zsl.org/sites/default/files/2024-09/PopClusterLnx_25_09_2024.zip

# Unzip
unzip PopClusterLnx_25_09_2024.zip

# Make the popcluster file executable
chmod u+x ./Bin/PopClusterLnx

# Test run, All Good!!
./Bin/PopClusterLnx INP:./Example/ant377NoScale.PcPjt

# renamed parameter file
mv ant377NoScale.PcPjt xipho.PcPjt

# Copy cluster vcf to PopCluster project folder
cp cluster_vcf_final.recode.vcf /scratch/a_monc/postdoc/xipho_project/PopCluster/popc_xipho

# sent my VCF file to Heru and then converted to the dataframe=0 (one ind per line) option using the windows GUI

# Running PopCluster
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q checkpt
#PBS -N PopCluster_xipho

cd /scratch/a_monc/postdoc/xipho_project/PopCluster/popc_xipho

/ddnA/work/a_monc/postdoc/xipho_project/PopCluster/Bin/PopClusterLnx INP:/ddnA/work/a_monc/postdoc/xipho_project/PopCluster/popc_xipho/xipho.PcPjt

# Fantastic, this worked really well!!
# Took like an hour to run
# The .K file has the good stuff

Method      Best_K
DLK2           3
FST/FIS        1


# Best run at K=3 (given by .K file):
/ddnA/work/a_monc/postdoc/xipho_project/PopCluster/popc_xipho/xipho_popc_K_3_R_1

# Added coordinates (lat and long columns) to the Q file for K3, Run 1

# Now opening in QGIS



