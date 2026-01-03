# 1/31/25 (created the out.idepth and out.imiss on 11/12/24)
# Alrighty!!! Filtering on the new and improved VCF dataset from snpArcher (with new Xipho reference)
# Goal is to get a good VCF for running IQ-TREE
# Main difference from Fst filtering is to thin by 100 bp (as done in tutorial below):
https://github.com/mmatschiner/tutorials/blob/master/species_tree_inference_with_snp_data/README.md
# Also, max-missing of 75%

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
# I am not masking out the the regions that ReLERNN could not get recomb rates for (fewer than 50 sites).
# Using stricter filtering for max-missing than with Fst

#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=8:00:00
#PBS -q checkpt
#PBS -N IQTREE_filtering_xipho_pt1

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
    --out IQTREE_filt_pt1
# After filtering, kept 25141922 out of a possible 46675869 Sites


# Ok, just thinning the VCF file
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N IQTREE_VCF_final_filtering_xipho

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /scratch/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_filt_pt1.recode.vcf \
    --thin 100 \
    --recode \
    --out IQTREE_vcf_final
# After filtering, kept 6242812 out of a possible 25141922 Sites


# IQ-TREE final vcf:
/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf


# Get final depth and missingness per individual
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -q single
#PBS -N xipho_depth_missingness_dataset1

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf \
    --missing-indv


# SNP count: 6242812
grep -v "^#" /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf | wc -l

# explore filter flags
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -f '%CHROM %POS %REF %FILTER %ALT\n' /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf | head -10000

# Installing IQ-TREE onto Super-Mike 3
# Install 64-bit Linux Intel
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz

tar -xvzf iqtree-2.3.6-Linux-intel.tar.gz

# IQTREE executable:
/scratch/a_monc/postdoc/iqtree-2.3.6-Linux-intel/bin/iqtree2


# to cite
S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nat. Methods, 14:587–589. DOI: 10.1038/nmeth.4285

# Command to run
iqtree -s SNP_data.phy -m MFP+ASC -B 1000 -T AUTO

# Download vcf2phylip
wget https://github.com/edgardomortiz/vcf2phylip/archive/refs/tags/v2.8.zip
unzip v2.8.zip

# Run vcf2phylip
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -q checkpt
#PBS -N vcf2phylip_xipho

python /scratch/a_monc/postdoc/vcf2phylip-2.8/vcf2phylip.py -i /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/IQTREE_vcf_final.recode.vcf --output-folder /scratch/a_monc/postdoc/xipho_project/vcf2phylip --output-prefix xipho -o Tap_xsPIME217_Mus


# Citation for vcf2phylip:
Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861

# phylip file:
/ddnA/work/a_monc/postdoc/xipho_project/vcf2phylip/xipho.min4.phy
cp /ddnA/work/a_monc/postdoc/xipho_project/vcf2phylip/xipho.min4.phy /scratch/a_monc/postdoc/xipho_project/IQTREE

# Actual IQ-TREE run
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q bigmem
#PBS -N IQTREE_xipho

cd /scratch/a_monc/postdoc/xipho_project/IQTREE

/scratch/a_monc/postdoc/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy -st DNA -m MFP+ASC -B 1000 -T AUTO

# Got the following error output in the "e" file, apparently I still had some invariant sites in the file according to IQTREE rules:
ERROR: For your convenience alignment with variable sites printed to /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy
ERROR: Invalid use of +ASC because of 4629913 invariant sites in the alignment

# IQ-TREE automatically created a variable sites only phylip:
/ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy

# Deleting all other IQ-TREE output files in IQTREE folder, except the varsites file
rm xipho.min4.phy.model.gz
rm xipho.min4.phy.log

# New IQ-TREE run with the truly variant site Phylip file
# Initially ran as checkpt, but ran out of memory after 29 minutes. So, rerunning here on bigmem.
# best number of threads?
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q bigmem
#PBS -N IQTREE_xipho_final

cd /scratch/a_monc/postdoc/xipho_project/IQTREE

/scratch/a_monc/postdoc/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy -st DNA -m MFP+ASC -B 1000 -T AUTO


#########################
# IQ-TREE ran successfully!! Yay!!
# Here are some run stats:
Best number of threads = 17
Total number of iterations: 256
CPU time used for tree search: 64888.870 sec (18h:1m:28s)
Wall-clock time used for tree search: 4295.170 sec (1h:11m:35s)
Total CPU time used: 93229.235 sec (25h:53m:49s)
Total wall-clock time used: 6291.955 sec (1h:44m:51s)

977701 parsimony-informative, 635198 singleton sites, 0 constant sites
Bayesian Information Criterion:         GTR+F+ASC+R2
Best-fit model: GTR+F+ASC+R2 chosen according to BIC

Consensus tree written to /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.contree
Reading input trees file /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.contree
Log-likelihood of consensus tree: -17996317.386

Analysis results written to: 
  IQ-TREE report:                /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.iqtree
  Maximum-likelihood tree:       /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.treefile
  Likelihood distances:          /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.splits.nex
  Consensus tree:                /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.contree
  Screen log file:               /ddnA/work/a_monc/postdoc/xipho_project/IQTREE/xipho.min4.phy.varsites.phy.log


  # Using R script with phytools
  # Collapsed nodes with <50 ultra bootstrap support

  # Readjusted location of bootstrap support values in illustrator
