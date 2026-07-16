# 4/24/26 
####
### Alrighty!!! Filtering on the new and improved VCF dataset from snpArcher (with Xipho scaffold-level reference)
# Goal

#### Main difference from Fst filtering is to thin by 100 bp (as done in tutorial below):
`https://github.com/mmatschiner/tutorials/blob/master/species_tree_inference_with_snp_data/README.md`
#### Also, max-missing of 75%. See "VCF_filtering_for_IQTREE.md" for full details.

## IQ-TREE final vcf:
```
/scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz
```

## Installing IQ-TREE onto Super-Mike 3
## Install 64-bit Linux Intel
```
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz

tar -xvzf iqtree-2.3.6-Linux-intel.tar.gz
```
## IQTREE executable:
```
/home/amonc/iqtree-2.3.6-Linux-intel/bin/iqtree2
```

## to cite
```
S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nat. Methods, 14:587–589. DOI: 10.1038/nmeth.4285
```

## Command to run
```
iqtree2 -s SNP_data.phy -m MFP+ASC -B 1000 -T AUTO
```

## Download vcf2phylip
```
wget https://github.com/edgardomortiz/vcf2phylip/archive/refs/tags/v2.8.zip
unzip v2.8.zip
```

## Run vcf2phylip
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=vcf2phylip_%J_stdout.txt
#SBATCH --error=vcf2phylip_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=vcf2phylip
#SBATCH --chdir=/scratch/amonc/xipho_revision/iqtree/vcf2phylip
#

python /home/amonc/vcf2phylip-2.8/vcf2phylip.py -i /scratch/amonc/xipho_revision/vcf_filtering/IQTREE.vcf.gz --output-folder /scratch/amonc/xipho_revision/iqtree/vcf2phylip --output-prefix xipho -o XELEGANS_MPEG75162_Mus
```

# Copy raw phylip output to ourdisk:
cp /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/iqtree

## Citation for vcf2phylip:
Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861


## Actual IQ-TREE run
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=iqtree_%J_stdout.txt
#SBATCH --error=iqtree_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=iqtree
#SBATCH --chdir=/scratch/amonc/xipho_revision/iqtree
#

iqtree2 -s /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy -st DNA -m MFP+ASC -B 1000 -T AUTO
```

ENDED HERE on FRIDAY PM 4/24/26, pick up below

#### Got the following error output in the "e" file, apparently I still had some invariant sites in the file according to IQTREE rules:
ERROR: For your convenience alignment with variable sites printed to /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy
ERROR: Invalid use of +ASC because of 1341895 invariant sites in the alignment

## IQ-TREE automatically created a variable sites only phylip:
```
/scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy
```
## Deleting all other IQ-TREE output files in IQTREE folder, except the varsites file
```
rm xipho.min4.phy.model.gz
rm xipho.min4.phy.log
```

## New IQ-TREE v2.3.6 run with the truly variant site Phylip file
### Initially ran as checkpt, but ran out of memory after 29 minutes. So, rerunning here on bigmem.
### best number of threads?
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=iqtree_final_%J_stdout.txt
#SBATCH --error=iqtree_final_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=iqtree_final
#SBATCH --chdir=/scratch/amonc/xipho_revision/iqtree
#

iqtree2 \
  -s /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy \
  -st DNA \
  -m MFP+ASC \
  -B 1000 \
  -T ${SLURM_CPUS_PER_TASK}
```

#########################
### IQ-TREE ran successfully!! Yay!!
## Here are some run stats:
```
Alignment has 34 sequences with 2778634 columns, 2775839 distinct patterns
1733177 parsimony-informative, 1045457 singleton sites, 0 constant sites

Akaike Information Criterion:           GTR+F+ASC+R2
Corrected Akaike Information Criterion: GTR+F+ASC+R2
Bayesian Information Criterion:         TVM+F+ASC+R2
Best-fit model: TVM+F+ASC+R2 chosen according to BIC

Total number of iterations: 120
CPU time used for tree search: 48769.380 sec (13h:32m:49s)
Wall-clock time used for tree search: 13213.018 sec (3h:40m:13s)
Total CPU time used: 85097.524 sec (23h:38m:17s)
Total wall-clock time used: 23168.309 sec (6h:26m:8s)

Computing bootstrap consensus tree...
Reading input file /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.splits.nex...
34 taxa and 98 splits.
Consensus tree written to /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.contree
Reading input trees file /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.contree
Log-likelihood of consensus tree: -31554816.563

Analysis results written to: 
  IQ-TREE report:                /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.iqtree
  Maximum-likelihood tree:       /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.treefile
  Likelihood distances:          /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.splits.nex
  Consensus tree:                /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.contree
  Screen log file:               /scratch/amonc/xipho_revision/iqtree/vcf2phylip/xipho.min4.phy.varsites.phy.log
```


  ### Using R script with phytools to collapse nodes with <50 ultra bootstrap support

  ### In FigTree v1.4.4 I uploaded the file to save as a pdf
   
  
  ### readjusted location of bootstrap support values in illustrator
