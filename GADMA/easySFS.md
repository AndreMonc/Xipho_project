## April 21, 2026

# Running easySFS to generate projections appropriate for my dataset
```
cd /scratch/amonc/xipho_revision/GADMA
```
## Install easySFS
## https://github.com/isaacovercast/easySFS
```
conda create -n easySFS
conda activate easySFS
conda install -c conda-forge numpy pandas scipy -y
git clone https://github.com/isaacovercast/easySFS.git
cd easySFS
chmod 777 easySFS.py
./easySFS.py
```


# Run easySFS
# 28 May 2026
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=easySFS_%J_stdout.txt
#SBATCH --error=easySFS_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=easySFS
#SBATCH --chdir=/scratch/amonc/xipho_revision/GADMA
#
#################################################

source ~/.bashrc
conda activate easySFS

## Create variables
VCF="/scratch/amonc/xipho_revision/GADMA/GADMA.vcf.gz"
pop_file="/scratch/amonc/xipho_revision/GADMA/xiph_pops.txt"

## estimate projections
/home/amonc/easySFS/easySFS.py -i $VCF -p $pop_file -a -f --preview
```

# Waiting here on 28 May 2026




# Results
# Text from the .out file after running above easySFS script
```
Running preview mode. We will print out the results for # of segregating sites
    for multiple values of projecting down for each population. The dadi
    manual recommends maximizing the # of seg sites for projections, but also
    a balance must be struck between # of seg sites and sample size.
Bel
(2, 11689)      (3, 17534)      (4, 21294)      (5, 24011)      (6, 26117)      (7, 27825)      (8, 29256)      (9, 30485)      (10, 31562)     (11, 32518)     (12, 33379)     (13, 34161)     (14, 34877)     (15, 35539)     (16, 36154)     (17, 36728)     (18, 37267)     (19, 30263)     (20, 30644)     

Tap
(2, 15089)      (3, 22634)      (4, 27842)      (5, 31882)      (6, 35222)      (7, 38098)      (8, 40641)      (9, 42935)      (10, 45035)     (11, 46977)     (12, 48789)     (13, 50493)     (14, 52103)     (15, 53632)     (16, 55089)     (17, 56484)     (18, 57822)     (19, 46042)     (20, 47000)     

Xin
(2, 13135)      (3, 19703)      (4, 24078)      (5, 27356)      (6, 29986)      (7, 32191)      (8, 34098)      (9, 35784)      (10, 37301)     (11, 38682)     (12, 39955)     (13, 41136)     (14, 42240)     (15, 43278)     (16, 44259)     (17, 45190)     (18, 46076)     (19, 46924)     (20, 47736)     (21, 48516)     (22, 49267)     (23, 49991)     (24, 50691)     (25, 30456)     (26, 30837)     
```

For Tap, Xin, Bel, I selected the following projection: 18, 24, 18 (seems pretty clearly the best--maximizes # seg sites with only dropping a few individuals)