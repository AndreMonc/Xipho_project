4/21/2026
# PopCluster analysis for Xiphorhynchus revision
## Downloading popcluster (first had to enter my information on webpage and agree to terms and conditions)
```
wget https://cms.zsl.org/sites/default/files/2025-07/PopClusterLnx_15_07_2025.zip # this command says "Forbidden" on terminal, so downloaded and then uploaded via FileZilla to OSCER HPC
```

## Make the popcluster file executable
```
chmod u+x ./Bin/PopClusterLnx
```

## Test run, All Good!!
```
/home/amonc/PopCluster/Bin/PopClusterLnx INP:./Example/ant377NoScale.PcPjt
```

## Make directory for my analysis
mkdir popcluster

/scratch/amonc/xipho_revision/popcluster

## upload parameter file to xipho folder
```
xipho.PcPjt
```

## Uncompress my VCF file (I think it needs to be uncompressed for conversion to genotype format with GUI)
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=popc_unzip_%J_stdout.txt
#SBATCH --error=popc_unzip_%J_stderr.txt
#SBATCH --time=06:00:00
#SBATCH --job-name=popc_unzip
#SBATCH --chdir=/scratch/amonc/xipho_revision/vcf_filtering
#
#################################################

gunzip -c PopCluster.vcf.gz > PopCluster.vcf
```

## Copy vcf to external drive to transfer to windows computer
## convert to the dataframe=0 (one ind per line) option using the windows GUI

## Add new genotype file to PopCluster xipho folder

## Running PopCluster
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=popcluster_run_%J_stdout.txt
#SBATCH --error=popcluster_run_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=popcluster_run
#SBATCH --chdir=/scratch/amonc/xipho_revision/popcluster
#
#################################################

/home/amonc/PopCluster/Bin/PopClusterLnx INP:/scratch/amonc/xipho_revision/popcluster/xipho.PcPjt
```

#### Took like an hour to run
## The .K file has the good stuff
```
K                                                              BestRun         LogL_Mean          LogL_Min          LogL_Max          DLK1          DLK2       FST/FIS
1  "/scratch/amonc/xipho_revision/popcluster/xipho_popcluster_K_1_R_1" -2.4825076104E+06 -2.4825076104E+06 -2.4825076104E+06             -             -  8.832830E-01
2  "/scratch/amonc/xipho_revision/popcluster/xipho_popcluster_K_2_R_1" -2.2433378296E+06 -2.2433378296E+06 -2.2433378296E+06  1.012178E-01  1.209846E+00  2.040531E-01
3  "/scratch/amonc/xipho_revision/popcluster/xipho_popcluster_K_3_R_2" -2.1881547991E+06 -2.1881548047E+06 -2.1881547978E+06  2.490494E-02  1.876816E+00  1.700943E-01
4  "/scratch/amonc/xipho_revision/popcluster/xipho_popcluster_K_4_R_1" -2.1868858661E+06 -2.1868858689E+06 -2.1868858658E+06  5.800778E-04  1.963999E+00  1.341037E-01
5  "/scratch/amonc/xipho_revision/popcluster/xipho_popcluster_K_5_R_7" -2.1981532591E+06 -2.2045882550E+06 -2.1944647071E+06 -3.459591E-03             -  4.518285E-02

      Method      Best_K
        DLK2           4
     FST/FIS           1
```

# Create figures from best runs for K3 and K4 (similar DLK2)





### Added coordinates (lat and long columns) to the Q file for K3, Run 1

### Now opening in QGIS



