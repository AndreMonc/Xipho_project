# Dsuite calculation of D statistics
# 4/21/26

# VCF_to_use
/scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz

# Make population file
```
bcftools query -l /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz | sort > vcf.samples.txt
```
```
awk 'BEGIN{OFS="\t"}
{
    if ($1 == "XELEGANS_MPEG75162_Mus")
        print $1, "Outgroup"
    else {
        split($1, a, "_")
        print $1, a[1]
    }
}' vcf.samples.txt > SETS.txt
```
/scratch/amonc/xipho_revision/dsuite/SETS.txt

# Download Dsuite Version: 0.5 r58
module load GCCcore/11.3.0
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make

# Run Dsuite Dtrios
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=Dtrios_%J_stdout.txt
#SBATCH --error=Dtrios_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=Dtrios
#SBATCH --chdir=/scratch/amonc/xipho_revision/dsuite
#
#################################################

module load GCCcore/11.3.0

Dsuite Dtrios /scratch/amonc/xipho_revision/vcf_filtering/Dstats.vcf.gz /scratch/amonc/xipho_revision/dsuite/SETS.txt
```


