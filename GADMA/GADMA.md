# Goal:
#### I need to run different demographic models to test whether there was a bottleneck in the Belem population
# 29 May 2026

# install GADMA - a huge pain in the butt
# Here is my solution

conda create -n gadma_env10 python=3.10 \
  dadi h5py=3.10 pandas=2.2.2 matplotlib pillow scikit-allel \
  ruamel.yaml=0.16.12 numpy=1.26 scipy=1.13 \
  -c conda-forge -c bioconda -y

conda activate gadma_env10

python -m pip install moments-popgen==1.3.1

python -m pip install gadma==2.0.3 --no-deps

# Verify
which gadma
gadma --version

python -c "import gadma, moments, dadi; print('all imports OK')"
python -c "import moments; print(moments.__version__)"
python -c "import numpy; print(numpy.__version__)"
conda list | egrep "gadma|moments|dadi|numpy|scipy|ruamel"


## Pop file
/scratch/amonc/xipho_revision/GADMA/xiph_pops.txt

# Length of reference genome
awk '
/^>/ {next}
{sum += length($0)}
END {print sum}
' /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/ref/xiph_elegans_ref.fa

Output: 1120096857
# Checking reference length one more way
samtools faidx /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/ref/xiph_elegans_ref.fa
awk '{sum += $2} END {print sum}' /ourdisk/hpc/moncriefflab/amonc/dont_archive/xipho_revision/ref/xiph_elegans_ref.fa.fai

# Ok, the reference checks worked well
Now, I am filtering the VCF to a raw version with 33 individuals (removing outgroup) and only variant sites
That VCF will give me my X variable for the equation below


## Calculation of effective sequence length
L = (X - Y) / X * Nseq

Nseq, Total length of sequence (length of the reference genome, xipho_elegans_ragtagRef_no_W.fa) = 1120096857
X, Total number of snps received from this data (my raw SNP count from snpArcher) = 40761548
Y, SNPs filtered out (X - SNPs in GADMA VCF file), (40761548-89028) = 40672520

## Final Sequence calculation
L = (40761548 - 40672520) / 40761548 * 1120096857
L = 89028 / 40761548 * 1120096857
L = 2446423

## 
Unzip my vcf file and update .yaml file
Apparently GADMA won't run a .vcf.gz file
##

## Full run of GADMA with .yaml
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=10
#SBATCH --output=GADMA_%J_stdout.txt
#SBATCH --error=GADMA_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=GADMA
#SBATCH --chdir=/scratch/amonc/xipho_revision/GADMA
#
#################################################

source ~/.bashrc
conda activate gadma_env10

gadma -p param_file_xipho.yaml -o /scratch/amonc/xipho_revision/GADMA/run1
```


# End of output .out file (GADMA_32012156_stdout.txt):
# Successfully finished
--Finish pipeline--


If you use GADMA in your research please cite:
[GADMA]
Noskova et al., 2020: https://doi.org/10.1093/gigascience/giaa005
Noskova et al., 2023: https://doi.org/10.1093/gigascience/giad059
[Engine moments]
Jouganous et al., 2017: https://doi.org/10.1534/genetics.117.200493

More information about citations: https://gadma.readthedocs.io/en/latest/citations.html


Thank you for using GADMA!

In case of any questions or problems, please contact: ekaterina.e.noskova@gmail.com

# Saving summary of 20 replicate runs from end of GADMA.log file (for supplementary table)
Run #16 is best under both log-likelihood and AIC criteria
# I want to plot the optimal model across all 20 replicate runs


python best_aic_model_moments_code.py
python best_aic_model_moments_code.py

# Due to Matplotlib incompatibility issue while trying to plot, I need to try patching as indicated below
# GADMA plotting fix (gadma_env10)
# The demographic model plotting scripts produced by GADMA were incompatible
# with the installed Matplotlib version. Two patches to moments/ModelPlot.py
# were required before plots could be generated successfully.
# My conda prefix: /home/amonc/miniforge3/envs/gadma_env10
conda activate gadma_env10

# Backup ModelPlot.py
cp $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py \
   $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py.bak

# Fix deprecated matplotlib grid() syntax
sed -i \
  -e 's/grid(b=True/grid(visible=True/g' \
  -e 's/grid(b=False/grid(visible=False/g' \
  $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py

# Prevent figsize from being passed to savefig()
sed -i '/plt.savefig(save_file, \*\*fig_kwargs)/i\        fig_kwargs.pop("figsize", None)' \
  $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py

# Verify patches (optional)
grep -n "grid(" \
  $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py

grep -n "fig_kwargs.pop" \
  $CONDA_PREFIX/lib/python3.10/site-packages/moments/ModelPlot.py

# Generate demographic model figure
cd /scratch/amonc/xipho_revision/GADMA/run1
python best_aic_model_moments_code.py

# Now this worked!
python best_aic_model_moments_code.py

# Output file: 
/scratch/amonc/xipho_revision/GADMA/run1/model_from_GADMA.png



