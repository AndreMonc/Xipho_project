
# Align new Xipho reference genome scaffolds to Chiroxiphia
# This will allow us to see what scaffolds align to sex chromosomes (Z and W), so that we can remove those from ARG analysis

## Download Chiroxiphia genome
## Number of Chromosomes 35 in the Chiroxiphia genome, however 93 scaffolds. I will only map to the 35 chromosomes
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/829/145/GCA_009829145.1_bChiLan1.pri/GCA_009829145.1_bChiLan1.pri_genomic.fna.gz
gzip -d GCA_009829145.1_bChiLan1.pri_genomic.fna.gz
```
GCA_009829145.1_bChiLan1.pri_genomic.fna # unzipped file name

## Rename Chiroxiphia scaffolds
#### Based on: https://tejashree1modak.github.io/bioblogs/fasta_rename/
```
cut -d ' ' -f1 your_file.fa > new_file.fa # trim everything from first space onwards
cut -d ' ' -f1 GCA_009829145.1_bChiLan1.pri_genomic.fna > GCA_009829145.1_bChiLan1.pri_genomic_trimmed.fna
python rename_fasta.py --mapping-file scaffold_rename.csv -i GCA_009829145.1_bChiLan1.pri_genomic_trimmed.fna -o /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/chiro_lanceo_ref.fa
```

## Rename Xiphorhynchus scaffolds to 1–477, in descending order of size
#### Based on: https://tejashree1modak.github.io/bioblogs/fasta_rename/
#### worked like a charm!
```
python rename_fasta.py --mapping-file scaffold_rename.csv -i xipele_purged.fasta.masked.mtDNAfiltered.fa -o /scratch/a_monc/postdoc/refs/Xiphorhynchus_elegans/xiph_elegans_ref.fa
```

## Create conda environment for ragtag
```
conda create -n "ragtag" 
conda activate ragtag
```

## Install minimap
```
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.28_x64-linux/minimap2

./minimap2-2.28_x64-linux/minimap2 # command to run the program
```

## Install unimap
```
git clone https://github.com/lh3/unimap
cd unimap && make

/scratch/a_monc/postdoc/unimap # command to run
```

## Install MUMmer
Downloaded from sourceforge, tar uploaded to cluster via filezilla
```
tar -xvzf MUMmer3.23.tar.gz
cd MUMmer3.23
make check #output says "check complete"
make install

/scratch/a_monc/postdoc/MUMmer3.23 # command to run
```
#### add paths of all these dependencies to bash profile

## Install ragtag v2.1.0
```
conda install -c bioconda ragtag
```

## HPC prompt example
```
#!/bin/bash
#PBS -A hpc_argweaver2
#PBS -l nodes=1:ppn=64
#PBS -l walltime=03:00:00
#PBS -q checkpt
#PBS -N ragtag

source activate ragtag

cd /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata

ragtag.py scaffold /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/chiro_lanceo_ref.fa /scratch/a_monc/postdoc/refs/Xiphorhynchus_elegans/xiph_elegans_ref.fa -e /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/exclude.txt -o /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/ragtag_output/
```

Inspect output .agp file:
```
ragtag.scaffold.agp
```

This shows that 33 scaffolds map to the sex chromosomes





