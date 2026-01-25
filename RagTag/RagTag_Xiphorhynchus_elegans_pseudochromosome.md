
# Align new Xipho reference genome scaffolds to Chiroxiphia

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

## Stats for final pseudochromosome genome assembly
#### Excluded all the W chromosome scaffolds (including the assembled W chromosome) from the Chiroxiphia reference
#### Retained all the Xipho scaffolds that did not map to Chiroxiphia
```
cd /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/ragtag_output/results_without_W/xipho_elegans_ragtagRef_no_W.fa
gfastats xipho_elegans_ragtagRef_no_W.fa
```

## Summaries
```
+++Assembly summary+++: 
# scaffolds: 272
Total scaffold length: 1120117357
Average scaffold length: 4118078.52
Scaffold N50: 66779567
Scaffold auN: 74116822.86
Scaffold L50: 6
Largest scaffold: 156659620
Smallest scaffold: 15061
# contigs: 477
Total contig length: 1120096857
Average contig length: 2348211.44
Contig N50: 13035591
Contig auN: 14134131.56
Contig L50: 29
Largest contig: 42926112
Smallest contig: 15061
# gaps in scaffolds: 205
Total gap length in scaffolds: 20500
Average gap length in scaffolds: 100.00
Gap N50 in scaffolds: 100
Gap auN in scaffolds: 100.00
Gap L50 in scaffolds: 103
Largest gap in scaffolds: 100
Smallest gap in scaffolds: 100
Base composition (A:C:G:T): 319173808:240553719:240711926:319657404
GC content %: 42.97
# soft-masked bases: 170758451
# segments: 477
Total segment length: 1120096857
Average segment length: 2348211.44
# gaps: 205
# paths: 272
```

#### Used this reference, xipho_elegans_ragtagRef_no_W.fa, for the snpArcher pipeline
```
cp /scratch/a_monc/postdoc/refs/Chiroxiphia_lanceolata/ragtag_output/results_without_W/xipho_elegans_ragtagRef_no_W.fa /scratch/a_monc/postdoc/refs/Xiphorhynchus_elegans/
```





