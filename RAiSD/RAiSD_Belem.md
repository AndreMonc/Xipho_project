### 5/6/26

# I want to have an independent method to identify sweeps across the genome
### RAiSD seems like a decent option

### Probably want a VCF without centromeric regions (so repeatmasked)
### I think I could use my ARG mask

### RAiSD v3.1 looks like the most up-to-date version
```
https://github.com/pephco/RAiSD
```
### However, following directions leads to installing v2.9, so going with v2.9. v2.9 seems like more official version
```
https://github.com/alachins/raisd


mkdir RAiSD
cd RAiSD
wget https://github.com/alachins/raisd/archive/master.zip
unzip master.zip
cd raisd-master
./install-RAiSD.sh
```

## Help command
./RAiSD -h

## VCF filtering for RAiSD
## As recommended in instructions, removing repetitive regions which could inflate estimates of mu.

## RAiSD command example
```
./RAiSD -n vcf_run -I input_file.vcf -s
```

## RAiSD command actual
#### Took about 45 min to run
```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=Belem_RAiSD_%J_stdout.txt
#SBATCH --error=Belem_RAiSD_%J_stderr.txt
#SBATCH --time=3:00:00
#SBATCH --job-name=Belem_RAiSD
#SBATCH --chdir=/scratch/amonc/xipho_revision/RAiSD/belem
#
#################################################
module load GSL/1.16-goolf-1.4.10

/home/amonc/RAiSD/raisd-master/RAiSD -n bel_run1 -I /scratch/amonc/xipho_revision/vcf_filtering/Belem_RAiSD.vcf -s
```

#### Perfect, this gave me a bunch of output report files with the mu value
#### Now I want to get those report files into bed format

## To run for all RAiSD report files in the same folder:
```
for infile in RAiSD_Report.bel_run1*; do
    # Strip any file extension (e.g., .txt, .vcf, etc.) if needed
    outfile="${infile}.bed"
    python RAiSD_to_bed.py --input_file "$infile" --output_file "$outfile"
done
```
## ----------------------------------------
## Now a different bash script to concatenate all non-empty bed files

## Name of the final combined BED file
```
output_file="combined_RAiSD_output_belem.bed"
```
##  Remove existing output file to avoid appending to an old version
```
rm -f "$output_file"
```
## Loop through each .bed file and check if it's non-empty before concatenating
```
for file in RAiSD_Report.bel_run1.*bed; do
  if [ -s "$file" ]; then
    cat "$file" >> "$output_file"
  else
    echo "Skipping empty file: $file"
  fi
done

echo "Concatenation complete. Output saved to: $output_file"
```
# ----------------------------------------
## Sort bed files 
```
sort -k1,1 -k2,2n windows.bed > windows.sorted.bed
sort -k1,1 -k2,2n combined_RAiSD_output_belem.bed > combined_RAiSD_output_belem.sorted.bed
```
## Now, just to have this option available (not used in the end), I want to merge all short intervals to 10-kb windows (averaging mu values therein)
```
bedtools map -a xipho_10kbwindows_112136.sorted.bed \
             -b combined_RAiSD_output.sorted.bed \
             -c 4 -o mean > RAiSD_windowed_output.bed
```

### Because RAiSD interprets regions with few snps as more likely to indicate selective sweeps, 
### I need to filter out the long masked repetitive regions absent from the VCF used to estimate RAiSD u statistics
### (Otherwise bizarre peaks in the middle of repetitive regions masked in all my other genomic stats)
```
bedtools subtract -a combined_RAiSD_output_belem.sorted.bed -b repeat_regions.bed > combined_RAiSD_output_belem.sorted.cleaned.bed
```