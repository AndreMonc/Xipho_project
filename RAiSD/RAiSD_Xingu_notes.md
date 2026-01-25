### 5/20/25

# I want to have an independent method to identify sweeps across the genome
### RAiSD seems like a decent option

## Probably want a VCF without centromeric regions (so repeatmasked)
## I think I could use my ARG mask

## RAiSD v3.1 looks like the most up-to-date version
```
https://github.com/pephco/RAiSD
```
## However, following directions leads to installing v2.9, so going with v2.9. v2.9 seems like more official version
```
https://github.com/alachins/raisd
```
```
mkdir RAiSD
cd RAiSD
wget https://github.com/alachins/raisd/archive/master.zip
unzip master.zip
cd raisd-master
./install-RAiSD.sh
```
## Run the following command from within /scratch/a_monc/postdoc/RAiSD/raisd-master every time before running RAiSD
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/gsl/lib
```
## Test Run (couldn't download test data)
```
wget 139.91.162.50/raisd_data/d1.tar.gz
tar -xvzf d1.tar.gz
./RAiSD -n test_run -I d1/msselection1.out -L 100000
```
## help command
```
./RAiSD -h
```
##  VCF filtering for RAiSD
#### As recommended in instructions, removing repetitive regions which could inflate estimates of mu.

## List individuals in VCF
```
/scratch/a_monc/postdoc/bcftools-1.3/bcftools query -l /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz
```
## Created list of Xingu individuals
```
/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xin_indiv.txt
```
## Xin indivs
```
Xin_xsCOUFT0424_Mus
Xin_xsFTA012_Mus
Xin_xsGAPTO037_Mus
Xin_xsGAPTO271_Mus
Xin_xsGAPX047_Mus
Xin_xsMAYA066_Mus
Xin_xsMOP011_Mus
Xin_xsMRJ498_Mus
Xin_xsTP36025_Toe
Xin_xsTP36276_Toe
Xin_xsTP48649_Toe
Xin_xsTP81164_Toe
Xin_xsUHE455_Mus
```

## Filter command
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N RAiSD_filt_xin

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --gzvcf /ddnA/work/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/FINAL_pseudo_XIPHO_raw.vcf.gz \
    --exclude-bed /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_GapRepeatIntervals_final.bed \
    --keep /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/xin_indiv.txt \
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
    --max-missing 0.5 \
    --recode \
    --out RAiSD_vcf_xin
```
After filtering, kept 14242451 out of a possible 46675869 Sites
/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/RAiSD_vcf_xin.recode.vcf


## RAiSD command example
```
./RAiSD -n vcf_run -I input_file.vcf -s
```

## RAiSD command actual
## Took about 45 min to run
```
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=72:00:00
#PBS -q checkpt
#PBS -N RAiSD_run2_xin

cd /scratch/a_monc/postdoc/RAiSD/raisd-master

./RAiSD -n vcf_run2_xin -I /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/RAiSD_vcf_xin.recode.vcf -s
```

### Perfect, this gave me a bunch of output report files with the mu value
### Now I want to get those report files into bed format
### For Xingu subset


## To run for all RAiSD report files in the same folder:
```
for infile in RAiSD_Report.vcf_run2_xin*; do
    # Strip any file extension (e.g., .txt, .vcf, etc.) if needed
    outfile="${infile}.bed"
    python RAiSD_to_bed.py --input_file "$infile" --output_file "$outfile"
done
```
## got a bunch of zsh: unknown file attribute, but I think that's ok. Everything looks ok.

## ----------------------------------------
## Now a different bash script to concatenate all non-empty bed files

## Name of the final combined BED file
```
output_file="combined_RAiSD_output.bed"
```
##  Remove existing output file to avoid appending to an old version
```
rm -f "$output_file"
```
## Loop through each .bed file and check if it's non-empty before concatenating
```
for file in *.bed; do
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
sort -k1,1 -k2,2n xipho_10kbwindows_112136.bed > xipho_10kbwindows_112136.sorted.bed
sort -k1,1 -k2,2n combined_RAiSD_output.bed > combined_RAiSD_output.sorted.bed
```
## Now, just to have this option available (not used in the end), I want to merge all short intervals to 10-kb windows (averaging mu values therein)
```
bedtools map -a xipho_10kbwindows_112136.sorted.bed \
             -b combined_RAiSD_output.sorted.bed \
             -c 4 -o mean > RAiSD_windowed_output.bed
```
## warning, but I think it's ok
***** WARNING: File combined_RAiSD_output.sorted.bed has inconsistent naming convention for record:
scaffold_113	398409	398410	6.741e+00

## Check depth and missingness for each individual to put in Table S1

## Filter command
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=64
#PBS -l walltime=6:00:00
#PBS -q checkpt
#PBS -N RAiSD_depth_missing_xin

cd /scratch/a_monc/postdoc/xipho_project/vcftools_filtering

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/RAiSD_vcf_xin.recode.vcf \
    --depth

vcftools \
    --vcf /ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/RAiSD_vcf_xin.recode.vcf \
    --missing-indv
```

### Because RAiSD interprets regions with few snps as more likely to indicate selective sweeps, 
### I need to filter out the long masked repetitive regions absent from the VCF used to estimate RAiSD u statistics
### (Otherwise bizarre peaks in the middle of repetitive regions masked in all my other genomic stats)
```
bedtools subtract -a combined_RAiSD_output.sorted.xingu.bed -b xipho_GapRepeatIntervals_final.bed > combined_RAiSD_output.sorted.cleaned.xingu.bed
```