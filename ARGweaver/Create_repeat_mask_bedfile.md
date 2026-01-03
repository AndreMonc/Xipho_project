# Background:
# Greg Thom, in Fall 2024, used RepeatMasker to soft-mask the new Xiphorhynchus elegans reference genome. He did this before I created the pseudochromosome reference (aligned to Chiroxiphia using RagTag)
# Thus, the bed file output from RepeatMasker will not have the right coordinates for the repeat regions in the pseudochromosome reference

# Goal: identify soft-masked regions and create bedfile of these regions for the new pseudochromosome reference: xipho_elegans_ragtagRef_no_W.fa

# Goal is to get bed file for all soft masked regions
1) Convert soft-masked to hardmasked
https://github.com/fulcrumgenomics/fgbio HardMaskFasta function
https://fulcrumgenomics.github.io/fgbio/tools/latest/HardMaskFasta.html

conda create -n soft_mask_intervals #on smic
conda activate soft_mask_intervals
conda install fgbio

fgbio HardMaskFasta --input=/ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa --output=/ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W_hardMasked.fa

2) Create bed file of hard-masked regions
# This bed file will also include any gap regions. Hitting two birds with one stone.
https://www.biostars.org/p/382034/ #helpful material here

wget https://github.com/hewm2008/Reseqtools/archive/refs/tags/v0.25.tar.gz
tar -xvzf v0.25.tar.gz
tar -xvzf iTools_Code20180520.tar.gz

# Before converting soft masked to hard-masked (curious how many Ns from Gaps there are)
/ddnA/work/a_monc/postdoc/Reseqtools-0.25/iTools_Code/iTools Fatools findN -InPut /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa -OutPut /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W_hardMaskedIntervalsGapsOnly.bed

# Now, getting intervals for hard-masked regions (due to both gaps and soft-masking)
/ddnA/work/a_monc/postdoc/Reseqtools-0.25/iTools_Code/iTools Fatools findN -InPut /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W_hardMasked.fa -OutPut /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W_hardMaskedIntervals_ALL.bed

3) Merge hard-masked intervals (due to original gaps + soft masking from repeatmaster) that are separated by 100 bp or less
bedtools merge -d 100 -i /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W_hardMaskedIntervals_ALL.bed > xipho_GapRepeatIntervals_final.bed

# Final gap + repeat mask to use with ARGweaver and for VCF filtering
/ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_GapRepeatIntervals_final.bed


