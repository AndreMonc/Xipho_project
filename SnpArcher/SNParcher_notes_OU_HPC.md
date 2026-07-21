## 7 March 2026
# For revision of Xiphorhynchus manuscript

Creating new assembly of WGS data including X. elegans outgroup downloaded from NCBI

Also, using the original PacBio reference assembly of X. elegans (rather than pseudochromosome approach with scaffolds mapped to Chiroxiphia)

## First, download the X. elegans outgroup (so I can run ABBA-BABA tests)
Searching for X. elegans on NCBI
Two samples, both from B10K project, Illumina reads; BUT SAME INDIVIDUAL
https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=12254008

## First sample:
```
SRR9853667
Accession: SRX6608273
SAMN12254008
BioSample: SAMN12254008; Sample name: OUT-0059; SRA: SRS5172553
isolate	OUT-0059
age	not collected
sex	male
tissue	muscle
biomaterial provider	Jason Weir
collection date	2012-05-08
geographic location	Brazil
latitude and longitude	11.08 S 54.33 W
sample type	tissue sample
specimen voucher	MPEG:75162

Library:
Name: JW-0006_S7_L007
Instrument: Illumina HiSeq X Ten
Strategy: WGS
Source: GENOMIC
Selection: RANDOM
Layout: PAIRED
Runs: 1 run, 131.5M spots, 40.8G bases, 16.5Gb
Run	# of Spots	# of Bases	Size	Published
SRR9853667	131,528,811	40.8G	16.5Gb	2020-03-31
```

## Second sample (same individual):
```
SRR9853670
Accession: SRX6608270
SAMN12254008
BioSample: SAMN12254008; Sample name: OUT-0059; SRA: SRS5172553

isolate	OUT-0059
age	not collected
sex	male
tissue	muscle
biomaterial provider	Jason Weir
collection date	2012-05-08
geographic location	Brazil
latitude and longitude	11.08 S 54.33 W
sample type	tissue sample
specimen voucher	MPEG:75162

Library:
Name: JW-0006_S7_L006
Instrument: Illumina HiSeq X Ten
Strategy: WGS
Source: GENOMIC
Selection: RANDOM
Layout: PAIRED
Runs: 1 run, 129.3M spots, 40.1G bases, 16.2Gb
Run	# of Spots	# of Bases	Size	Published
SRR9853670	129,322,941	40.1G	16.2Gb	2020-03-31
```

## Find Runs for each sample
```
SRR9853667
SRR9853670
```
## First, Download SRA toolkit
```
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-alma_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz
```
I added bin folder to my bash profile path

## Command to prefetch
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -q single
#PBS -N Download_xipho_reads

cd /scratch/a_monc/postdoc/xipho_revision/SRA_download_outgroups

prefetch SRR9853667
prefetch SRR9853670
```

## Command to split, first run
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -q single
#PBS -N split_xipho_SRR9853667

cd /scratch/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853667

fasterq-dump SRR9853667 --split-files --skip-technical
```

## Command to split, second run
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -q single
#PBS -N split_xipho_SRR9853670

cd /scratch/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853670

fasterq-dump SRR9853670 --split-files --skip-technical
```

## zip my fastq files
```
#!/bin/bash
#PBS -A hpc_argweaver4
#PBS -l nodes=1:ppn=8
#PBS -l walltime=6:00:00
#PBS -q single
#PBS -N zip_fastqs

cd /scratch/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853667
gzip *.fastq

cd /scratch/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853670
gzip *.fastq
```

## Read locations (all same individual), I also renamed files 2 and 3 as 1 and 2
```
/ddnA/work/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853667/SRR9853667_1.fastq.gz
/ddnA/work/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853667/SRR9853667_2.fastq.gz
/ddnA/work/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853670/SRR9853670_1.fastq.gz
/ddnA/work/a_monc/postdoc/xipho_revision/SRA_download_outgroups/SRR9853670/SRR9853670_2.fastq.gz
```

## For reference genome: rename Xiphorhynchus scaffolds to 1–477, in descending order of size
Based on: https://tejashree1modak.github.io/bioblogs/fasta_rename/
### worked like a charm!

```
python rename_fasta.py --mapping-file scaffold_rename.csv -i xipele_purged.fasta.masked.mtDNAfiltered.fa -o /ddnA/work/a_monc/postdoc/xipho_revision/reference/xiph_elegans_ref.fa
```
## Check ref genome stats
```
conda install bioconda::gfastats # gfastats-1.3.11 
gfastats xiph_elegans_ref.fa

+++Assembly summary+++: 
# scaffolds: 477
Total scaffold length: 1120096857
Average scaffold length: 2348211.44
Scaffold N50: 13035591
Scaffold auN: 14134131.56
Scaffold L50: 29
Largest scaffold: 42926112
Smallest scaffold: 15061
# contigs: 477
Total contig length: 1120096857
Average contig length: 2348211.44
Contig N50: 13035591
Contig auN: 14134131.56
Contig L50: 29
Largest contig: 42926112
Smallest contig: 15061
# gaps in scaffolds: 0
Total gap length in scaffolds: 0
Average gap length in scaffolds: 0.00
Gap N50 in scaffolds: 0
Gap auN in scaffolds: 0.00
Gap L50 in scaffolds: 0
Largest gap in scaffolds: 0
Smallest gap in scaffolds: 0
Base composition (A:C:G:T): 319341574:240594133:240671512:319489638
GC content %: 42.97
# soft-masked bases: 170758451
# segments: 477
Total segment length: 1120096857
Average segment length: 2348211.44
# gaps: 0
# paths: 477
```

## snpArcher working directory:
Switching to OSCER HPC at University of Oklahoma due to memory issue error on LSU cluster
```
/scratch/amonc/xipho_revision/snpArcher_wd
```

## Location of data
```
/scratch/amonc/sra_project/fastq # X. spixii data
/scratch/amonc/xipho_revision/outgroup_reads # X. elegans data
```

## Location of reference
```
/scratch/amonc/xipho_revision/ref/xiph_elegans_ref.fa
```
## Sample sheet
```
/scratch/amonc/xipho_revision/snpArcher_wd/Xipho_sample_sheet.csv
```
## Tmp storage
```
/scratch/amonc/xipho_revision/snpArcher_wd/tmp_files
```
## slurm folder
```
/scratch/amonc/xipho_revision/snpArcher_wd/slurm
```
## config folder
```
/scratch/amonc/xipho_revision/snpArcher_wd/config
```

# Download Miniforge3
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh #saved to my home folder: /home/amonc/miniforge3
```
```
mamba create -c conda-forge -c bioconda -n snparcher "snakemake>=8" "python==3.11.4"
mamba activate snparcher
conda config --set channel_priority strict
pip install snakemake-executor-plugin-slurm
mamba deactivate
```
## Downloand snpArcher (version 1.0)
wget https://github.com/harvardinformatics/snparcher/archive/refs/tags/v1.0.zip
unzip v1.0.zip
cd snparcher-1.0

## All-sites update
Added an alternate bam2vcf_gatk_intervals.smk file, with rule to create all-sites VCF
This alternate file version is included in the snpArcher directory along with this notes file
```
/scratch/amonc/snparcher-1.0/workflow/rules/bam2vcf_gatk_intervals.smk
```

# Steps to run snparcher
Running tmux on schooner1 head node
# Using Snakemake 9.16.3
```
cd /scratch/amonc/xipho_revision/snpArcher_wd
tmux new -s xipho_revi #give whatever session name you want
mamba activate snparcher
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm 
```
Ctrl-b d #use to detach from tmux window/session while program within continues to run
tmux a -t xipho_revi #attaches to tmux session with specified name

Once jobs finish (check with squeue), then need to create new window in the session, unlock directory and then resubmit job to keep things going

```
cd /scratch/amonc/xipho_revision/snpArcher_wd
tmux a -t xipho_revi #attaches to tmux session with specified name
Ctrl-b c #makes new window in tmux session
Ctrl-b n #switches between windows in tmux session
```
First unlock
```
mamba activate snparcher
cd /scratch/amonc/xipho_revision/snpArcher_wd
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete --unlock
```
Then run
```
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm
Ctrl-b d #detaches from tmux window/session while program within continues to run
```
Tmux end session only once snpArcher run is finally completed
```
tmux kill-session -t session_name #kills tmux session with specified name
```


Helpful commands
```
scontrol show job 30135561 | grep -E 'TimeLimit|RunTime|JobState'
```

Everything ran super well from friday pm to monday, but then scattered issues of the next couple of weeks while snpArcher was running. Almost all issues traced to memory limits.

```
Status query timing (cumulative): calls=3225, failures=0, min=0.109s, avg=0.201s, max=2.152s
[Mon Mar 16 14:14:32 2026]
Finished jobid: 1952 (Rule: bam2gvcf)
800 of 2115 steps (38%) done
[Mon Mar 16 14:15:12 2026]
Finished jobid: 2008 (Rule: bam2gvcf)
801 of 2115 steps (38%) done
[Mon Mar 16 14:45:51 2026]
Finished jobid: 2005 (Rule: bam2gvcf)
802 of 2115 steps (38%) done
Shutting down, this might take some time.
Cleaning up SLURM log files older than 10 day(s).
Exiting because a job execution failed. Look below for error messages
[Mon Mar 16 14:46:32 2026]
Error in rule bam2gvcf:
    message: None
    jobid: 1776
    input: results/xiph_elegans_ref.fa/bams/Tap_xsPIME217_Mus_final.bam, results/xiph_elegans_ref.fa/bams/Tap_xsPIME217_Mus_final.bam.bai, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.sa, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.pac, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.bwt, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.ann, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.amb, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.fai, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.dict, results/xiph_elegans_ref.fa/intervals/gvcf_intervals/0022-scattered.interval_list
    output: results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsPIME217_Mus/0022.raw.g.vcf.gz, results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsPIME217_Mus/0022.raw.g.vcf.gz.tbi
    log: logs/xiph_elegans_ref.fa/gatk_hc/Tap_xsPIME217_Mus/0022.txt (check log file(s) for error details)
    conda-env: /scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/conda/8ebe92f6ac1c4c917cd5b24f85e7cb72_
    shell:
        
        gatk HaplotypeCaller         --java-options -Xmx56000m         -R results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna         -I results/xiph_elegans_ref.fa/bams/Tap_xsPIME217_Mus_final.bam         -O results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsPIME217_Mus/0022.raw.g.vcf.gz         -L results/xiph_elegans_ref.fa/intervals/gvcf_intervals/0022-scattered.interval_list         -ploidy 2         --emit-ref-confidence GVCF --min-pruning 2 --min-dangling-branch-length 4 &> logs/xiph_elegans_ref.fa/gatk_hc/Tap_xsPIME217_Mus/0022.txt
        
        (command exited with non-zero exit code)
Removing -L we added from /home/amonc/.curlrc...
Complete log(s): /scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/log/2026-03-13T142618.295082.snakemake.log
WorkflowError:
At least one job did not complete successfully.
```
# Will run snpArcher again with the --rerun-incomplete
```
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```
# March 20, 2026, all jobs completed
```
Status query timing (cumulative): calls=3000, failures=0, min=0.111s, avg=0.175s, max=1.651s
[Fri Mar 20 01:59:43 2026]
Finished jobid: 51 (Rule: bcftools_norm)
465 of 1313 steps (35%) done
[Fri Mar 20 02:06:19 2026]
Finished jobid: 39 (Rule: bcftools_norm)
466 of 1313 steps (35%) done
Status query timing (cumulative): calls=3025, failures=0, min=0.111s, avg=0.175s, max=1.651s
[Fri Mar 20 02:49:02 2026]
Finished jobid: 1205 (Rule: bam2gvcf)
467 of 1313 steps (36%) done
[Fri Mar 20 02:52:45 2026]
Finished jobid: 1134 (Rule: bam2gvcf)
468 of 1313 steps (36%) done
Shutting down, this might take some time.
Cleaning up SLURM log files older than 10 day(s).
Exiting because a job execution failed. Look below for error messages
[Fri Mar 20 02:53:26 2026]
Error in rule dedup:
    message: None
    jobid: 211
    input: results/xiph_elegans_ref.fa/bams/postMerge/XELEGANS_MPEG75162_Mus.bam, results/xiph_elegans_ref.fa/bams/postMerge/XELEGANS_MPEG75162_Mus.bam.bai
    output: results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam, results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam.bai
    log: logs/xiph_elegans_ref.fa/sambamba_dedup/XELEGANS_MPEG75162_Mus.txt (check log file(s) for error details)
    conda-env: /scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/conda/343f93b935638e95c142069f676165b7_
    shell:
        sambamba markdup -t 16 results/xiph_elegans_ref.fa/bams/postMerge/XELEGANS_MPEG75162_Mus.bam results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam 2> logs/xiph_elegans_ref.fa/sambamba_dedup/XELEGANS_MPEG75162_Mus.txt
        (command exited with non-zero exit code)
Removing -L we added from /home/amonc/.curlrc...
Complete log(s): /scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/log/2026-03-16T181054.415384.snakemake.log
WorkflowError:
At least one job did not complete successfully.
```
# Then I looked inside the log file:
```
sambamba 0.8.0
 by Artem Tarasov and Pjotr Prins (C) 2012-2020
    LDC 1.20.0 / DMD v2.090.1 / LLVM7.0.0 / bootstrap LDC - the LLVM D compiler (0.17.6)

finding positions of the duplicate reads in the file...
sambamba-markdup: Cannot open file `/tmp/sambamba-pid16156-markdup-emmy/PairedEndsInfoogra15' in mode `w+' (Too many open files)
XELEGANS_MPEG75162_Mus.txt (END)
```

I think I can just restart this job and see what happens.
First, raise the number of CPUs per task for rule dedup?
Reduce number of threads
If that doesn't work, I can update the dedup rule in slum config file to a bigger memory partition: largemem
sambamba 
```
"sambamba markdup --overflow-list-size 600000 -t {threads} {input.bam} {output.dedupBam} 2> {log}" 
```
# Running snpArcher again with the --rerun-incomplete
```
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```
# Helpful commands
```
scontrol show job 30214515 | grep -E 'TimeLimit|RunTime|JobState'
scontrol show job 30214516 | grep -E 'TimeLimit|RunTime|JobState'
scontrol show job 30214693 | grep -E 'TimeLimit|RunTime|JobState'
```

```
  dedup:
    mem_mb: 64000 # old: attempt * 64000
    mem_mb_reduced: 60000 # old: (attempt * 64000) * 0.9 # Mem allocated to java (tries to prevent OOM errors)
    slurm_partition: "normal"
    slurm_account: "general"
    cpus_per_task: 2 # old: 16
    runtime: "48h"
```
```
"sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"
```
I am having trouble running sambamba markdup for one particularly large BAM file.

Here is my specific issue in the following log file (logs/xiph_elegans_ref.fa/sambamba_dedup/XELEGANS_MPEG75162_Mus.txt file):
```
"sambamba 0.8.0
 by Artem Tarasov and Pjotr Prins (C) 2012-2020
    LDC 1.20.0 / DMD v2.090.1 / LLVM7.0.0 / bootstrap LDC - the LLVM D compiler (0.17.6)

finding positions of the duplicate reads in the file...
sambamba-markdup: Cannot open file `/tmp/sambamba-pid64344-markdup-nejx/PairedEndsInfocbin15' in mode `w+' (Too many open files)"

These are my current slurm settings in snpArcher for the dedup rule (which runs "sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"):   
dedup:
    mem_mb: 64000 # old: attempt * 64000
    mem_mb_reduced: 60000 # old: (attempt * 64000) * 0.9 # Mem allocated to java (tries to prevent OOM errors)
    slurm_partition: "normal"
    slurm_account: "general"
    cpus_per_task: 2 # old: 16
    runtime: "48h"
```
Here is the log of what is actually running:
```
rule dedup:
    input: results/xiph_elegans_ref.fa/bams/postMerge/XELEGANS_MPEG75162_Mus.bam, results/xiph_elegans_ref.fa/bams/postMerge/XELEGANS_MPEG75162_Mus.bam.bai
    output: results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam, results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam.bai
    log: logs/xiph_elegans_ref.fa/sambamba_dedup/XELEGANS_MPEG75162_Mus.txt
    jobid: 211
    benchmark: benchmarks/xiph_elegans_ref.fa/sambamba_dedup/XELEGANS_MPEG75162_Mus.txt
    reason: Missing output files: results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam.bai, results/xiph_elegans_ref.fa/bams/XELEGANS_MPEG75162_Mus_final.bam
    wildcards: refGenome=xiph_elegans_ref.fa, sample=XELEGANS_MPEG75162_Mus
    threads: 16
    resources: tmpdir=<TBD>, disk_mb=79986, disk=79.99 GB, disk_mib=76281, mem_mb=64000, mem=64 GB, mem_mib=61036, mem_mb_reduced=60000, slurm_partition=normal, slurm_account=general, cpus_per_task=2, runtime=2880
```
Here is some interesting text online I found about optimizing sambamba (although not specifically sambamba markdup):
"Ok, for what I expect will be my final update for this answer, I compared the following:

samtools 1.14 + zlib
samtools 1.15 + libdeflate (this is a different version, but shouldn't have a major effect. They released a new version just before I asked cluster admins to compile against libdeflate and didn't notice until later.)
sambamba v0.8.2
Note: This time I used the following settings:

Only two samples instead of three
Max Nextflow queue size of 30 to avoid too many threads reading from the same two files.
I tested CPU requests from 1 to 7 (step size 1) and then from 9 to 17 (step size 2).
For jobs where I allotted few CPUs and low memory, I provided a 200% buffer to prevent the jobs from failing with OUT OF MEMORY errors. i.e., I gave the SLURM job 2x the amount of memory than I told the tools they could use. sambamba was the worst offender for exceeding memory in these conditions, but I gave them all 200% buffer to keep comparisons equal.
Seems pretty clear to me that sambamba is the fastest and uses requested CPU resources most efficiently (i.e., is able to use all of the CPUs requested up until ~8 or 9 CPUs). For reference, the difference between sambamba and samtools + libdeflate with 9 CPUs is ~10 minutes on average.

We will certainly be using sambamba going forward, probably with ~9 CPUs and ~9GB per thread since I had to provide a 200% buffer in lower CPU/Mem jobs to avoid job failure, anyway.

Again, all code and results are available on my GitHub repository (Samtools sort optimization test)."

Also, I found some text suggesting that "--overflow-list-size 600000" might be good to add to my sambamba markdup command.

# My update to dedup rule
```
rule dedup:
    input:
	unpack(dedup_input)
    output:
	dedupBam = "results/{refGenome}/bams/{sample}_final.bam",
        dedupBai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    conda:
	"../envs/sambamba.yml"
    log:
        "logs/{refGenome}/sambamba_dedup/{sample}.txt"
    benchmark:
	"benchmarks/{refGenome}/sambamba_dedup/{sample}.txt"
    params:
         tmpdir = lambda wildcards: f"/scratch/amonc/sambamba_markdup/{wildcards.sample}",
         overflow = 1000000,
         hashsize = 1000000
    shell:
         r"""
         set -euo pipefail
         mkdir -p {params.tmpdir}
         sambamba markdup \
             -t {threads} \
             --tmpdir {params.tmpdir} \
             --overflow-list-size {params.overflow} \
             --hash-table-size {params.hashsize} \
             {input.bam} {output.dedupBam} \
             2> {log}
         """
```

# Running snpArcher again with the --rerun-incomplete
```
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```

# 24 March 2026--Had a bunch of errors in bam2gvcf haplotype caller, appears to be a memory issue
```
Error in rule bam2gvcf:
    message: None
    jobid: 817
    input: results/xiph_elegans_ref.fa/bams/Tap_xsMPDS1217_Mus_final.bam, results/xiph_elegans_ref.fa/bams/Tap_xsMPDS1217_Mus_final.bam.bai, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.sa, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.pac, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.bwt, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.ann, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.amb, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna.fai, results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.dict, results/xiph_elegans_ref.fa/intervals/gvcf_intervals/0018-scattered.interval_list
    output: results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsMPDS1217_Mus/0018.raw.g.vcf.gz, results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsMPDS1217_Mus/0018.raw.g.vcf.gz.tbi
    log: logs/xiph_elegans_ref.fa/gatk_hc/Tap_xsMPDS1217_Mus/0018.txt (check log file(s) for error details)
    conda-env: /scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/conda/8ebe92f6ac1c4c917cd5b24f85e7cb72_
    shell:
        
        gatk HaplotypeCaller         --java-options -Xmx56000m         -R results/xiph_elegans_ref.fa/data/genome/xiph_elegans_ref.fa.fna         -I results/xiph_elegans_ref.fa/bams/Tap_xsMPDS1217_Mus_final.bam         -O results/xiph_elegans_ref.fa/interval_gvcfs/Tap_xsMPDS1217_Mus/0018.raw.g.vcf.gz         -L results/xiph_elegans_ref.fa/intervals/gvcf_intervals/0018-scattered.interval_list         -ploidy 2         --emit-ref-confidence GVCF --min-pruning 2 --min-dangling-branch-length 4 &> logs/xiph_elegans_ref.fa/gatk_hc/Tap_xsMPDS1217_Mus/0018.txt
```

# gatk Haplotype caller command in rule bam2gvcf for reference:
  ```
        gatk HaplotypeCaller \
        --java-options -Xmx{resources.mem_mb_reduced}m \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -L {input.l} \
        -ploidy {params.ploidy} \
        --emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang} &> {log}

```
# Helpful commands
```
scontrol show job 30263252
```

# Appears to be memory issue
```
Xin_xsTP36276_Toe/0002.raw.g.vcf.gz
logs/xiph_elegans_ref.fa/gatk_hc/Tap_xsMSF111_Mus/0012.txt
logs/xiph_elegans_ref.fa/gatk_hc/Xin_xsFTA012_Mus/0001.txq
```

# Updated slurm for various rules including bam2gvcf to bump up the memory
# Running snpArcher again without run complete (I saw in the log that potentially corrupted files were removed)
```
tmux a -t xipho_revi #attaches to tmux session with specified name
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm
```

# 26 March 2026
I canceled all jobs (only four bam2gvcf running). I bumped down resources required--especially the amount of time for this rule (48 hrs to 10 hrs). Hopefully that will allow more jobs to run simultaneously. Here are the intervals that did not finish, I want to see if snpArcher automatically wipes them:
```
xiph_elegans_ref.fa_Tap_xsPIME217_Mus_0004
xiph_elegans_ref.fa_Bel_xsFRC041_Mus_0004
xiph_elegans_ref.fa_Bel_xsGUR156_Mus_0011
xiph_elegans_ref.fa_Xin_xsMRJ498_Mus_0011
```
```
tmux a -t xipho_revi #attaches to tmux session with specified name
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --unlock
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```

I need to check that the four files above are either removed automatically or that I remove them ***
They should get flagged and removed 
```
xiph_elegans_ref.fa_Tap_xsA08267_Mus_0001
```

Got email from OU hpc staff "Horst", that few public nodes have more than 64 GB Ram
I am dropping bam2gvcf to one thread and 64 GB ram, 4 CPUs per task

# Rerunning
```
mamba activate snparcher
tmux a -t xipho_revi #attaches to tmux session with specified name
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --unlock
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```

# I cancelled a bunch of jobs to see if I can run bam2gvcf faster (night of March 26)
# Rerunning
```
tmux a -t xipho_revi #attaches to tmux session with specified name
mamba activate snparcher
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --unlock
# gonna try without rerun incomplete
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm
```
# I am only running 5 jobs at a time. Not fast enough. Dropping resources for bam2gvcf to see if that helps.

# Ok, I think bam2gvcf finished!!!! 29 March 2026
# Actually a couple jobs appear to have failed
```
/scratch/amonc/xipho_revision/snpArcher_wd/.snakemake/log/2026-03-27T181714.100516.snakemake.log
interval_gvcfs/Tap_xsA08267_Mus/0024.raw.g.vcf.gz
interval_gvcfs/Xin_xsTP36276_Toe/0027.raw.g.vcf.gz
```
# Bumping up bam2gvcf resources and time
```
snakemake -s /scratch/amonc/snparcher-1.0/workflow/Snakefile -d /scratch/amonc/xipho_revision/snpArcher_wd --workflow-profile /scratch/amonc/xipho_revision/snpArcher_wd/slurm --rerun-incomplete
```

# snpArcher run completed April 1. 
Note to future self--don't update snakemake rules unless absolutely necessary (snpArcher run will get set back). Better to update # the resource rules only, just to get any stubborn samples through the bottleneck. That is almost always sufficient.