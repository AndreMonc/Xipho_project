## 26 Sept 2024
## Ok, creating folder for new Xiphorhynchus snpArcher run with pseudochromosome approach (mapping scaffolds to Chiroxiphia with RagTag)

working directory:
```
/scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom
```
## Transfer files to new WD
```
cp Xiphorhynchus_sample_sheet.csv /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/
cp config.yaml /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/config/
cp config.yaml /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/slurm/
```

## make necessary edits to config files

## Location of data
/scratch/a_monc/postdoc/xipho_project/raw_data

## Location of reference
/scratch/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa

## Sample sheet
/scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/Xiphorhynchus_sample_sheet.csv

## Tmp storage
/scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/tmp_files

## slurm folder
/scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/slurm


# Steps to run snparcher
### running tmux on mike1 head node (note that you will not be able to access tmux session if you log into mike2 head node)
```
cd /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom
tmux new -s xipho_pseudo #give whatever session name you want
mamba activate snparcher
cd /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom
snakemake -s /scratch/a_monc/postdoc/snpArcher/workflow/Snakefile -d /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom --workflow-profile /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/slurm
```
# Ctrl-b d # use to detach from tmux window/session while program within continues to run
#tmux a -t xipho_pseudo #attaches to tmux session with specified name

# Once jobs finish (check with squeue), then need to create new window in the session, unlock directory and then resubmit job to keep things going
# New window
```
tmux a -t xipho_pseudo #attaches to tmux session with specified name
Ctrl-b c #makes new window in tmux session
Ctrl-b n #switches between windows in tmux session
```
# First unlock
```
mamba activate snparcher
cd /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom
snakemake -s /scratch/a_monc/postdoc/snpArcher/workflow/Snakefile -d /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom --workflow-profile /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/slurm --unlock
```
# Then run
```
snakemake -s /scratch/a_monc/postdoc/snpArcher/workflow/Snakefile -d /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom --workflow-profile /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/slurm
Ctrl-b d #detaches from tmux window/session while program within continues to run
```
# Tmux end session only once snpArcher run is finally completed
```
tmux kill-session -t session_name #kills tmux session with specified name
```
