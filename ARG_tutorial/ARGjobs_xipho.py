# !/usr/bin/env python
# encoding: utf-8

"""
Create the ARG qsub files
Copyright 2024 Andre E. Moncrieff. All rights reserved.
"""

import argparse
import pandas
import numpy
import csv
import itertools
import re


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ARG_info", required=True,
                        help="Enter the file name (including .txt extension) for job name list file",
                        type=str)
    args = parser.parse_args()
    return args


def read_in_csv(txt_file_dataframe):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str)
    return raw_dataframe


def build_job(jobs_file):
    job_name = jobs_file['job_name'].tolist()
    vcf_file = jobs_file['vcf_file'].tolist()
    region = jobs_file['region'].tolist()
    output_base = jobs_file['output_base'].tolist()
    
    for num in  range(len(job_name)):
        part1 = ['#!/bin/bash',
                    '#PBS -A hpc_argweaver3',
                    '#PBS -l nodes=1:ppn=4',
                    '#PBS -l walltime=168:00:00',
                    '#PBS -q single']
        part2 = ['#PBS -N ' + str(job_name[num])]
        part3 = ['\n', 'cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/output', '\n']
        part4 = ['/scratch/a_monc/postdoc/ARGweaver/bin/arg-sample --vcf /scratch/a_monc/postdoc/xipho_project/ARGweaver/vcf2MBwindows/' + str(vcf_file[num] + ' \\')]
        part5 = ['--region ' + str(region[num]) + ' \\']
        part6 = ['--vcf-min-qual 30 \\',
                '--vcf-genotype-filter "DP<5;DP>50;GQ<20" \\',
                '--ind-maskmap /ddnA/work/a_monc/postdoc/xipho_project/Depth_from_BAMs_500bp_windows/ind_mask_file.txt \\',
                '--maskmap /ddnA/work/a_monc/postdoc/xipho_project/ReLERNN/xipho_Argweaver_truly_final_mask.bed \\',
                '--mask-cluster 2,5 \\',
                '--mutrate 6.9e-9 \\',
                '--recombmap /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARGweaver_recomb_map_ReLERNN.bed \\',
                '--popsize 391476 \\',
                '--compress-seq 5 \\',
                '--ntimes 20 \\',
                '--maxtime 1e7 \\',
                '--delta 0.005 \\',
                '--sample-step 50 \\',
                '--iters 2000 \\']
        part7 = ['-o ' + str(output_base[num])]
        part1.extend(part2)
        part1.extend(part3)
        part1.extend(part4)
        part1.extend(part5)
        part1.extend(part6) 
        part1.extend(part7)

        with open('/Users/andremoncrieff/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Running_ARGweaver/qsub_file_prep/qsubs_without_resume_argweaver3/' + str(job_name[num]) + '.qsub', 'w') as f:
            for item in part1:
                f.write("%s\n" % item)
        
    
def main():
    #create args object
    args = parser()
    #read in dataframes
    jobs_file = read_in_csv(args.ARG_info)
    build_job(jobs_file)

if __name__ == '__main__':
    main()