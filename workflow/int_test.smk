
#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###             Integration testing snakefile for ATAC-seq                   ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

##################################
###   Load Required Packages   ###
##################################

import numpy as np
import os
import pandas as pd
import re

###########################
###   Variable Naming   ###
###########################

# Names directly from configuration YAML
atac_repo =      config['atac_repo']
configfile: atac_repo + "/config/int_common.yaml"
atac_peak_cut =  config['atac_peak_cut']
autosome_bed =   config['autosome_bed']
blacklist_bed =      config['blacklist_bed']
data_dir =       config['data_dir']
genome_fasta =   config['genome_fasta']
threads =        config['threads']

# Names derived from configuration YAML base
atac_bam_dir = data_dir + "/analysis/atac/bams"
atac_fastq_dir = data_dir + "/analysis/atac/fastqs"
atac_atac_keep_bed = data_dir + "/inputs/mm10chrs.bed"
atac_macs2_dir = data_dir + "/analysis/atac/macs2"
atac_script_dir = atac_repo + "/scripts"
bowtie2_dir =   config["data_dir"] + "/ref/hg38"
bowtie2_index = config["data_dir"] + "/ref/hg38/hg38"
atac_atac_keep_bed = config["data_dir"] + "/ref/atac_keep.bed"
qc_dir = config['data_dir'] + "/analysis/qc"
log_dir =        config['data_dir'] + "/logs"
sample_sheet = data_dir + "/inputs/libraries.tsv"

#####################
###   Functions   ###
#####################

def make_library_dictionary(tsv_path, col):
    lib_in = pd.read_table(tsv_path)
    readable = []
    for x in lib_in[col]:
        readable.append(os.access(x, os.R_OK))
    lib_in['readable']=readable
    lib = lib_in[lib_in.readable == True]
    lib_dict = dict(zip(lib['library'], lib[col]))
    return lib_dict

lib_dict = make_library_dictionary(sample_sheet, "path")

ATAC_LIBS = list(lib_dict.keys())

ATAC_GROUPS = pd.read_table(sample_sheet)['group'].unique().tolist()

###   Rules   ###

rule all:
    input:
        expand(atac_fastq_dir + "/{library}_{proc}_{read}.fastq.gz",
               library = ATAC_LIBS,
               proc = ['raw', 'proc'],
               read = ['R1', 'R2']),
        bowtie2_index + ".1.bt2",
        expand(atac_bam_dir + "/{library}_raw.bam",
               library = ATAC_LIBS),
        # expand(qc_dir +
        #        "/{library}_{processing}_{read}_fastqc.html",
        #        library = ATAC_LIBS,
        #        processing = ["raw","proc"],
        #        read = ["R1", "R2"]),
        # expand(qc_dir + "/{library}_{processing}_{stat}.txt",
        #        library = ATAC_LIBS,
        #        processing = ["raw","dedup","filt"],
        #        stat = ["flagstat", "samstats"]),
        # qc_dir + "/insert_sizes.tsv",
        # qc_dir + "/insert_sizes.pdf",
        # qc_dir + "/atac_qc.rdata",
        # expand(atac_macs2_dir +
        #        "/{library}_multi_peaks.narrowPeak",
        #        library = ATAC_LIBS),
        # expand(atac_macs2_dir +
        #        "/{library}_peaks.broadPeak",
        #        library = ATAC_LIBS),
        # expand(atac_macs2_dir +
        #        "/{library}_single_peaks.narrowPeak",
        #        library = ATAC_LIBS),
        # expand(atac_macs2_dir + "/{group}_consensus.bed",
        #        group = ATAC_GROUPS),
        # expand(atac_macs2_dir + "/{library}_naive.bed",
        #        library = ATAC_LIBS),
        # data_dir + "/ref/txdb",

###   Benchmark   ###
#+begin_src snakemake

# onsuccess:
#     shell("""
#         bash {cfdna_wgs_scriptdir}/agg_bench.sh {benchdir} {qc_dir}/agg_bench.tsv
#         """)

rule symlink_inputs:
    input: lambda wildcards: lib_dict[wildcards.library],
    output:
        read1 = atac_fastq_dir + "/{library}_raw_R1.fastq.gz",
        read2 = atac_fastq_dir + "/{library}_raw_R2.fastq.gz",
    params:
        out_dir = atac_fastq_dir,
        script = atac_script_dir + "/symlink.sh",
    shell:
        """
        {params.script} \
        {input} \
        {output.read1} \
        {output.read2} \
        {params.out_dir}
        """

include: atac_repo + "/workflow/atac_read_processing.smk"
include: atac_repo + "/workflow/atac_qc.smk"
#include: atac_repo + "/workflow/atac_peaks.smk"
