print("Integration testing snakefile for ATAC-seq read and alignment processing\n")

# Import common packages
import pandas as pd
import re
import numpy as np

atac_peak_cut=2
samplesheet = config["datadir"] + "/inputs/libraries.tsv"
macs2dir = config["datadir"] + "/macs2"
atac_bams = config["datadir"] + "/atac_bams"
atac_blacklist = config["atac_blacklist"]
default_threads = config["threads"]["default"]
atac_scripts = config["scriptdir"]["atac"]
logdir =             config["logdir"]
qcdir = config["qcdir"]
atac_fasta = config["atac_fasta"]
atac_fastq =     config["datadir"] + "/atac-fastq"
atac_container =     config["container"]["atac"]
atac_bowtie2_index = config["datadir"] + "/ref/ucsc_mm10_chr19/ucsc_mm10_chr19"
atac_bowtie2_dir =   config["datadir"] + "/ref/ucsc_mm10_chr19"
atac_bam_raw =       config["datadir"] + "/bam_atac_raw"
atac_bam_dedup =     config["datadir"] + "/bam_atac_dedup"
atac_bam_filt =      config["datadir"] + "/bam_atac_filt"
atac_keep_bed =      config["datadir"] + "/inputs/mm10chrs.bed"
atac_bam_tn5 =       config["datadir"] + "/bam_atac_tn5"
atac_bam_open =      config["datadir"] + "/bam_atac_open"

libraries = pd.read_table("test/inputs/libraries.tsv")
libraries["r1_path"]="test/inputs/" + libraries["basename"]

readable = []
for x in libraries.r1_path:
    readable.append(os.access(x, os.R_OK))
libraries['readable']=readable

libraries = libraries[libraries.readable == True]

library_indict = libraries["library"].tolist()
file_indict = libraries["r1_path"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

ATAC_LIBRARIES = list(lib_dict.keys())

rule all:
    input:
        expand(atac_fastq + "/{library}_{processing}_{read}.fastq.gz",
            library = ATAC_LIBRARIES,
            processing = ["raw", "proc"],
            read = ["R1", "R2"]),
        expand(atac_fastq + "/{library}_proc_{read}.fastq.gz",
            library = ATAC_LIBRARIES,
            read = ["R1", "R2"]),
	atac_bowtie2_dir,
        expand(atac_bam_raw + "/{library}.bam",
            library = ATAC_LIBRARIES),
        expand(atac_bam_dedup + "/{library}_dedup.bam",
            library = ATAC_LIBRARIES),
        expand(atac_bam_filt + "/{library}_filt.bam",
            library = ATAC_LIBRARIES),
        qcdir + "/deeptools_frag_lengths.png",
        qcdir + "/deeptools_frag_lengths.txt",
        expand(qcdir + "/{library}_{processing}_{read}_fastqc.html",
            library = ATAC_LIBRARIES,
            processing = ["raw", "proc"],
            read = ["R1", "R2"]),
        expand(qcdir + "/{library}_{processing}_flagstat.txt",
            library = ATAC_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(qcdir + "/{library}_{processing}_samstats.txt",
            library = ATAC_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        qcdir + "/atac_qc.rdata",
        expand(macs2dir + "/{library}_peaks.broadPeak",
            library = glob_wildcards("test/qcpass/{library}_filt.bam")),
        #expand(macs2dir + "/{library}_{peakagg}_peaks.narrowPeak",
        #    library = get_qc_pass2,
        #    peakagg = ["single","multi"]),
        #expand(macs2dir + "/{group}_consensus.bed", group = ["ctrl", "exp"]),
        #expand(macs2dir + "/{library}_naive.bed",
        #    library = ATAC_LIBRARIES),
        #qcdir + "/insert_sizes.tsv",
        #qcdir + "/insert_sizes.pdf",
        #expand(atac_bam_tn5 + "/{library}_tn5.bam", library = ATAC_LIBRARIES),
        #expand(atac_bam_open + "/{library}_open.bam", library = ATAC_LIBRARIES),
        #config["datadir"] + "/ref/txdb",
        #config["datadir"] + "/qc/atac_qc.rdata",
        #expand(config["datadir"] + "/macs2/{library_id}_{bam_process}_{macs_broad}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_broad = MACS_BROAD_EXT),
        #expand(config["datadir"] + "/macs2/{library_id}_{bam_process}_{macs_narrow}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_narrow = MACS_NARROW_EXT),
        #expand(config["datadir"] + "/csaw/background_counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/csaw/counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/open_chrom/{library_id}_open_chrom.txt", library_id = LIBRARY_IDS),
	#expand(config["datadir"] + "/csaw/norm_counts_rse_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/csaw/dge_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/dca/dca_granges_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/dca/{bam_process}_dca.csv", bam_process = BAM_PROCESS),
        #expand(config["datadir"] + "/dca/{bam_process}_chipseek.rds", bam_process = BAM_PROCESS),

rule symlink_read_align_inputs:
    container:
        atac_container,
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        r1 = atac_fastq + "/{library}_raw_R1.fastq.gz",
        r2 = atac_fastq + "/{library}_raw_R2.fastq.gz",
    shell:
        """
        r2=$(echo {input} | sed "s/_R1/_R2/g")
        ln -sf --relative {input} {output.r1}
        ln -sf --relative $r2 {output.r2}
        """

include: config["atac_repo"] + "/workflow/atac.smk"
