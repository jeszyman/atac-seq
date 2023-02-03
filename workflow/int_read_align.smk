print("Integration testing snakefile for ATAC-seq read and alignment processing\n")

# Import common packages
import pandas as pd
import re
import numpy as np

sample_sheet = config["data_dir"] + "/inputs/libraries.tsv"
atac_macs2_dir = config["data_dir"] + "/macs2"
atac_bam_dir = config["data_dir"] + "/atac_bam_dir"
blacklist_bed = config["blacklist_bed"]
threads = config["threads"]["default"]
atac_scripts = config["scriptdir"]["atac"]
log_dir =             config["log_dir"]
qc_dir = config["qc_dir"]
genome_fasta = config["genome_fasta"]
atac_fastq_dir =     config["data_dir"] + "/atac-fastq"
atac_container =     config["container"]["atac"]
atac_bam_raw =       config["data_dir"] + "/bam_atac_raw"
atac_bam_dedup =     config["data_dir"] + "/bam_atac_dedup"
atac_bam_filt =      config["data_dir"] + "/bam_atac_filt"

atac_bam_tn5 =       config["data_dir"] + "/bam_atac_tn5"
atac_bam_open =      config["data_dir"] + "/bam_atac_open"

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

ATAC_LIBS = list(lib_dict.keys())

rule all:
    input:
        expand(atac_fastq_dir + "/{library}_{processing}_{read}.fastq.gz",
            library = ATAC_LIBS,
            processing = ["raw", "proc"],
            read = ["R1", "R2"]),
        expand(atac_fastq_dir + "/{library}_proc_{read}.fastq.gz",
            library = ATAC_LIBS,
            read = ["R1", "R2"]),
	bowtie2_dir,
        expand(atac_bam_raw + "/{library}.bam",
            library = ATAC_LIBS),
        expand(atac_bam_dedup + "/{library}_dedup.bam",
            library = ATAC_LIBS),
        expand(atac_bam_filt + "/{library}_filt.bam",
            library = ATAC_LIBS),
        qc_dir + "/deeptools_frag_lengths.png",
        qc_dir + "/deeptools_frag_lengths.txt",
        expand(qc_dir + "/{library}_{processing}_{read}_fastqc.html",
            library = ATAC_LIBS,
            processing = ["raw", "proc"],
            read = ["R1", "R2"]),
        expand(qc_dir + "/{library}_{processing}_flagstat.txt",
            library = ATAC_LIBS,
            processing = ["raw", "dedup", "filt"]),
        expand(qc_dir + "/{library}_{processing}_samstats.txt",
            library = ATAC_LIBS,
            processing = ["raw", "dedup", "filt"]),
        qc_dir + "/atac_qc.rdata",
        #expand(atac_macs2_dir + "/{library}_peaks.broadPeak",
        #    library = glob_wildcards("test/qcpass/{library}_filt.bam")),
        #expand(atac_macs2_dir + "/{library}_{peakagg}_peaks.narrowPeak",
        #    library = get_qc_pass2,
        #    peakagg = ["single","multi"]),
        #expand(atac_macs2_dir + "/{group}_consensus.bed", group = ["ctrl", "exp"]),
        #expand(atac_macs2_dir + "/{library}_naive.bed",
        #    library = ATAC_LIBS),
        #qc_dir + "/insert_sizes.tsv",
        #qc_dir + "/insert_sizes.pdf",
        #expand(atac_bam_tn5 + "/{library}_tn5.bam", library = ATAC_LIBS),
        #expand(atac_bam_open + "/{library}_open.bam", library = ATAC_LIBS),
        #config["data_dir"] + "/ref/txdb",
        #config["data_dir"] + "/qc/atac_qc.rdata",
        #expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_broad}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_broad = MACS_BROAD_EXT),
        #expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_narrow}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_narrow = MACS_NARROW_EXT),
        #expand(config["data_dir"] + "/csaw/background_counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/csaw/counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/open_chrom/{library_id}_open_chrom.txt", library_id = LIBRARY_IDS),
	#expand(config["data_dir"] + "/csaw/norm_counts_rse_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/csaw/dge_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/dca_granges_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/{bam_process}_dca.csv", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/{bam_process}_chipseek.rds", bam_process = BAM_PROCESS),

rule symlink_read_align_inputs:
    container:
        atac_container,
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        r1 = atac_fastq_dir + "/{library}_raw_R1.fastq.gz",
        r2 = atac_fastq_dir + "/{library}_raw_R2.fastq.gz",
    shell:
        """
        r2=$(echo {input} | sed "s/_R1/_R2/g")
        ln -sf --relative {input} {output.r1}
        ln -sf --relative $r2 {output.r2}
        """

include: config["atac_repo"] + "/workflow/atac.smk"
