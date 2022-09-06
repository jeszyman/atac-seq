atac_fastq_raw =     config["datadir"] + "/fastq_atac_raw"
atac_fastq_proc =    config["datadir"] + "/fastq_atac_proc"
logdir =             config["logdir"]
atac_container =     config["container"]["atac"]
atac_bowtie2_index = config["datadir"] + "/ref/ucsc_mm10_chr19/ucsc_mm10_chr19"
atac_bowtie2_dir =   config["datadir"] + "/ref/ucsc_mm10_chr19"
atac_bam_raw =       config["datadir"] + "/bam_atac_raw"
atac_bam_dedup =     config["datadir"] + "/bam_atac_dedup"
atac_bam_filt =      config["datadir"] + "/bam_atac_filt"
atac_keep_bed =      config["datadir"] + "/inputs/mm10chrs.bed"
atac_bam_tn5 =       config["datadir"] + "/bam_atac_tn5"
atac_bam_open =      config["datadir"] + "/bam_atac_open"

import pandas as pd
import re
import numpy as np

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
        expand(atac_fastq_raw + "/{library}_R1.fastq.gz", library = ATAC_LIBRARIES),
        expand(atac_fastq_raw + "/{library}_R2.fastq.gz", library = ATAC_LIBRARIES),
        expand(atac_fastq_proc + "/{library}_flex_1.fastq.gz", library = ATAC_LIBRARIES),
        expand(atac_fastq_proc + "/{library}_flex_2.fastq.gz", library = ATAC_LIBRARIES),
	atac_bowtie2_dir,
        expand(atac_bam_raw + "/{library}.bam",	library = ATAC_LIBRARIES),
        expand(atac_bam_dedup + "/{library}_dedup.bam", library = ATAC_LIBRARIES),
        expand(atac_bam_filt + "/{library}_filt.bam", library = ATAC_LIBRARIES),
        expand(atac_bam_tn5 + "/{library}_tn5.bam", library = ATAC_LIBRARIES),
        expand(atac_bam_open + "/{library}_open.bam", library = ATAC_LIBRARIES),
        config["datadir"] + "/ref/txdb",
        config["datadir"] + "/qc/atac_qc.rdata",
        expand(config["datadir"] + "/qc/{library}_R{read}_fastqc.html", library = ATAC_LIBRARIES, read=["1","2"]),
        expand(config["datadir"] + "/qc/{library}_flex_{read}_fastqc.html", library = ATAC_LIBRARIES, read=["1","2"]),
        expand(config["datadir"] + "/qc/{library}_filt_stat.txt", library = ATAC_LIBRARIES),
        expand(config["datadir"] + "/qc/{library}_filt_flagstat.txt", library = ATAC_LIBRARIES),

        #expand(config["qc_dir"] + "/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1", "R2"]),
        #expand(config["qc_dir"] + "/{library_id}_stat.txt", library_id = LIBRARY_IDS),
        #expand(config["qc_dir"] + "/{library_id}_flagstat.txt", library_id = LIBRARY_IDS),
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

rule symlink:
    container:
        atac_container,
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        r1 = atac_fastq_raw + "/{library}_R1.fastq.gz",
        r2 = atac_fastq_raw + "/{library}_R2.fastq.gz",
    shell:
        """
        r2=$(echo {input} | sed "s/_R1/_R2/g")
        ln -sf --relative {input} {output.r1}
        ln -sf --relative $r2 {output.r2}
        """

include: "atac.smk"
