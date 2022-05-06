container: config["container"]

LIBRARY_IDS = ["atac1","atac2","atac3","atac4"]

rule all:
    input:
        expand(config["fastq_dir"] + "/{library_id}_flex_1.fastq.gz", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}.bam", library_id = LIBRARY_IDS),

include: "atac_read_process.smk"
