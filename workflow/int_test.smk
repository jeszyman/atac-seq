container: config["container"]

LIBRARY_IDS = ["atac1","atac2","atac3","atac4"]

rule all:
    input:
        expand(config["bam_dir"] + "/{library_id}_dedup.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt.bam.bai", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_open.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt_tn5.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_open_tn5.bam", library_id = LIBRARY_IDS),
        expand(config["qc_dir"] + "/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1", "R2"]),
        expand(config["qc_dir"] + "/{library_id}_stat.txt", library_id = LIBRARY_IDS),
        expand(config["qc_dir"] + "/{library_id}_flagstat.txt", library_id = LIBRARY_IDS),

include: "atac_read_process.smk"
