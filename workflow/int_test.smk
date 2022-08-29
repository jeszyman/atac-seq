container: config["container"]

MACS_BROAD_EXT = ["peaks.broadPeak", "peaks.gappedPeak", "peaks.xls"]

MACS_NARROW_EXT = ["peaks.narrowPeak", "summits.bed"]

BAM_PROCESS = ["regfilt", "open"]

LIBRARY_IDS = ["atac1", "atac2", "atac3", "atac4"]

rule all:
    input:
        expand(config["fq_sym_dir"] + "/{library_id}_{read}.fastq.gz", library_id = LIBRARY_IDS, read = ["R1", "R2"]),
        expand(config["fq_sym_dir"] + "/{library_id}_flex_1.fastq.gz", library_id = LIBRARY_IDS),
        expand(config["fq_sym_dir"] + "/{library_id}_flex_2.fastq.gz", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_dedup.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt.bam.bai", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_open.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_regfilt_tn5.bam", library_id = LIBRARY_IDS),
        expand(config["bam_dir"] + "/{library_id}_open_tn5.bam", library_id = LIBRARY_IDS),
        expand(config["qc_dir"] + "/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1", "R2"]),
        expand(config["qc_dir"] + "/{library_id}_stat.txt", library_id = LIBRARY_IDS),
        expand(config["qc_dir"] + "/{library_id}_flagstat.txt", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_broad}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_broad = MACS_BROAD_EXT),
        expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_narrow}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_narrow = MACS_NARROW_EXT),
        expand(config["data_dir"] + "/csaw/background_counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        expand(config["data_dir"] + "/csaw/counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        expand(config["data_dir"] + "/open_chrom/{library_id}_open_chrom.txt", library_id = LIBRARY_IDS),
	#expand(config["data_dir"] + "/csaw/norm_counts_rse_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/csaw/dge_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/dca_granges_{bam_process}.rds", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/{bam_process}_dca.csv", bam_process = BAM_PROCESS),
        #expand(config["data_dir"] + "/dca/{bam_process}_chipseek.rds", bam_process = BAM_PROCESS),


include: "atac_read_process.smk"
include: "peak_call_and_norm.smk"
include: "dca_and_annot.smk"
