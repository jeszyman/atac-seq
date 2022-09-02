rule all:
    input:

container: config["container"]

MACS_BROAD_EXT = ["peaks.broadPeak", "peaks.gappedPeak", "peaks.xls"]

MACS_NARROW_EXT = ["peaks.narrowPeak", "summits.bed"]

LIBRARY_IDS = ["atac1","atac2","atac3","atac4"]

BAM_PROCESS = ["regfilt", "open"]

rule all:
    input:
        expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_broad}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_broad = MACS_BROAD_EXT),
        expand(config["data_dir"] + "/macs2/{library_id}_{bam_process}_{macs_narrow}", library_id = LIBRARY_IDS, bam_process = ["open", "regfilt"], macs_narrow = MACS_NARROW_EXT),
        expand(config["data_dir"] + "/csaw/background_counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),
        expand(config["data_dir"] + "/csaw/counts_all_{bam_process}_rse.rds", bam_process = BAM_PROCESS),

include: "peak_call_and_dif.smk"
