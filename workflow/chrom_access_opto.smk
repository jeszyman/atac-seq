CHROM_FILT =  ["regfilt", "open"]
COHORT = ["sham", "ir48h", "ir6w"]
CONTRAST = ["all", "ir6w_sham", "ir48h_sham"]
JOIN = ["union", "intersect", "naive"]
IR48H_LIBS = ["lib008", "lib009", "lib010", "lib012"]
IR6W_LIBS = ["lib003", "lib004", "lib005", "lib006", "lib017", "lib019", "lib021", "lib023", "lib025"]
RUNSAMPLES =  ["lib001", "lib002", "lib003", "lib004", "lib005", "lib006", "lib007", "lib008", "lib009", "lib010", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib025"]
SHAM_LIBS = ["lib001", "lib002", "lib007", "lib013", "lib014", "lib015", "lib016", "lib018", "lib020", "lib022"]
SHAM_LIBS_FILT = ["lib013", "lib014", "lib015", "lib016", "lib018", "lib020", "lib022"]
IR6W_LIBS_FILT = ["lib017", "lib019", "lib021", "lib023", "lib025"]
WIDTH = ["broad", "narrow"]
FILTSAMPLES =  ["lib008", "lib009", "lib010", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib025"]
CALLER = ["csaw", "macs2"]

rule all:
    input:
        expand(config["data_dir"] + "/atac/bam/{cohort}_{chrom_filt}_merged_tn5.bam", cohort = COHORT, chrom_filt = CHROM_FILT),
        expand(config["data_dir"] + "/atac/bam/sham_{chrom_filt}_merged_tn5_filt.bam", chrom_filt = CHROM_FILT),
        expand(config["data_dir"] + "/atac/bam/ir6w_{chrom_filt}_merged_tn5_filt.bam", chrom_filt = CHROM_FILT),
        expand(config["data_dir"] + "/atac/macs2/{cohort}_{chrom_filt}_{width}_filt_peaks.xls", cohort = ["sham", "ir6w"], chrom_filt = CHROM_FILT, width = WIDTH)
        expand(config["data_dir"] + "/atac/macs2_consensus_beds/all_{join}_{chrom_filt}_{width}_filt.bed", join=JOIN, chrom_filt=CHROM_FILT, width=WIDTH),
        expand(config["data_dir"] + "/atac/macs2_consensus_beds/ir6w_sham_{join}_{chrom_filt}_{width}_filt.bed", join=JOIN, chrom_filt=CHROM_FILT, width=WIDTH),
        expand(config["data_dir"] + "/atac/macs2_consensus_beds/ir48h_sham_{join}_{chrom_filt}_{width}_filt.bed",  join=JOIN, chrom_filt=CHROM_FILT, width=WIDTH),
        expand(config["data_dir"] + "/atac/macs2_consensus_granges/{contrast}_{join}_{chrom_filt}_{width}_filt.rds", contrast = CONTRAST, join=JOIN, chrom_filt=CHROM_FILT, width=WIDTH),

rule make_merged_bams:
    input:
        ir48h =     expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = IR48H_LIBS),
        ir6w =      expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = IR6W_LIBS),
        ir6w_filt = expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = IR6W_LIBS_FILT),
        sham =      expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = SHAM_LIBS),
        sham_filt = expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = SHAM_LIBS_FILT),
    output:
        ir48h = config["data_dir"] + "/atac/bam/ir48h_{chrom_filt}_merged_tn5.bam",
        ir6w = config["data_dir"] + "/atac/bam/ir6w_{chrom_filt}_merged_tn5.bam",
        ir6w_filt = config["data_dir"] + "/atac/bam/ir6w_{chrom_filt}_merged_tn5_filt.bam",
        sham = config["data_dir"] + "/atac/bam/sham_{chrom_filt}_merged_tn5.bam",
        sham_filt = config["data_dir"] + "/atac/bam/sham_{chrom_filt}_merged_tn5_filt.bam",
    shell:
        """
        samtools merge -@ {config[threads]} {output.sham} {input.sham}
        samtools merge -@ {config[threads]} {output.ir6w} {input.ir6w}
        samtools merge -@ {config[threads]} {output.ir48h} {input.ir48h}
        samtools merge -@ {config[threads]} {output.sham_filt} {input.sham_filt}
        samtools merge -@ {config[threads]} {output.ir6w_filt} {input.ir6w_filt}
        """
