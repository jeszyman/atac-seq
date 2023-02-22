rule diff_chrom_accessibility:
    input:
        dge = config["data_dir"] + "/csaw/dge_{bam_process}.rds",
        norm_counts_rse = config["data_dir"] + "/csaw/norm_counts_rse_{bam_process}.rds",
    params:
        groups_str = "ir48h ir48h sham sham",
        contrast = "ir48h-sham",
        script = config["atac_script_dir_dir"] + "/diff_chrom_accessibility.R",
    output:
        config["data_dir"] + "/dca/dca_granges_{bam_process}.rds"
    log:
        config["log_dir"] + "/diff_chrom_accessibility_{bam_process}.log",
    shell:
        """
        Rscript {params.script} \
        {input.dge} \
        "{params.groups_str}" \
        "{params.contrast}" \
        {input.norm_counts_rse} \
        {output} \
        >& {log}
        """

rule peak_annotation:
    input:
        config["data_dir"] + "/dca/dca_granges_{bam_process}.rds",
    params:
        script = config["atac_script_dir_dir"] + "/peak_annotation.R"
    output:
        csv = config["data_dir"] + "/dca/{bam_process}_dca.csv",
        chipseek = config["data_dir"] + "/dca/{bam_process}_chipseek.rds"
    log:
        config["data_dir"] + "/logs/{bam_process}_peak_annotation.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.csv} \
        {output.chipseek} \
        >& {log}
        """
