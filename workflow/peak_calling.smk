

rule all:
    input:
        config["data_dir"] + "/atac/bk_rse.rds",	
        config["data_dir"] + "/atac/counts_rse.rds",

rule make_peak_counts:
    params:
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_regfilt_tn5.bam$",	
        lib_str = {lib_str}

        expand(config["data_dir"] + "/atac/bam/{library_id}.bam", library_id=RUNSAMPLES),

        lib_str = config["IR48H_V_SHAM"],	
        script = config["repo"] + "/workflow/scripts/make_peak_counts.R",
    output:
        background_counts = config["data_dir"] + "/atac/{c}background_counts_rse.rds"
        counts_rse = config["data_dir"] + "/atac/counts_rse.rds"
    log:
        config["data_dir"] + "/logs/make_peak_counts.log",
    shell:
        """
        lib_str="{params.lib_str}"
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        "${{lib_str}}" \
        {config.threads} \
        {output.background_counts} \
        {output.counts_rse} \
        >& {log}
        """

rule make_peak_counts:
    params:
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_regfilt_tn5.bam$",	
        lib_str = config["IR48H_V_SHAM"],	
        script = config["repo"] + "/workflow/scripts/make_peak_counts.R",
    output:
        background_counts = config["data_dir"] + "/atac/background_counts_rse.rds"
        counts_rse = config["data_dir"] + "/atac/counts_rse.rds"
	window_size = config["data_dir"] + "/atac/window_size.rds",
    log:
        config["data_dir"] + "/logs/make_peak_counts.log",
    shell:
        """
        lib_str="{params.lib_str}"
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        "${{lib_str}}" \
        {config.threads} \
        {output.background_counts} \
        {output.counts_rse} \
        >& {log}
        """

rule differential_accessibility:
    input:
        background_rds = config["data_dir"] + "/atac/background_counts_rse.rds"
        counts_rds = config["data_dir"] + "/atac/counts_rse.rds",
        data_model = config["data_dir"] + "/data_model/data_model.RData",
    params:
        script = config["repo"] + "/workflow/scripts/differential_accessibility.R",
    output:
        config["data_dir"] + "/atac/dca.rds",
    log:
        config["data_dir"] + "/logs/differential_accessibility.log"
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.background} \
	{input.data_model} \
	{output}
        >& {log}
        """

rule peak_annotation:
    input:
        config["data_dir"] + "/atac/dca.rds"
    params:
        script = config["repo"] + "/workflow/scripts/peak_annotation.R"
    output:
        annotated_counts = config["data_dir"] + "/atac/annotated_counts.rds",
    log:
        config["data_dir"] + "/logs/peak_annotation.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.annot} \
        >& {log}
        """
