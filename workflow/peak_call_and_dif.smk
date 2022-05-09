rule make_peak_counts:
    params:
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_regfilt_tn5.bam$",	
        lib_str = config["IR48H_V_SHAM"],	
        script = config["repo"] + "/workflow/scripts/make_peak_counts.R",
    output:
        background_counts = config["data_dir"] + "/atac/background_counts_all_rse.rds"
        counts_rse = config["data_dir"] + "/atac/counts_all_rse.rds"
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
