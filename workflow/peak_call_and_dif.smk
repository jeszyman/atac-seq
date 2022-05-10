rule make_peak_counts:
    input:
        expand(config["bam_dir"] + "/{library_id}_open_tn5.bam", library_id = LIBRARY_IDS)
    params:
        script = config["atac_scripts_dir"] + "/select_window_size.R",
    output:
        background_counts = config["dca_dir"] + "/background_counts_all_rse.rds",
        counts_rse = config["dca_dir"] + "/counts_all_rse.rds",	
    log:
        config["log_dir"] + "/make_peak_counts.log",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {config[threads]} \
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
