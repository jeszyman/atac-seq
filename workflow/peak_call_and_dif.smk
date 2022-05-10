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
