rule make_txdb:
    container: atac_container,
    log: log_dir + "/make_txdb.log",
    output: config["data_dir"] + "/ref/txdb",
    params:
        gtf = config["gtf"],
        script = atac_script_dir + "/make_txdb.R",
    shell:
        """
        Rscript {params.script} \
        {params.gtf} \
        {output} \
        > {log} 2>&1
        """
