rule all:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_stat.txt", library_id = RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id = RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_libcomplex.rds", library_id = RUNSAMPLES),
        config["data_dir"] + "/qc/frag_dist.rds",
