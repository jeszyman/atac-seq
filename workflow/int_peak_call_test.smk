container: config["container"]
	   
LIBRARY_IDS = ["atac1","atac2","atac3","atac4"]

rule all:
    input:
        config["dca_dir"] + "/background_counts_all_rse.rds",
        config["dca_dir"] + "/counts_all_rse.rds",

include: "peak_call_and_dif.smk"