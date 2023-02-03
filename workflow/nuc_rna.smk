

rule all:
    input:
        config["data_dir"] + "/ms_nuc_rna/nuc_deseq.rdata",
        config["data_dir"] + "/ms_nuc_rna/nuc_full_deseq.rdata",
        nuc_res = config["data_dir"] + "/ms_nuc_rna/nuc_res.csv",

rule make_txdb:
    output:
        txdb = config["data_dir"] + "/ref/ucsc_mm10_txdb.rds"
    script:
        "scripts/make_txdb.R"

rule make_counts:
    input:
        txdb = config["data_dir"] + "/ref/ucsc_mm10_txdb.rds",
        salmon_5469 = config["data_dir"] + "/inputs/Rentschler_s5469_MGI2048",
        salmon_5708 = config["data_dir"] + "/inputs/Rentschler_s5708_MGI2548",
    output:
        nuc_gene_cnts = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds"
    script:
        "scripts/make_counts.R"

rule make_full_nuc_deseq2_object:
    input:
        data_model = config["data_dir"] + "/data_model/data_model.RData",
        txi = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds",
    params:
        script = config["repo"] + "/workflow/scripts/make_full_nuc_deseq2_object.R"
    output:
        nuc_deseq = config["data_dir"] + "/ms_nuc_rna/nuc_full_deseq.rdata",
    log:
        config["data_dir"] + "/logs/make_full_nuc_deseq2_object.log"
    shell:
        """
        Rscript {params.script} \
	{input.data_model} \
	{input.txi} \
	{output.nuc_deseq} \
        >& {log}
        """

rule make_nuc_deseq2_object:
    input:
        data_model = config["data_dir"] + "/data_model/data_model.RData",
        txi = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds",
    params:
        script = config["repo"] + "/workflow/scripts/make_nuc_deseq2_object.R"
    output:
        nuc_deseq = config["data_dir"] + "/ms_nuc_rna/nuc_deseq.rdata",
        nuc_res = config["data_dir"] + "/ms_nuc_rna/nuc_res.csv",
    log:
        config["data_dir"] + "/logs/make_nuc_deseq2_object.log"
    shell:
        """
        Rscript {params.script} \
	{input.data_model} \
	{input.txi} \
	{output.nuc_deseq} \
	{output.nuc_res}
        >& {log}
        """
