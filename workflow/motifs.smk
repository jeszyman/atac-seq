<<smk_preamble>>

rule all:
    input:

rule extract_gene_list:
    input:
    params:
        script = config["repo"] + "/workflow/scripts/extract_gene_list.R"
    output:
    log:
        config["data_dir"] + "/logs/extract_gene_list.log"
    shell:
        """
        Rscript {params.script} \
        >& {log}
        """

<<smk_preamble>>

rule all:
    input:

rule extract_gene_list:
    input:
    params:
        script = config["repo"] + "/workflow/scripts/extract_gene_list.R"
    output:
    log:
        config["data_dir"] + "/logs/extract_gene_list.log"
    shell:
        """
        Rscript {params.script} \
        >& {log}
        """
