

rule all:
    input:



# extract ensembl ID lists from csaw-EdgeR DCA workflow

# - Snakemake

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
