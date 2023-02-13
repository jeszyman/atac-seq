RUNSAMPLES =  ["lib001", "lib002", "lib003", "lib004", "lib005", "lib006", "lib007", "lib008", "lib009", "lib010", "lib011", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib024", "lib025"]

rule all:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_stat.txt", library_id=RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id=RUNSAMPLES),

rule samstats:
    input:
        bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
    output:
        stat = config["data_dir"] + "/qc/{library_id}_stat.txt",
        flagstat = config["data_dir"] + "/qc/{library_id}_flagstat.txt",
    log:
        config["data_dir"] + "/logs/{library_id}_samstats.log",
    shell:
        """
        "workflow/scripts/samstats.sh" {config[threads]} {input.bam} {output.stat} {output.flagstat} 2>&1 >> {log}
        """

rule fastqc:
    input:
    output:
    shell:
        """
        scripts/fastqc.sh
        """

rule multiqc:
    input:
    output:
    shell:
        """
        scripts/multiqc.sh
        """

rule atac-seq_qc:
    input:
    params:
        script = config["repo"] + "workflow/scripts/atac-seq_qc.R"
    output:
    log:
        config["data_dir"] + "/logs/atac-seq_qc.log"
    shell:
        """
        Rscript {params.script} \
        >& {log}
        """

rule library_complexity:
input:
    bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
params:
    script = config["repo"] + "/workflow/scripts/library_complexity.R",
output:
    bam = config["data_dir"] + "/qc/{library_id}_libcomplex.rds",
log:
    config["data_dir"] + "/logs/{library_id}_library_complexity.log",
shell:
    """
    Rscript {params.script} \
    {input.bam} \
    {output.bam} \
    >& {log}
    """

rule make_frag_distribution_mat:
    input:
        bam_dir = config["data_dir"] + "/atac/bam",
    params:
        script = config["repo"] + "/workflow/scripts/make_frag_distribution_mat.R",
    output:
        frag_dist = config["data_dir"] + "/qc/frag_dist.rds",
    log:
        config["data_dir"] + "/logs/make_frag_distribution_mat.log"
    shell:
        """
        Rscript {params.script} \
	{input.bam} \
	{output.frag_dist}
        >& {log}
        """
