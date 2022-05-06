IDS, = glob_wildcards(config["data_dir"] + "/atac/fastq/{id}.fastq.gz")
RUNSAMPLES =  ["lib001", "lib002", "lib003", "lib004", "lib005", "lib006", "lib007", "lib008", "lib009", "lib010", "lib011", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib024", "lib025"]

rule all:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_stat.txt", library_id = RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id = RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_libcomplex.rds", library_id = RUNSAMPLES),
        config["data_dir"] + "/qc/frag_dist.rds",

rule fastqc:
    input:
        config["data_dir"] + "/atac/fastq/{fq_id}.fastq.gz",
    output:
        config["data_dir"] + "/qc/{fq_id}_fastqc.html",
    params:
        outdir = config["data_dir"] + "/qc",
    shell:
        """
        fastqc -t 6 -o /mnt/ris/jschwarz/cardiac-radiobiology/qc/ {input}
        """

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
        workflow/scripts/samstats.sh {config[threads]} {input.bam} {output.stat} {output.flagstat} 2>&1 >> {log}
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
