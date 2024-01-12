container: config["container"]
RUNSAMPLES =  ["lib001", "lib002", "lib003", "lib004", "lib005", "lib006", "lib007", "lib008", "lib009", "lib010", "lib011", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib024", "lib025"]

rule all:
    input:
        expand(config["data_dir"] + "/atac/bam/{library_id}.bam", library_id=RUNSAMPLES),
        config["data_dir"] + "/ref/keep.bed",
        expand(config["data_dir"] + "/atac/bam/{library_id}_regfilt.bam", library_id=RUNSAMPLES),
        expand(config["data_dir"] + "/atac/bam/{library_id}_open.bam", library_id=RUNSAMPLES),
        expand(config["data_dir"] + "/atac/bam/{library_id}_regfilt_tn5.bam", library_id=RUNSAMPLES),
        expand(config["data_dir"] + "/atac/bam/{library_id}_open_tn5.bam", library_id=RUNSAMPLES),


# - Snakemake

rule read_trim:
    input:
        r1 = config["data_dir"] + "/atac/fastq/{library_id}_R1.fastq.gz",
        r2 = config["data_dir"] + "/atac/fastq/{library_id}_R2.fastq.gz",
    params:
        outdir = config["data_dir"] + "/atac/fastq",
        threads = config["threads"],
    output:
        config["data_dir"] + "/atac/fastq/{library_id}_flex_1.fastq.gz",
        config["data_dir"] + "/atac/fastq/{library_id}_flex_2.fastq.gz",
    resources:
        mem_mb=5000
    shell:
        """
        workflow/scripts/read_trim.sh {input.r1} {input.r2} {params.outdir} {params.threads}
        """


# - Snakemake

rule make_bowtie_index:
    input:
        fa = config["data_dir"] + "/ref/mm10.fa",
    params:
        prefix = config["data_dir"] + "/ref/ucsc_mm10_bt2/ucsc_mm10_bt2",
        threads = config["threads"]
    output:
        config["data_dir"] + "/ref/ucsc_mm10_bt2/ucsc_mm10_bt2.1.bt2",
    shell:
        """
        workflow/scripts/make_bowtie_index.sh {input.fa} {params.prefix} {params.threads}
        """


# :LOGBOOK:
# - State "WAITING"    from "TODO"       [2021-12-23 Thu 12:41]
# :END:
# - Snakemake

rule align_bt2:
    input:
        r1 = config["data_dir"] + "/atac/fastq/{library_id}_flex_1.fastq.gz",
        r2 = config["data_dir"] + "/atac/fastq/{library_id}_flex_2.fastq.gz",
    params:
        prefix = config["data_dir"] + "/ref/ucsc_mm10_bt2/ucsc_mm10_bt2",
        threads = config["threads"],
    output:
        bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
    shell:
        """
        workflow/scripts/align_bt2.sh {input.r1} {input.r2} {params.prefix} {params.threads} {output.bam}
        """
