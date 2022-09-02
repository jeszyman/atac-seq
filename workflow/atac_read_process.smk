# Uses flexbar to trim and quality-filter fastq reads
rule trim:
    container:
        atac_container,
    input:
        r1 = atac_fastq_raw + "/{library}_R1.fastq.gz",
        r2 = atac_fastq_raw + "/{library}_R2.fastq.gz",
    log:
        config["logdir"] + "/{library}_atac_trim.log",
    output:
        atac_fastq_proc + "/{library}_flex_1.fastq.gz",
        atac_fastq_proc + "/{library}_flex_2.fastq.gz",
    params:
        outdir = atac_fastq_proc,
        script = config["scriptdir"]["atac"] + "/trim.sh",
        threads = config["threads"],
    resources:
        mem_mb=5000,
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Make bowtie2 index
rule atac_index:
    input:
        config["fasta"]
    params:
        base = atac_bowtie2_index,
        script = config["scriptdir"]["atac"] + "/index.sh",
    output:
        directory(atac_bowtie2_dir),
    log:
        config["logdir"] + "/atac_index.log",
    container:
        atac_container,
    shell:
        """
        {params.script} \
        {input} \
        {params.base} \
        {output} &> {log}
        """

rule align_bt2:
    container:
        atac_container,
    input:
        r1 = atac_fastq_proc + "/{library}_flex_1.fastq.gz",
        r2 = atac_fastq_proc + "/{library}_flex_2.fastq.gz",
    log:
        config["logdir"] + "/{library}_align_bt2.log",
    params:
        prefix = atac_bowtie2_index,
        script = config["scriptdir"]["atac"] + "/align_bt2.sh",
        threads = config["threads"],
    output:
        atac_bam_raw + "/{library}.bam",
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {params.prefix} \
        {params.threads} \
        {output}
        """
