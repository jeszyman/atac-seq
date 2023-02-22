
# Uses fastp to trim and quality-filter fastq reads
rule fastp:
    input:
        r1 = atac_fastq_dir + "/{library}_raw_R1.fastq.gz",
        r2 = atac_fastq_dir + "/{library}_raw_R2.fastq.gz",
    log:
        cmd  = log_dir  + "/{library}_atac_fastp.log",
        json = log_dir  + "/{library}_atac_fastp.json",
        html = log_dir  + "/{library}_atac_fastp.html",
    output:
        r1 = atac_fastq_dir + "/{library}_proc_R1.fastq.gz",
        r2 = atac_fastq_dir + "/{library}_proc_R2.fastq.gz",
    params:
        script  = atac_script_dir + "/trim.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {log.json} \
        {log.html} \
        {output.r1} \
        {output.r2} \
        {params.threads} \
        &> {log.cmd}
        """

# Make bowtie2 index
rule atac_index:
    input:     genome_fasta,
    log:       log_dir + "/atac_index.log",
    output:
        directory(bowtie2_dir),
        bowtie2_index + ".1.bt2",
    params:
        base = bowtie2_index,
        script = atac_script_dir + "/index.sh",
    shell:
        """
        {params.script} \
        {input} \
        {params.base} \
        {output} &> {log}
        """

rule align_bt2:
    input:
        r1 = atac_fastq_dir + "/{library}_proc_R1.fastq.gz",
        r2 = atac_fastq_dir + "/{library}_proc_R2.fastq.gz",
        index = bowtie2_index + ".1.bt2",
    log: log_dir + "/{library}_align_bt2.log",
    params:
        prefix = bowtie2_index,
        script = atac_script_dir + "/align_bt2.sh",
        threads = 8,
    output:
        atac_bam_dir + "/{library}_raw.bam",
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {params.prefix} \
        {params.threads} \
        {output}
        """
