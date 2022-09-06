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
	atac_bowtie2_index + ".1.bt2",
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
        index = atac_bowtie2_index + ".1.bt2",
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

# De-duplicate alignments
rule dedup:
    container:
        atac_container,
    input:
        atac_bam_raw + "/{library}.bam",
    log:
        config["logdir"] + "/{library}_atac_dedup.log",
    output:
        atac_bam_dedup + "/{library}_dedup.bam",
    params:
        script = config["scriptdir"]["atac"] + "/dedup.sh",
	threads = config["threads"],
    resources:
        mem_mb=5000
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

# Filter alignments by quality and reference position
rule filter_bam:
    container:
        atac_container,
    input:
        atac_bam_dedup + "/{library}_dedup.bam",
    log:
        config["logdir"] + "/{library}_atac_filter_bam.log",
    output:
        atac_bam_filt + "/{library}_filt.bam",
    params:
        keep_bed = atac_keep_bed,
        script = config["scriptdir"]["atac"] + "/filter_bam.sh",
	threads = config["threads"],
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.keep_bed} \
        {params.threads} &> {log}
        """

rule tn5_shift:
    container:
        atac_container,
    input:
        atac_bam_filt + "/{library}_filt.bam",
    log:
        config["logdir"] + "/{library}_tn5_shift.log",
    output:
        tmp = temp(atac_bam_tn5 + "/{library}_tn5_tmp.bam"),
        tn5 =      atac_bam_tn5 + "/{library}_tn5.bam",
    params:
        script = config["scriptdir"]["atac"] + "/tn5_shift.sh",
        threads = config["threads"],
    shell:
        """
        {params.script} \
        {input} \
        {output.tmp} \
        {output.tn5} \
        {params.threads} &> {log}
        """

rule open_chrom:
    container:
        atac_container,
    input:
        atac_bam_tn5 + "/{library}_tn5.bam",
    log:
        config["logdir"] + "/{library}_open_chrom.log",
    output:
        tmp = temp(atac_bam_open + "/{library}_open_tmp.bam"),
        open = atac_bam_open + "/{library}_open.bam",
    params:
        script = config["scriptdir"]["atac"] + "/open_chrom.sh",
        threads = config["threads"],
    shell:
        """
        {params.script} \
        {input} \
        {output.tmp} \
        {output.open} \
        {params.threads} &> {log}
        """

rule make_txdb:
    container:
        atac_container,
    log:
        config["logdir"] + "/make_txdb.log",
    output:
        config["datadir"] + "/ref/txdb",
    params:
        gtf = config["gtf"],
        script = config["scriptdir"]["atac"] + "/make_txdb.R",
    shell:
        """
        Rscript {params.script} \
        {params.gtf} \
        {output} \
        > {log} 2>&1
        """

rule atacseq_qc:
    container:
        atac_container,
    input:
        dup_bams = expand(atac_bam_raw + "/{library}.bam", library = ATAC_LIBRARIES),
        processed_bams = expand(atac_bam_tn5 + "/{library}_tn5.bam", library = ATAC_LIBRARIES),
        txdb = config["datadir"] + "/ref/txdb",
    log:
        config["logdir"] + "/atacseq_qc.log",
    output:
        config["datadir"] + "/qc/atac_qc.rdata",
    params:
        script = config["scriptdir"]["atac"] + "/atacseq_qc.R",
    shell:
        """
        Rscript {params.script} \
        "{input.dup_bams}" \
        "{input.processed_bams}" \
        {input.txdb} \
        {output} > {log} 2>&1
        """

# <DESCRIPTIVE COMMENT>
rule fastqc:
    container:
        atac_container,
    input:
        raw = atac_fastq_raw + "/{library}_R{read}.fastq.gz",
        filt = atac_fastq_proc + "/{library}_flex_{read}.fastq.gz",
    log:
        config["logdir"] + "/{library}_R{read}_atac_fastqc.log",
    output:
        raw = config["datadir"] + "/qc/{library}_R{read}_fastqc.html",
        filt = config["datadir"] + "/qc/{library}_flex_{read}_fastqc.html",
    params:
        outdir = config["datadir"] + "/qc",
    shell:
        """
        fastqc --outdir {params.outdir} --quiet --threads config[threads] {input.raw}
        fastqc --outdir {params.outdir} --quiet --threads config[threads] {input.filt}
        """

rule samstats:
    container:
        atac_container,
    input:
        atac_bam_filt + "/{library}_filt.bam",
    output:
        stat = config["datadir"] + "/qc/{library}_filt_stat.txt",
        flagstat = config["datadir"] + "/qc/{library}_filt_flagstat.txt",
    log:
        config["logdir"] + "/{library}_samstats.log",
    shell:
        """
        workflow/scripts/samstats.sh {config[threads]} {input} {output.stat} {output.flagstat} 2>&1 >> {log}
        """