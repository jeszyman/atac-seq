# Uses flexbar to trim and quality-filter fastq reads
rule fastp:
    container: atac_container,
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
        script  = atac_scripts + "/trim.sh",
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
    container: atac_container,
    input:     genome_fasta,
    log:       log_dir + "/atac_index.log",
    output:
        directory(bowtie2_dir),
        bowtie2_index + ".1.bt2",
    params:
        base = bowtie2_index,
        script = config["scriptdir"]["atac"] + "/index.sh",
    shell:
        """
        {params.script} \
        {input} \
        {params.base} \
        {output} &> {log}
        """

rule align_bt2:
    container: atac_container,
    input:
        r1 = atac_fastq_dir + "/{library}_proc_R1.fastq.gz",
        r2 = atac_fastq_dir + "/{library}_proc_R2.fastq.gz",
        index = bowtie2_index + ".1.bt2",
    log: log_dir + "/{library}_align_bt2.log",
    params:
        prefix = bowtie2_index,
        script = atac_scripts + "/align_bt2.sh",
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

# De-duplicate alignments
rule dedup:
    container: atac_container,
    input:     atac_bam_dir + "/{library}_raw.bam",
    log:       log_dir + "/{library}_atac_dedup.log",
    output:    atac_bam_dir + "/{library}_dedup.bam",
    params:
        script  = atac_scripts + "/dedup.sh",
	threads = threads,
    resources: mem_mb = 5000
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

rule make_atac_keep_bed:
    input:
        autosome_bed = autosome_bed,
        blacklist_bed_bed = blacklist_bed,
    output: config["data_dir"] + "/ref/atac_keep.bed",
    shell:
        """
        bedtools subtract -a {input.autosome_bed} -b {input.blacklist_bed_bed} > {output}
        """

# Filter alignments by quality and reference position
rule filter_bam:
    container: atac_container,
    input:
        bam = atac_bam_dir + "/{library}_dedup.bam",
        bed = data_dir + "/ref/atac_keep.bed",
    log: log_dir + "/{library}_atac_filter_bam.log",
    output: atac_bam_dir + "/{library}_filt.bam",
    params:
        script = atac_scripts + "/filter_bam.sh",
	threads = threads,
    shell:
        """
        {params.script} \
        {input.bam} \
        {output} \
        {input.bed} \
        {params.threads} &> {log}
        """

rule fastqc:
    container: atac_container,
    input: atac_fastq_dir + "/{library}_{processing}_{read}.fastq.gz",
    log: log_dir + "/{library}_{processing}_{read}_fastqc.log",
    output: qc_dir + "/{library}_{processing}_{read}_fastqc.html",
    params:
        outdir = qc_dir,
        script = atac_scripts + "/fastqc_wrapper.sh",
	threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Alignment samtools QC
rule alignment_qc:
    container: atac_container,
    input: atac_bam_dir + "/{library}_{processing}.bam",
    log:
        flagstat = log_dir + "/{library}_{processing}_flagstat.log",
        samstat = log_dir + "/{library}_{processing}_samstat.log",
    output:
        flagstat = qc_dir + "/{library}_{processing}_flagstat.txt",
        samstat = qc_dir + "/{library}_{processing}_samstats.txt",
    params:
        script = atac_scripts + "/alignment_qc.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {params.threads}
        """

checkpoint insert_size:
    container: atac_container,
    input: expand(atac_bam_dir + "/{library}_filt.bam", library = ATAC_LIBS),
    log: log_dir + "/insert_size.log",
    output:
        tsv = qc_dir + "/insert_sizes.tsv",
        plot = qc_dir + "/insert_sizes.pdf",
    params:
        peak_cut = atac_peak_cut,
        script = atac_scripts + "/insert_size.R",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.tsv} \
        {output.plot} \
        {params.peak_cut} > {log} 2>&1
        """

rule atacseq_qc:
    container: atac_container,
    input:
        dup_bams = expand(atac_bam_dir + "/{library}_raw.bam", library = ATAC_LIBS),
        processed_bams = expand(atac_bam_dir + "/{library}_filt.bam", library = ATAC_LIBS),
        txdb = config["data_dir"] + "/ref/txdb",
    log: log_dir + "/atacseq_qc.log",
    output: qc_dir + "/atac_qc.rdata",
    params:
        script = atac_scripts + "/atacseq_qc.R",
    shell:
        """
        Rscript {params.script} \
        "{input.dup_bams}" \
        "{input.processed_bams}" \
        {input.txdb} \
        {output} > {log} 2>&1
        """
