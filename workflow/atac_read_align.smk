# Uses flexbar to trim and quality-filter fastq reads
rule fastp:
    container: atac_container,
    input:
        r1 = atac_fastqs + "/{library}_raw_R1.fastq.gz",
        r2 = atac_fastqs + "/{library}_raw_R2.fastq.gz",
    log:
        cmd  = logdir  + "/{library}_atac_fastp.log",
        json = logdir  + "/{library}_atac_fastp.json",
        html = logdir  + "/{library}_atac_fastp.html",
    output:
        r1 = atac_fastqs + "/{library}_proc_R1.fastq.gz",
        r2 = atac_fastqs + "/{library}_proc_R2.fastq.gz",
    params:
        script  = atac_scripts + "/trim.sh",
        threads = default_threads,
    resources:
        mem_mb = 5000,
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
    input:     atac_fasta,
    log:       logdir + "/atac_index.log",
    output:
        directory(atac_bowtie2_dir),
        atac_bowtie2_index + ".1.bt2",
    params:
        base = atac_bowtie2_index,
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
        r1 = atac_fastqs + "/{library}_proc_R1.fastq.gz",
        r2 = atac_fastqs + "/{library}_proc_R2.fastq.gz",
        index = atac_bowtie2_index + ".1.bt2",
    log: logdir + "/{library}_align_bt2.log",
    params:
        prefix = atac_bowtie2_index,
        script = atac_scripts + "/align_bt2.sh",
        threads = 8,
    output:
        atac_bams + "/{library}_raw.bam",
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
    input:     atac_bams + "/{library}_raw.bam",
    log:       logdir + "/{library}_atac_dedup.log",
    output:    atac_bams + "/{library}_dedup.bam",
    params:
        script  = atac_scripts + "/dedup.sh",
	threads = default_threads,
    resources: mem_mb = 5000
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

rule make_keep_bed:
    input:
        autosome_bed = config["datadir"] + "/ref/grcm38_primary_assembly_chr.bed",
        blacklist_bed = config["datadir"] + "/inputs/mm10-blacklist.v2_ENSEMBL_chr.bed",
    output:
        keep_bed = config["datadir"] + "/ref/keep.bed",
    shell:
        """
        bedtools subtract -a {input.autosome_bed} -b {input.blacklist_bed} > {output.keep_bed}
        """

# Filter alignments by quality and reference position
rule filter_bam:
    container: atac_container,
    input: atac_bams + "/{library}_dedup.bam",
    log: logdir + "/{library}_atac_filter_bam.log",
    output: atac_bams + "/{library}_filt.bam",
    params:
        keep_bed = atac_keep_bed,
        script = atac_scripts + "/filter_bam.sh",
	threads = default_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.keep_bed} \
        {params.threads} &> {log}
        """

#
rule bampefragsize:
    container: atac_container,
    input:     expand(atac_bams + "/{library}_filt.bam", library = ATAC_LIBRARIES),
    log:       logdir + "/bampefragsize.log",
    output:
        hist = qcdir + "/deeptools_frag_lengths.png",
        raw = qcdir + "/deeptools_frag_lengths.txt",
    params:
        blacklist = atac_blacklist,
        script  = atac_scripts + "/bamPEFragmentSize_wrapper.sh",
	threads = default_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        {params.threads}
        """

rule fastqc:
    container: atac_container,
    input: atac_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_fastqc.log",
    output: qcdir + "/{library}_{processing}_{read}_fastqc.html",
    params:
        outdir = qcdir,
        script = atac_scripts + "/fastqc_wrapper.sh",
	threads = default_threads,
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
    input: atac_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat.log",
        samstat = logdir + "/{library}_{processing}_samstat.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = atac_scripts + "/alignment_qc.sh",
        threads = default_threads,
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

rule make_txdb:
    container: atac_container,
    log: logdir + "/make_txdb.log",
    output: config["datadir"] + "/ref/txdb",
    params:
        gtf = config["gtf"],
        script = atac_scripts + "/make_txdb.R",
    shell:
        """
        Rscript {params.script} \
        {params.gtf} \
        {output} \
        > {log} 2>&1
        """

rule atacseq_qc:
    container: atac_container,
    input:
        dup_bams = expand(atac_bams + "/{library}_raw.bam", library = ATAC_LIBRARIES),
        processed_bams = expand(atac_bams + "/{library}_filt.bam", library = ATAC_LIBRARIES),
        txdb = config["datadir"] + "/ref/txdb",
    log: logdir + "/atacseq_qc.log",
    output: qcdir + "/atac_qc.rdata",
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

checkpoint insert_size:
    container: atac_container,
    input: expand(atac_bams + "/{library}_filt.bam", library = ATAC_LIBRARIES),
    log: logdir + "/insert_size.log",
    output:
        tsv = qcdir + "/insert_sizes.tsv",
        plot = qcdir + "/insert_sizes.pdf",
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

rule macs2_narrow:
    container: atac_container,
    input: atac_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_macs2_narrow.log",
    output: macs2dir + "/{library}_multi_peaks.narrowPeak",
    params:
        gsize = config["gsize"],
        outdir = macs2dir,
        script = atac_scripts + "/macs2_narrow.sh",
    shell:
        """
        name=$(basename -s _filt.bam {input})
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

rule macs2_broad:
    container: atac_container,
    input: atac_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_macs2_broad.log",
    output: macs2dir + "/{library}_peaks.broadPeak",
    params:
        gsize = config["gsize"],
        outdir = macs2dir,
        script = atac_scripts + "/macs2_broad.sh",
    shell:
        """
        name=$(basename -s _filt.bam {input})
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

#
checkpoint macs2_single_summit:
    container: atac_container,
    input: atac_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_macs2_single_summit.log",
    output: macs2dir + "/{library}_single_peaks.narrowPeak",
    params:
        gsize = config["gsize"],
        outdir = macs2dir,
        script  = atac_scripts + "/run_macs2_corces_onesummit.sh",
    shell:
        """
        base=$(basename -s _filt.bam {input})
        name=${{base}}_single
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

rule macs2_consensus:
    container: atac_container,
    input:
        libraries = config["datadir"] + "/inputs/libraries.tsv",
    output: macs2dir + "/{group}_consensus.bed",
    params:
        groups = atac_groups,
        logdir = logdir,
        macs2dir = macs2dir,
        script = atac_scripts + "/macs2_consensus.R",
    shell:
        """
        Rscript {params.script} \
        {input.libraries} \
        "{params.groups}" \
        {params.macs2dir} \
        > {params.logdir}/macs2_consensus.log 2>&1
        """

#
rule naive_overlap :
    container: atac_container,
    input:
        peaks = macs2dir + "/{library}_single_peaks.narrowPeak",
        samplesheet = samplesheet,
    log: logdir + "/{library}_naive_overlap.log",
    output: macs2dir + "/{library}_naive.bed",
    params:
        macs2dir = macs2dir,
        script = atac_scripts + "/naive_overlap.R",
    shell:
        """
        Rscript {params.script} \
        {input.peaks} \
        {input.samplesheet} \
        {params.macs2dir} \
        {output} \
        > {log} 2>&1
        """
