rule read_trim:
    input:
        r1 = config["fastq_dir"] + "/{library_id}_R1.fastq.gz",
        r2 = config["fastq_dir"] + "/{library_id}_R2.fastq.gz",
    params:
        outdir = config["fastq_dir"],
        threads = config["threads"],
    output:
        config["fastq_dir"] + "/{library_id}_flex_1.fastq.gz",
        config["fastq_dir"] + "/{library_id}_flex_2.fastq.gz",
    resources: 
        mem_mb=5000
    shell:
        """
        workflow/scripts/read_trim.sh {input.r1} {input.r2} {params.outdir} {params.threads}
        """

rule align_bt2:
    input:
        r1 = config["fastq_dir"] + "/{library_id}_flex_1.fastq.gz",
        r2 = config["fastq_dir"] + "/{library_id}_flex_2.fastq.gz",	
    params:
        prefix = config["bowtie_prefix"],
        threads = config["threads"],
    output:
        bam = config["bam_dir"] + "/{library_id}.bam",
    shell:
        """
        workflow/scripts/align_bt2.sh {input.r1} {input.r2} {params.prefix} {params.threads} {output.bam}
        """

rule filter_and_dedup:
    input:
        bam = config["bam_dir"] + "/{library_id}.bam",
    output:
        dedup_bam = config["bam_dir"] + "/{library_id}_dedup.bam",
        qfilt_bam = temp(config["bam_dir"] + "/{library_id}_qfilt.bam"),
        regfilt_bam = config["bam_dir"] + "/{library_id}_regfilt.bam",
        regfilt_index = config["bam_dir"] + "/{library_id}_regfilt.bam.bai",
    resources: 
        mem_mb=5000
    shell:
        """
        workflow/scripts/filter_and_dedup.sh {input.bam} \
	                                     {config[keep_bed]} \
	                                     {config[threads]} \
	                                     {output.dedup_bam} \
	                                     {output.qfilt_bam} \
	                                     {output.regfilt_bam} 
        """

rule get_open_chrom:
    input:
        regfilt_bam = config["bam_dir"] + "/{library_id}_regfilt.bam",
    output:
        unsort_open_bam = temp(config["bam_dir"] + "/{library_id}_unsort_open.bam"),
        open_bam = config["bam_dir"] + "/{library_id}_open.bam",
    shell:
        """
        workflow/scripts/get_open_chrom.sh {input.regfilt_bam} \
                                           {config[threads]} \
                                           {output.unsort_open_bam} \
                                           {output.open_bam}
        """

rule tn5_shift:
    input:
        config["bam_dir"] + "/{library_id}_regfilt.bam",
    output:
        tmp_bam = temp(config["bam_dir"] + "/{library_id}_regfilt_tmp.bam"),
        tn5_bam =      config["bam_dir"] + "/{library_id}_regfilt_tn5.bam",
    log:
        config["log_dir"] + "/tn5_shift_and_open_{library_id}_regfilt.log",
    shell:
        """
        workflow/scripts/tn5_shift.sh {input} \
	                              {config[threads]} \
	                              {output.tmp_bam} \
                                      {output.tn5_bam} > {log} 2>&1
        """

rule tn5_shift_open:
    input:
        config["bam_dir"] + "/{library_id}_open.bam",
    output:
        tmp_bam = temp(config["bam_dir"] + "/{library_id}_open_tmp.bam"),
        tn5_bam =      config["bam_dir"] + "/{library_id}_open_tn5.bam",
    log:
        config["log_dir"] + "/tn5_shift_and_open_{library_id}_open.log",
    shell:
        """
        workflow/scripts/tn5_shift.sh {input} \
	                              {config[threads]} \
	                              {output.tmp_bam} \
                                      {output.tn5_bam} > {log} 2>&1
        """

rule fastqc:
    input:
        raw = config["fastq_dir"] + "/{library_id}_{read}.fastq.gz",
    output:
        raw_html = config["qc_dir"] + "/{library_id}_{read}_fastqc.html",
    log:
        raw = config["log_dir"] + "/fastqc_raw_{library_id}_{read}.log",	
    shell:
        """
        fastqc --outdir {config[qc_dir]} \
        --quiet \
        --threads {config[threads]} {input.raw} &> {log}
        """

rule samstats:
    input:
        bam = config["bam_dir"] + "/{library_id}.bam",
    output:
        stat = config["qc_dir"] + "/{library_id}_stat.txt",
        flagstat = config["qc_dir"] + "/{library_id}_flagstat.txt",
    log:
        config["log_dir"] + "/{library_id}_samstats.log",
    shell:
        """
        workflow/scripts/samstats.sh {config[threads]} {input.bam} {output.stat} {output.flagstat} 2>&1 >> {log}
        """
