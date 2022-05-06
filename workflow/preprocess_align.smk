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

rule make_keep_bed:
    input:
        autosome_bed = config["data_dir"] + "/ref/grcm38_primary_assembly_chr.bed",
        blacklist_bed = config["data_dir"] + "/ref/mm10-blacklist.v2_ENSEMBL_chr.bed",
    output:
        keep_bed = config["data_dir"] + "/ref/keep.bed",
    shell:
        """
        bedtools subtract -a {input.autosome_bed} -b {input.blacklist_bed} > {output.keep_bed}
        """

rule filter_and_dedup:
    input:
        bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
    params:
        keep_bed = config["data_dir"] + "/ref/keep.bed",
        threads = config["threads"],	
    output:
        dedup_bam = config["data_dir"] + "/atac/bam/{library_id}_dedup.bam",
        qfilt_bam = temp(config["data_dir"] + "/atac/bam/{library_id}_qfilt.bam"),
        regfilt_bam = config["data_dir"] + "/atac/bam/{library_id}_regfilt.bam",
        regfilt_index = config["data_dir"] + "/atac/bam/{library_id}_regfilt.bam.bai",
    resources: 
        mem_mb=5000
    shell:
        """
        workflow/scripts/filter_and_dedup.sh {input.bam} \
	                                     {params.keep_bed} \
	                                     {params.threads} \
	                                     {output.dedup_bam} \
	                                     {output.qfilt_bam} \
	                                     {output.regfilt_bam} 
        """

rule get_open_chrom:
    input:
        regfilt_bam = config["data_dir"] + "/atac/bam/{library_id}_regfilt.bam",
    output:
        unsort_open_bam = temp(config["data_dir"] + "/atac/bam/{library_id}_unsort_open.bam"),
        open_bam = config["data_dir"] + "/atac/bam/{library_id}_open.bam",
    shell:
        """
        workflow/scripts/get_open_chrom.sh {input.regfilt_bam} \
                                           {config[threads]} \
                                           {output.unsort_open_bam} \
                                           {output.open_bam}
        """

rule tn5_shift_and_open:
    input:
        atac_bam =         config["data_dir"] + "/atac/bam/{library_id}_regfilt.bam",
    output:
        tmp_bam = temp(config["data_dir"] + "/atac/bam/{library_id}_regfilt_tmp.bam"),
        tn5_bam =      config["data_dir"] + "/atac/bam/{library_id}_regfilt_tn5.bam",
    log:
        config["data_dir"] + "/logs/tn5_shift_and_open_{library_id}_regfilt.log",
    shell:
        """
        workflow/scripts/tn5_shift.sh {input.atac_bam} \
	                              {config[threads]} \
	                              {output.tmp_bam} \
                                      {output.tn5_bam} > {log} 2>&1
        """

rule tn5_shift_open:
    input:
        atac_bam =         config["data_dir"] + "/atac/bam/{library_id}_open.bam",
    output:
        tmp_bam = temp(config["data_dir"] + "/atac/bam/{library_id}_open_tmp.bam"),
        tn5_bam =      config["data_dir"] + "/atac/bam/{library_id}_open_tn5.bam",
    log:
        config["data_dir"] + "/logs/tn5_shift_and_open_{library_id}_open.log",
    shell:
        """
        workflow/scripts/tn5_shift.sh {input.atac_bam} \
	                              {config[threads]} \
	                              {output.tmp_bam} \
                                      {output.tn5_bam} > {log} 2>&1
        """
