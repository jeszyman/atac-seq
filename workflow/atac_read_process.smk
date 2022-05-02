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
