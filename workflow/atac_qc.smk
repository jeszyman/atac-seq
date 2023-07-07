rule fastqc:
    input: atac_fastq_dir + "/{library}_{processing}_{read}.fastq.gz",
    log: log_dir + "/{library}_{processing}_{read}_fastqc.log",
    output: qc_dir + "/{library}_{processing}_{read}_fastqc.html",
    params:
        outdir = qc_dir,
        script = atac_script_dir + "/fastqc_wrapper.sh",
	threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """
