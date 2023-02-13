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
