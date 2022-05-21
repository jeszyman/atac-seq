rule call_macs2_narrow:
    input:
        config["data_dir"] + "/bam/{library_id}_{bam_process}_tn5.bam",
    output:
        expand(config["data_dir"] + "/macs2/{{library_id}}_{{bam_process}}_{ext}", ext = MACS_NARROW_EXT)
    shell:
        """
        base=$(echo $(basename {input}) | sed 's/_tn5.*$//g')
        workflow/scripts/call_macs2_narrow.sh {input} ${{base}} "{config[data_dir]}/macs2"
        """

rule call_macs2_broad:
    input:
        config["data_dir"] + "/bam/{library_id}_{bam_process}_tn5.bam",
    output:
        expand(config["data_dir"] + "/macs2/{{library_id}}_{{bam_process}}_{ext}", ext = MACS_BROAD_EXT)
    shell:
        """
        base=$(echo $(basename {input}) | sed 's/_tn5.*$//g')
        workflow/scripts/call_macs2_broad.sh {input} ${{base}} "{config[data_dir]}/macs2"
        """

rule make_peak_counts:
    input:
        expand(config["bam_dir"] + "/{library_id}_{{bam_process}}_tn5.bam", library_id = LIBRARY_IDS)
    params:
        script = config["atac_scripts_dir"] + "/select_window_size.R",
    output:
        background_counts = config["data_dir"] + "/csaw/background_counts_all_{bam_process}_rse.rds",
        counts_rse = config["data_dir"] + "/csaw/counts_all_{bam_process}_rse.rds",	
    log:
        config["log_dir"] + "/make_peak_counts_{bam_process}.log",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {config[threads]} \
        {output.background_counts} \
        {output.counts_rse} \
        >& {log}
        """

rule open_genome:
    input:
        config["data_dir"] + "/bam/{library_id}_open_tn5.bam",
    params:
        genome_bed = "resources/mm10.bed",
    output:
        config["data_dir"] + "/open_chrom/{library_id}_open_chrom.txt"
    shell:
        """
        bedmap --echo --bases-uniq --delim '\t' {params.genome_bed} {input} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' > {output}

        """
