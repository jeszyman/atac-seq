rule macs2_narrow:
    container: atac_container,
    input: atac_bam_dir + "/{library}_filt.bam",
    log: log_dir + "/{library}_macs2_narrow.log",
    output: atac_macs2_dir + "/{library}_multi_peaks.narrowPeak",
    params:
        gsize = config["gsize"],
        outdir = atac_macs2_dir,
        script = atac_script_dir + "/macs2_narrow.sh",
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
    input: atac_bam_dir + "/{library}_filt.bam",
    log: log_dir + "/{library}_macs2_broad.log",
    output: atac_macs2_dir + "/{library}_peaks.broadPeak",
    params:
        gsize = config["gsize"],
        outdir = atac_macs2_dir,
        script = atac_script_dir + "/macs2_broad.sh",
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
    input: atac_bam_dir + "/{library}_filt.bam",
    log: log_dir + "/{library}_macs2_single_summit.log",
    output: atac_macs2_dir + "/{library}_single_peaks.narrowPeak",
    params:
        gsize = config["gsize"],
        outdir = atac_macs2_dir,
        script  = atac_script_dir + "/run_macs2_corces_onesummit.sh",
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
        libraries = sample_sheet,
    output: atac_macs2_dir + "/{group}_consensus.bed",
    params:
        log_dir = log_dir,
        atac_macs2_dir = atac_macs2_dir,
        script = atac_script_dir + "/macs2_consensus.R",
    shell:
        """
        Rscript {params.script} \
        {input.libraries} \
        {params.atac_macs2_dir} \
        > {params.log_dir}/macs2_consensus.log 2>&1
        """

#
rule naive_overlap :
    container: atac_container,
    input:
        peaks = atac_macs2_dir + "/{library}_single_peaks.narrowPeak",
        sample_sheet = sample_sheet,
    log: log_dir + "/{library}_naive_overlap.log",
    output: atac_macs2_dir + "/{library}_naive.bed",
    params:
        atac_macs2_dir = atac_macs2_dir,
        script = atac_script_dir + "/naive_overlap.R",
    shell:
        """
        Rscript {params.script} \
        {input.peaks} \
        {input.sample_sheet} \
        {params.atac_macs2_dir} \
        {output} \
        > {log} 2>&1
        """

rule peak_annotation:
    input:
        config["data_dir"] + "/atac/dca.rds"
    params:
        script = config["repo"] + "/workflow/scripts/peak_annotation.R"
    output:
        annotated_counts = config["data_dir"] + "/atac/annotated_counts.rds",
    log:
        config["data_dir"] + "/logs/peak_annotation.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.annot} \
        >& {log}
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
