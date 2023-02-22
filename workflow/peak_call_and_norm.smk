
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
        script = config["atac_script_dir_dir"] + "/select_window_size.R",
        groups_str = "ir48h ir48h sham sham"
    output:
        norm_counts_rse = config["data_dir"] + "/csaw/norm_counts_rse_{bam_process}.rds",
        dge = config["data_dir"] + "/csaw/dge_{bam_process}.rds",
    log:
        config["log_dir"] + "/make_peak_counts_{bam_process}.log",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {config[threads]} \
        {output.norm_counts_rse} \
        {output.dge} \
        {params.groups_str} \
        >& {log}
        """

rule make_merged_bams:
    input:
        ir48h =     expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = IR48H_LIBS),
        sham =      expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = SHAM_LIBS),
    output:
        ir48h = config["data_dir"] + "/atac/bam/ir48h_{chrom_filt}_merged_tn5.bam",
        sham = config["data_dir"] + "/atac/bam/sham_{chrom_filt}_merged_tn5.bam",
    shell:
        """
        samtools merge -@ {config[threads]} {output.sham} {input.sham}
        samtools merge -@ {config[threads]} {output.ir48h} {input.ir48h}
        """

rule call_macs2:
    input:
        config["data_dir"] + "/atac/bam/{library_id}_{chrom_filt}_tn5.bam",
    params:
        outdir = config["data_dir"] + "/atac/macs2"
    output:
        config["data_dir"] + "/atac/macs2/{library_id}_{chrom_filt}_{width}_peaks.xls",
    shell:
        """
        macs2 callpeak --treatment {input} \
              --bdg \
              --call-summits \
              --extsize 150 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.library_id}_{wildcards.chrom_filt}_narrow \
              --nolambda \
              --nomodel \
              --outdir {params.outdir} \
              --SPMR
        #
        macs2 callpeak --treatment {input} \
              --broad \
              --broad-cutoff 0.05 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.library_id}_{wildcards.chrom_filt}_broad \
              --outdir {params.outdir}
        """

rule call_macs2_merged:
    input:
        config["data_dir"] + "/atac/bam/{cohort}_{chrom_filt}_merged_tn5.bam",
    params:
        outdir = config["data_dir"] + "/atac/macs2"
    output:
        config["data_dir"] + "/atac/macs2/{cohort}_{chrom_filt}_{width}_peaks.xls",
    shell:
        """
        macs2 callpeak --treatment {input} \
              --bdg \
              --call-summits \
              --extsize 150 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.cohort}_{wildcards.chrom_filt}_narrow \
              --nolambda \
              --nomodel \
              --outdir {params.outdir} \
              --SPMR
        #
        macs2 callpeak --treatment {input} \
              --broad \
              --broad-cutoff 0.05 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.cohort}_{wildcards.chrom_filt}_broad \
              --outdir {params.outdir}
        """

rule make_macs2_union_consensus_peaks:
    input:
        sham = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = SHAM_LIBS),
        ir48h = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = IR48H_LIBS),
    output:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/union_sham_{chrom_filt}_{width}.bed",
        ir48h = config["data_dir"] + "/atac/macs2_consensus_beds/union_ir48h_{chrom_filt}_{width}.bed",
    shell:
        """
        bedops -m {input.sham} > {output.sham}
        bedops -m {input.ir48h} > {output.ir48h}
        """

rule make_macs2_naive_consensus_peaks:
    input:
        sham = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = SHAM_LIBS),
        ir48h = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = IR48H_LIBS),
        sham_merge = config["data_dir"] + "/atac/macs2/sham_{chrom_filt}_{width}_peaks.{width}Peak",
        ir48h_merge = config["data_dir"] + "/atac/macs2/ir48h_{chrom_filt}_{width}_peaks.{width}Peak",
    output:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/naive_sham_{chrom_filt}_{width}.bed",
        ir48h = config["data_dir"] + "/atac/macs2_consensus_beds/naive_ir48h_{chrom_filt}_{width}.bed",
    shell:
        """
        bedops --element-of 50% {input.sham_merge} {input.sham} > {output.sham}
        bedops --element-of 50% {input.ir48h_merge} {input.ir48h} > {output.ir48h}
        """

rule make_cross_cohort_consensus:
    input:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/{join}_sham_{chrom_filt}_{width}.bed",
        ir48h = config["data_dir"] + "/atac/macs2_consensus_beds/{join}_ir48h_{chrom_filt}_{width}.bed",
    output:
        ir48h_sham = config["data_dir"] + "/atac/macs2_consensus_beds/ir48h_sham_{join}_{chrom_filt}_{width}.bed",
    shell:
        """
        bedops --merge {input.sham} {input.ir48h} > {output.ir48h_sham}
        """

rule bed_to_granges:
    input:
        config["data_dir"] + "/atac/macs2_consensus_beds/ir48h_sham_{join}_{chrom_filt}_{width}.bed",
    params:
        script = config["repo"] + "/workflow/scripts/bed_to_granges.R"
    output:
        config["data_dir"] + "/atac/macs2_consensus_granges/ir48h_sham_{join}_{chrom_filt}_{width}.rds",
    log:
        config["data_dir"] + "/logs/bed_to_granges_{contrast}_{join}_{chrom_filt}_{width}.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        >& {log}
        """

rule count_from_macs2_consensus:
    input:
        consensus_file = config["data_dir"] + "/atac/macs2_consensus_granges/all_{join}_{chrom_filt}_{width}.rds",
    params:
        script = config["repo"] + "/workflow/scripts/count_from_macs2_consensus.R",
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_{chrom_filt}_tn5.bam$",
    output:
        rse = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_rse.rds",
        dge = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_dge.rds",
    log:
        config["data_dir"] + "/logs/count_from_macs2_consensus_{join}_{chrom_filt}_{width}.log"
    shell:
        """
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        {input.consensus_file} \
        {output.rse} \
        {output.dge} \
        >& {log}
        """

rule call_csaw_local_peaks:
    input:
        expand(config["data_dir"] + "/atac/bam/{library_id}_{{chrom_filt}}_tn5.bam", library_id = RUNSAMPLES)
    params:
        script = config["repo"] + "/workflow/scripts/call_csaw_local_peaks.R",
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_{chrom_filt}_tn5.bam$",
    output:
        rse = config["data_dir"] + "/atac/counts/csaw_all_csaw_{chrom_filt}_csaw_peaks_rse.rds",
        dge = config["data_dir"] + "/atac/counts/csaw_all_csaw_{chrom_filt}_csaw_peaks_dge.rds",
    log:
        config["data_dir"] + "/logs/call_csaw_local_peaks_{chrom_filt}.log"
    shell:
        """
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        {output.rse} \
        {output.dge} \
        >& {log}
        """

rule make_background_bins:
    params:
        script = config["repo"] + "/workflow/scripts/make_background_bins.R",
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_{chrom_filt}_tn5.bam$",
    output:
        config["data_dir"] + "/atac/counts/bkbin_{chrom_filt}_rse.rds"
    log:
        config["data_dir"] + "/logs/make_background_bins_{chrom_filt}.log"
    shell:
        """
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        {output} \
        >& {log}
        """

rule normalize:
    input:
        counts = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_rse.rds",
        bk = config["data_dir"] + "/atac/counts/bkbin_{chrom_filt}_rse.rds",
    params:
        script = config["repo"] + "/workflow/scripts/normalize.R"
    output:
        tmm = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_tmm_rse.rds",
        loess = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_loess_rse.rds",
    log:
        config["data_dir"] + "/logs/normalize_macs2_{join}_{chrom_filt}_{width}_norm.log"
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.bk} \
        {output.tmm} \
        {output.loess} \
        >& {log}
        """

rule call_macs2_merged_filtered:
    input:
        config["data_dir"] + "/atac/bam/{cohort}_{chrom_filt}_merged_tn5_filt.bam",
    params:
        outdir = config["data_dir"] + "/atac/macs2"
    output:
        config["data_dir"] + "/atac/macs2/{cohort}_{chrom_filt}_{width}_filt_peaks.{width}Peak",
    shell:
        """
        macs2 callpeak --treatment {input} \
              --bdg \
              --call-summits \
              --extsize 150 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.cohort}_{wildcards.chrom_filt}_narrow_filt \
              --nolambda \
              --nomodel \
              --outdir {params.outdir} \
              --SPMR
        #
        macs2 callpeak --treatment {input} \
              --broad \
              --broad-cutoff 0.05 \
              --format BAMPE \
              --gsize mm \
              --keep-dup all \
              --name {wildcards.cohort}_{wildcards.chrom_filt}_broad_filt \
              --outdir {params.outdir}
        """

rule make_macs2_union_filtered_consensus_peaks:
    input:
        sham = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = SHAM_LIBS_FILT),
        ir6w = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = IR6W_LIBS_FILT),
    output:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/union_sham_{chrom_filt}_{width}_filt.bed",
        ir6w = config["data_dir"] + "/atac/macs2_consensus_beds/union_ir6w_{chrom_filt}_{width}_filt.bed",
    shell:
        """
        bedops -m {input.sham} > {output.sham}
        bedops -m {input.ir6w} > {output.ir6w}
        """

rule make_macs2_intersect_filtered_consensus_peaks:
    input:
        sham = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = SHAM_LIBS_FILT),
        ir6w = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = IR6W_LIBS_FILT),
    output:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/intersect_sham_{chrom_filt}_{width}_filt.bed",
        ir6w = config["data_dir"] + "/atac/macs2_consensus_beds/intersect_ir6w_{chrom_filt}_{width}_filt.bed",
    shell:
        """
        bedops --intersect {input.sham} > {output.sham}
        bedops --intersect {input.ir6w} > {output.ir6w}
        """

rule make_macs2_naive_filt_consensus_peaks:
    input:
        sham = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = SHAM_LIBS_FILT),
        ir6w = expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_peaks.{{width}}Peak", library_id = IR6W_LIBS_FILT),
        sham_merge = config["data_dir"] + "/atac/macs2/sham_{chrom_filt}_{width}_filt_peaks.{width}Peak",
        ir6w_merge = config["data_dir"] + "/atac/macs2/ir6w_{chrom_filt}_{width}_filt_peaks.{width}Peak",
    output:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/naive_sham_{chrom_filt}_{width}_filt.bed",
        ir6w = config["data_dir"] + "/atac/macs2_consensus_beds/naive_ir6w_{chrom_filt}_{width}_filt.bed",
    shell:
        """
        bedops --element-of 50% {input.sham_merge} {input.sham} > {output.sham}
        bedops --element-of 50% {input.ir6w_merge} {input.ir6w} > {output.ir6w}
        """

rule make_cross_cohort_filt_consensus:
    input:
        sham = config["data_dir"] + "/atac/macs2_consensus_beds/{join}_sham_{chrom_filt}_{width}_filt.bed",
        ir48h = config["data_dir"] + "/atac/macs2_consensus_beds/{join}_ir48h_{chrom_filt}_{width}.bed",
        ir6w = config["data_dir"] + "/atac/macs2_consensus_beds/{join}_ir6w_{chrom_filt}_{width}_filt.bed",
    output:
        all = config["data_dir"] + "/atac/macs2_consensus_beds/all_{join}_{chrom_filt}_{width}_filt.bed",
        ir6w_sham = config["data_dir"] + "/atac/macs2_consensus_beds/ir6w_sham_{join}_{chrom_filt}_{width}_filt.bed",
        ir48h_sham = config["data_dir"] + "/atac/macs2_consensus_beds/ir48h_sham_{join}_{chrom_filt}_{width}_filt.bed",
    shell:
        """
        bedops --merge {input.sham} {input.ir48h} {input.ir6w} > {output.all}
        bedops --merge {input.sham} {input.ir6w} > {output.ir6w_sham}
        bedops --merge {input.sham} {input.ir48h} > {output.ir48h_sham}
        """

rule bed_to_granges_filt:
    input:
        config["data_dir"] + "/atac/macs2_consensus_beds/{contrast}_{join}_{chrom_filt}_{width}_filt.bed",
    params:
        script = config["repo"] + "/workflow/scripts/bed_to_granges.R"
    output:
        config["data_dir"] + "/atac/macs2_consensus_granges/{contrast}_{join}_{chrom_filt}_{width}_filt.rds",
    log:
        config["data_dir"] + "/logs/bed_to_granges_{contrast}_{join}_{chrom_filt}_{width}_filt.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        >& {log}
        """

rule count_from_filtered_macs2_consensus:
    input:
        consensus_file = config["data_dir"] + "/atac/macs2_consensus_granges/all_{join}_{chrom_filt}_{width}_filt.rds",
    params:
        script = config["repo"] + "/workflow/scripts/count_from_macs2_consensus.R",
        bam_dir = config["data_dir"] + "/atac/bam",
        bam_pattern = "_{chrom_filt}_tn5.bam$",
    output:
        rse = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_filt_rse.rds",
        dge = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_filt_dge.rds",
    log:
        config["data_dir"] + "/logs/count_from_macs2_consensus_{join}_{chrom_filt}_{width}_filt.log"
    shell:
        """
        Rscript {params.script} \
        {params.bam_dir} \
        {params.bam_pattern} \
        {input.consensus_file} \
        {output.rse} \
        {output.dge} \
        >& {log}
        """

# Normalize each library-filtered count matrix by tmm and loess
#
rule normalize_filt:
    input:
        counts = config["data_dir"] + "/atac/counts/csaw_all_csaw_open_csaw_peaks_filt_rse.rds",
        bk =     config["data_dir"] + "/atac/counts/bkbin_{chrom_filt}_filt_rse.rds",
    params:
        script = config["repo"] + "/workflow/scripts/normalize.R"
    output:
        tmm =    config["data_dir"] + "/atac/counts/csaw_all_csaw_{chrom_filt}_filt_tmm_rse.rds",
        loess =  config["data_dir"] + "/atac/counts/csaw_all_csaw_{chrom_filt}_filt_tmm_rse.rds",
    log:
        config["data_dir"] + "/logs/normalize_filt_csaw_{chrom_filt}.log"
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.bk} \
        {output.tmm} \
        {output.loess} \
        >& {log}
        """

# Normalize each library-filtered count matrix by tmm and loess
#
rule normalize_filt:
    input:
        counts = config["data_dir"] + "/atac/counts/csaw_all_csaw_open_csaw_peaks_filt_rse.rds",
        bk =     config["data_dir"] + "/atac/counts/bkbin_{chrom_filt}_filt_rse.rds",
    params:
        script = config["repo"] + "/workflow/scripts/normalize.R"
    output:
        tmm =    config["data_dir"] + "/atac/norm/csaw_all_csaw_{chrom_filt}_filt_tmm_rse.rds",
        loess =  config["data_dir"] + "/atac/norm/csaw_all_csaw_{chrom_filt}_filt_tmm_rse.rds",
    log:
        config["data_dir"] + "/logs/normalize_filt_csaw_{chrom_filt}.log"
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.bk} \
        {output.tmm} \
        {output.loess} \
        >& {log}
        """

rule normalize:
    input:
        counts = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_peaks_rse.rds",
        bk = config["data_dir"] + "/atac/counts/bkbin_{chrom_filt}_rse.rds",
    params:
        script = config["repo"] + "/workflow/scripts/normalize.R"
    output:
        tmm = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_tmm_rse.rds",
        loess = config["data_dir"] + "/atac/counts/macs2_all_{join}_{chrom_filt}_{width}_loess_rse.rds",
    log:
        config["data_dir"] + "/logs/normalize_macs2_{join}_{chrom_filt}_{width}_norm.log"
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.bk} \
        {output.tmm} \
        {output.loess} \
        >& {log}
        """

rule tn5_shift:
    container:
        atac_container,
    input:
        atac_bam_filt + "/{library}_filt.bam",
    log:
        config["log_dir"] + "/{library}_tn5_shift.log",
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
        config["log_dir"] + "/{library}_open_chrom.log",
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

rule library_complexity:
    input:
        config["bam_dir"] + "/{library_id}.bam",
    params:
        script = config["atac_script_dir_dir"] + "/library_complexity.R",
    output:
        config["qc_dir"] + "/{library_id}_libcomplex.rds",
    log:
        config["log_dir"] + "/{library_id}_library_complexity.log",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        >& {log}
        """

rule symlink_fastqs:
    params:
        fastq = lambda w: libraries[libraries.library_id == w.library_id].fq_basename.tolist()
    output:
        r1 = config["fq_sym_dir"] + "/{library_id}_R1.fastq.gz",
        r2 = config["fq_sym_dir"] + "/{library_id}_R2.fastq.gz",
    shell:
        """
        ln -sf --relative {config[fq_src_dir]}/{params.fastq}_R1.fastq.gz {output.r1}
        ln -sf --relative {config[fq_src_dir]}/{params.fastq}_R2.fastq.gz {output.r2}
        """

rule filter_and_dedup:
    input:
        bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
    params:
        atac_keep_bed = config["data_dir"] + "/ref/keep.bed",
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
                                             {params.atac_keep_bed} \
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

rule diff_chrom_accessibility:
    input:
        dge = config["data_dir"] + "/csaw/dge_{bam_process}.rds",
        norm_counts_rse = config["data_dir"] + "/csaw/norm_counts_rse_{bam_process}.rds",
    params:
        groups_str = "ir48h ir48h sham sham",
        contrast = "ir48h-sham",
        script = config["atac_script_dir_dir"] + "/diff_chrom_accessibility.R",
    output:
        config["data_dir"] + "/dca/dca_granges_{bam_process}.rds"
    log:
        config["log_dir"] + "/diff_chrom_accessibility_{bam_process}.log",
    shell:
        """
        Rscript {params.script} \
        {input.dge} \
        "{params.groups_str}" \
        "{params.contrast}" \
        {input.norm_counts_rse} \
        {output} \
        >& {log}
        """

rule peak_annotation:
    input:
        config["data_dir"] + "/dca/dca_granges_{bam_process}.rds",
    params:
        script = config["atac_script_dir_dir"] + "/peak_annotation.R"
    output:
        csv = config["data_dir"] + "/dca/{bam_process}_dca.csv",
        chipseek = config["data_dir"] + "/dca/{bam_process}_chipseek.rds"
    log:
        config["data_dir"] + "/logs/{bam_process}_peak_annotation.log"
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.csv} \
        {output.chipseek} \
        >& {log}
        """

RUNSAMPLES =  ["lib001", "lib002", "lib003", "lib004", "lib005", "lib006", "lib007", "lib008", "lib009", "lib010", "lib011", "lib012", "lib013", "lib014", "lib015", "lib016", "lib017", "lib018", "lib019", "lib020", "lib021", "lib022", "lib023", "lib024", "lib025"]

rule all:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_stat.txt", library_id=RUNSAMPLES),
        expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id=RUNSAMPLES),

rule samstats:
    input:
        bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
    output:
        stat = config["data_dir"] + "/qc/{library_id}_stat.txt",
        flagstat = config["data_dir"] + "/qc/{library_id}_flagstat.txt",
    log:
        config["data_dir"] + "/logs/{library_id}_samstats.log",
    shell:
        """
        "workflow/scripts/samstats.sh" {config[threads]} {input.bam} {output.stat} {output.flagstat} 2>&1 >> {log}
        """

rule fastqc:
    input:
    output:
    shell:
        """
        scripts/fastqc.sh
        """

rule multiqc:
    input:
    output:
    shell:
        """
        scripts/multiqc.sh
        """

rule atac-seq_qc:
    input:
    params:
        script = config["repo"] + "workflow/scripts/atac-seq_qc.R"
    output:
    log:
        config["data_dir"] + "/logs/atac-seq_qc.log"
    shell:
        """
        Rscript {params.script} \
        >& {log}
        """

rule library_complexity:
input:
    bam = config["data_dir"] + "/atac/bam/{library_id}.bam",
params:
    script = config["repo"] + "/workflow/scripts/library_complexity.R",
output:
    bam = config["data_dir"] + "/qc/{library_id}_libcomplex.rds",
log:
    config["data_dir"] + "/logs/{library_id}_library_complexity.log",
shell:
    """
    Rscript {params.script} \
    {input.bam} \
    {output.bam} \
    >& {log}
    """

rule make_frag_distribution_mat:
    input:
        bam_dir = config["data_dir"] + "/atac/bam",
    params:
        script = config["repo"] + "/workflow/scripts/make_frag_distribution_mat.R",
    output:
        frag_dist = config["data_dir"] + "/qc/frag_dist.rds",
    log:
        config["data_dir"] + "/logs/make_frag_distribution_mat.log"
    shell:
        """
        Rscript {params.script} \
        {input.bam} \
        {output.frag_dist}
        >& {log}
        """



rule all:
    input:
        config["data_dir"] + "/ms_nuc_rna/nuc_deseq.rdata",
        config["data_dir"] + "/ms_nuc_rna/nuc_full_deseq.rdata",
        nuc_res = config["data_dir"] + "/ms_nuc_rna/nuc_res.csv",

rule make_txdb:
    output:
        txdb = config["data_dir"] + "/ref/ucsc_mm10_txdb.rds"
    script:
        "scripts/make_txdb.R"

rule make_counts:
    input:
        txdb = config["data_dir"] + "/ref/ucsc_mm10_txdb.rds",
        salmon_5469 = config["data_dir"] + "/inputs/Rentschler_s5469_MGI2048",
        salmon_5708 = config["data_dir"] + "/inputs/Rentschler_s5708_MGI2548",
    output:
        nuc_gene_cnts = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds"
    script:
        "scripts/make_counts.R"

rule make_full_nuc_deseq2_object:
    input:
        data_model = config["data_dir"] + "/data_model/data_model.RData",
        txi = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds",
    params:
        script = config["repo"] + "/workflow/scripts/make_full_nuc_deseq2_object.R"
    output:
        nuc_deseq = config["data_dir"] + "/ms_nuc_rna/nuc_full_deseq.rdata",
    log:
        config["data_dir"] + "/logs/make_full_nuc_deseq2_object.log"
    shell:
        """
        Rscript {params.script} \
        {input.data_model} \
        {input.txi} \
        {output.nuc_deseq} \
        >& {log}
        """

rule make_nuc_deseq2_object:
    input:
        data_model = config["data_dir"] + "/data_model/data_model.RData",
        txi = config["data_dir"] + "/ms_nuc_rna/nuc_gene_cnts.rds",
    params:
        script = config["repo"] + "/workflow/scripts/make_nuc_deseq2_object.R"
    output:
        nuc_deseq = config["data_dir"] + "/ms_nuc_rna/nuc_deseq.rdata",
        nuc_res = config["data_dir"] + "/ms_nuc_rna/nuc_res.csv",
    log:
        config["data_dir"] + "/logs/make_nuc_deseq2_object.log"
    shell:
        """
        Rscript {params.script} \
        {input.data_model} \
        {input.txi} \
        {output.nuc_deseq} \
        {output.nuc_res}
        >& {log}
        """

rule make_macs2_consensus_peaks:
    input:
        expand(config["data_dir"] + "/atac/macs2/{library_id}_{{chrom_filt}}_{{width}}_granges.rds")
    params:
        script = config["repo"] + "/workflow/scripts/make_macs2_consensus_peaks.R"
    output:
    log:
        config["data_dir"] + "/logs/make_macs2_consensus_peaks.log"
    shell:
        """
        Rscript {params.script} \
        >& {log}
        """
