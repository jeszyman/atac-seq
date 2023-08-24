#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###                       ATAC-seq Snakemake File                            ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

rule fastp:
    conda: "atac",
    input:
        r1 = f"{atac_fastq_dir}/{{library}}_raw_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_raw_R2.fastq.gz",
    log:
        cmd  = f"{log_dir}/{{library}}_atac_fastp.log",
        json = f"{log_dir}/{{library}}_atac_fastp.json",
        html = f"{log_dir}/{{library}}_atac_fastp.html",
    output:
        r1 = f"{atac_fastq_dir}/{{library}}_proc_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_proc_R2.fastq.gz",
    params:
        script  = f"{atac_script_dir}/trim.sh",
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
        &> {log.cmd} && [[ -s {log.html} ]]
        """

rule atac_index:
    conda: "atac",
    input: f"{ref_dir}/{{build}}.fna.gz",
    log:   f"{log_dir}/{{build}}_atac_index.log",
    output:
        f"{ref_dir}/{{build}}_bowtie2/{{build}}.1.bt2",
    params:
        base = f"{ref_dir}/{{build}}_bowtie2/{{build}}",
        dir = f"{ref_dir}/{{build}}_bowtie2",
        script = atac_script_dir + "/index.sh",
        threads = threads
    shell:
        """
        {params.script} \
        {input} \
        {params.base} \
        {params.dir} \
        {params.threads} &> {log}
        """

rule align_bt2:
    conda: "atac",
    input:
        r1 = f"{atac_fastq_dir}/{{library}}_proc_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_proc_R2.fastq.gz",
        index = f"{ref_dir}/{{build}}_bowtie2/{{build}}.1.bt2",
    log: f"{log_dir}/{{library}}_{{build}}_align_bt2.log",
    params:
        prefix = f"{ref_dir}/{{build}}_bowtie2/{{build}}",
        script = atac_script_dir + "/align_bt2.sh",
        threads = 4,
    output:
        f"{atac_dir}/bams/{{library}}_{{build}}_raw.bam",
    resources:
        align_load = 50,
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {params.prefix} \
        {params.threads} \
        {output}
        """

rule atac_dedup:
    conda: "atac",
    input: f"{atac_dir}/bams/{{library}}_{{build}}_raw.bam",
    log: f"{log_dir}/{{library}}_{{build}}_atac_dedup.log",
    output: f"{atac_dir}/bams/{{library}}_{{build}}_dedup.bam",
    params:
        script = f"{atac_script_dir}/dedup.sh",
        threads = 4,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

rule make_atac_keep_bed:
    conda: "atac",
    input: f"{ref_dir}/{{build}}_chrome_sizes.txt",
    log: f"{log_dir}/{{build}}_make_atac_keep_bed.log",
    output: f"{ref_dir}/{{build}}_atac_keep.bed",
    params:
        script = f"{atac_script_dir}/make_atac_keep_bed.sh"
    shell:
        """
        {params.script} \
        {input} \
        {output} &> {log}
        """

rule filter_atac_bams:
    conda: "atac",
    input:
        bam = f"{atac_dir}/bams/{{library}}_{{build}}_dedup.bam",
        keep = f"{ref_dir}/{{build}}_atac_keep.bed",
    log: f"{log_dir}/{{library}}_{{build}}_filter_atac_bams.log",
    output:
        f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam",
        f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam.bai",
    params:
        script = f"{atac_script_dir}/filter_atac_bams.sh",
        threads = 4,
    shell:
        """
        {params.script} \
        {input.bam} \
        {input.keep} \
        {output} {params.threads} &> {log}
        """

rule atac_fastqc:
    conda: "atac"
    input: f"{atac_fastq_dir}/{{library}}_{{processing}}_{{read}}.fastq.gz",
    log: f"{log_dir}/{{library}}_{{processing}}_{{read}}_fastqc.log",
    output: f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = qc_dir,
        script = f"{atac_script_dir}/fastqc_wrapper.sh",
	threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

rule atac_idx:
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam"
    output: f"{qc_dir}/{{library}}_{{build}}_{{species}}_idxstat.txt"
    shell: "samtools idxstats {input} > {output}"

#input: f"{atac_dir}/{{species}}/bams/{{library}}_{{processing}}.bam",
rule samtools_stats:
    input:
        f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_{{processing}}.bam",
    log: f"{log_dir}/{{library}}_{{build}}_{{processing}}_{{species}}_samtool_stats.log",
    output:
        stat = f"{atac_dir}/{{species}}/qc/{{library}}_{{build}}_{{processing}}_samstats.txt",
        flagstat = f"{atac_dir}/{{species}}/qc/{{library}}_{{build}}_{{processing}}_flagstat.txt",
    params:
        script = f"{atac_script_dir}/samtools_stats.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {output.stat} \
        {output.flagstat} \
        {params.threads} 2>&1 >> {log}
        """

rule atacseq_qc:
    input:
        dup_bams = lambda wildcards: expand(f"{atac_bam_dir}/{{library}}_{{build}}_raw.bam",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 build = atac_map[wildcards.atac_set]['build']),
        proc_bams = lambda wildcards: expand(f"{atac_bam_dir}/{{library}}_{{build}}_dedup.bam",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 build = atac_map[wildcards.atac_set]['build']),
        txdb = f"{{build}}_ensembl_txdb",
    log: f"{log_dir}/{{build}}_atacseq_qc.log",
    output: f"{qc_dir}/{{build}}_atac_qc.rdata",
    params:
        script = f"{atac_script_dir}/atacseq_qc.R",
    shell:
        """
        Rscript {params.script} \
        "{input.dup_bams}" \
        "{input.proc_bams}" \
        {input.txdb} \
        {output} > {log} 2>&1
        """

rule atac_multiqc:
    input:
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.html",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 processing = ["raw", "proc"],
                                 read = ["R1","R2"]),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{processing}}_samstats.txt",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 processing = ["raw", "dedup", "filt"]),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{processing}}_flagstat.txt",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 processing = ["raw", "dedup", "filt"]),
    log: f"{log_dir}/atac_multiqc.log",
    output: f"{atac_dir}/mouse/qc/mouse_atac_multiqc.html",
    params:
        out_dir = f"{atac_dir}/mouse/qc",
        script = f"{atac_script_dir}/multiqc.sh",
    shell:
        """
        {params.script} \
        {input} {params.out_dir} &> {log}
        """

rule macs2:
    input:
        f"{atac_dir}/bams/{{library}}_{{build}}_{{bam_set}}.bam",
    log:
        f"{log_dir}/{{library}}_{{build}}_{{bam_set}}_macs2.log",
    output:
        broad = f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.broadPeak",
        narrow = f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.narrowPeak",
    params:
        gsize = lambda wildcards: build_map[wildcards.build]['gsize'],
        outdir = f"{atac_dir}/peaks",
        script = f"{atac_script_dir}/macs2.sh",
    shell:
        """
        name=$(basename -s .bam {input})
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

rule peak_annotation:
    input:
        f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.{{peaktype}}Peak",
    log:
        f"{log_dir}/{{library}}_{{build}}_{{bam_set}}_{{peaktype}}_peak_annotation.log",
    output:
        f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.{{peaktype}}Peak_anno.bed",
    params:
        script = f"{atac_script_dir}/peak_annotation.R",
        txdb = lambda wildcards: build_map[wildcards.build]['txdb'],
    shell:
        """
        Rscript {params.script} {input} "{params.txdb}" {output} > {log} 2>&1
        """

rule downsample_bam:
    input:
        f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam",
    log:
        f"{log_dir}/{{library}}_{{build}}_ds{{milreads}}_bam.log",
    output:
        ds = f"{atac_dir}/bams/{{library}}_{{build}}_ds{{milreads}}.bam",
        index = f"{atac_dir}/bams/{{library}}_{{build}}_ds{{milreads}}.bam.bai",
    params:
        script = f"{atac_script_dir}/downsample_bam.sh",
        threads = threads,
        milreads = 9,
    shell:
        """
        {params.script} \
        {input} {params.milreads} {params.threads} {output.ds} &> {log}
        """

rule chr_state_open_genome:
    input:
        lambda wildcards: expand(f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.{{peaktype}}_anno.bed",
                                 library=atac_map[wildcards.atac_set]['libs'],
                                 build=atac_map[wildcards.atac_set]['build'],
                                 bam_set=atac_map[wildcards.atac_set]['bam_set'],
                                 peaktype=atac_map[wildcards.atac_set]['peaktype']),
    log:
        f"{log_dir}/{{atac_set}}_{{state}}_{{qval}}_chr_state_open_genome.log",
    output:
        f"{atac_dir}/models/{{atac_set}}/open/{{state}}_q{{qval}}_open_chrom.txt"
    params:
        genome_bed=lambda wildcards: f"{ref_dir}/{atac_map[wildcards.atac_set]['build']}_sorted_autosomes.bed",
        script = f"{atac_script_dir}/chr_state_open_genome.sh",
        threads = 4
    shell:
        """
        {params.script} \
        "{input}" \
        {params.genome_bed} \
        {wildcards.state} \
        {wildcards.qval} \
        {params.threads} \
        {output} > {log} 2>&1
        """

rule make_dca_design:
    input: libraries_full_rds,
    log: f"{log_dir}/{{atac_set}}_make_dca_design.log",
    output: f"{atac_dir}/models/{{atac_set}}/design.rds",
    params:
        formula = lambda wildcards: atac_map[wildcards.atac_set]['formula'],
        libs = lambda wildcards: atac_map[wildcards.atac_set]['libs'],
        script = f"{atac_script_dir}/make_dca_design.R",
    shell:
        """
        Rscript {params.script} {input} "{params.formula}" "{params.libs}" \
        {output} \
        > {log} 2>&1
        """

rule peak_filtering:
    input:
        chrs = lambda wildcards: f"{ref_dir}/{atac_map[wildcards.atac_set]['species']}_peak_chrs.txt",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
        peaks = lambda wildcards: expand(f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.narrowPeak_anno.bed",
                                         library = atac_map[wildcards.atac_set]['libs'],
                                         build = atac_map[wildcards.atac_set]['build'],
                                         bam_set = atac_map[wildcards.atac_set]['bam_set']),
    log: f"{log_dir}/{{atac_set}}_peak_filtering.log",
    output:
        all = temp(f"{atac_dir}/models/{{atac_set}}/corces_peaks_all.bed"),
        clust = temp(f"{atac_dir}/models/{{atac_set}}/corces_peaks_clust.bed"),
        keep = f"{atac_dir}/models/{{atac_set}}/corces_peaks_keep.bed",
    params:
        corces_min = lambda wildcards: atac_map[wildcards.atac_set]['corces_min'],
        lib_peaks_min = lambda wildcards: atac_map[wildcards.atac_set]['lib_peaks_min'],
        out_dir = f"{atac_dir}/models/{{atac_set}}/",
        script = f"{atac_script_dir}/peak_filtering.R",
    shell:
        """
        mkdir -p {params.out_dir} &&
        Rscript {params.script} \
        {input.chrs} \
        {input.libs} \
        "{input.peaks}" \
        {params.corces_min} \
        {output.all} \
        {output.clust} \
        {params.lib_peaks_min} \
        {output.keep} > {log} 2>&1
        """

#bed = lambda wildcards: f"{atac_dir}/{{atac_set}}/{{atac_set}}_union.bed",

rule bamscale:
    input:
        bams = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam",
                                        build = atac_map[wildcards.atac_set]['build'],
                                        library = atac_map[wildcards.atac_set]['libs']),
        bais = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam.bai",
                                        build = atac_map[wildcards.atac_set]['build'],
                                        library = atac_map[wildcards.atac_set]['libs']),
        bed= f"{atac_dir}/models/{{atac_set}}/corces_peaks_keep.bed",
    log: f"{log_dir}/{{atac_set}}_bamscale.log",
    params:
        out_dir = f"{atac_dir}/models/{{atac_set}}/bamscale",
        tmp_dir = f"/tmp/{{atac_set}}",
        script = f"{atac_script_dir}/bamscale.sh",
    output:
        f"{atac_dir}/models/{{atac_set}}/bamscale/FPKM_normalized_coverages.tsv",
        f"{atac_dir}/models/{{atac_set}}/bamscale/raw_coverages.tsv",
    shell:
        """
        {params.script} \
        "{input.bams}" \
        "{input.bais}" \
        {input.bed} \
        {params.tmp_dir} \
        {wildcards.atac_set} \
        {params.out_dir} > {log} 2>&1
        """

rule atac_edger_fit:
    input:
        counts = f"{atac_dir}/models/{{atac_set}}/bamscale/raw_coverages.tsv",
        design = f"{atac_dir}/models/{{atac_set}}/design.rds",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
    log:
        f"{log_dir}/{{atac_set}}_atac_edger_fit.log",
    output:
        dge = f"{atac_dir}/models/{{atac_set}}/dge.rds",
        fit = f"{atac_dir}/models/{{atac_set}}/fit.rds",
    params:
        script = f"{atac_script_dir}/atac_edger_fit.R",
    shell:
        """
        Rscript {params.script} {input} {output} > {log} 2>&1
        """

rule atac_ruv:
    input:
        counts = f"{atac_dir}/models/{{atac_set}}/bamscale/raw_coverages.tsv",
        datmod = f"{datamodel_dir}/lists/libraries_full.rds",
        design = f"{atac_dir}/models/{{atac_set}}/design.rds",
    log: f"{log_dir}/{{atac_set}}_ruvk{{ruv_k}}.log",
    output:
        counts = f"{atac_dir}/models/{{atac_set}}/ruv/ruv_{{ruv_k}}_counts.rds",
        fit = f"{atac_dir}/models/{{atac_set}}/ruv/ruv_{{ruv_k}}_fit.rds",
    params:
        ruv_k = lambda wildcards: wildcards.ruv_k,
        script = f"{atac_script_dir}/atac_ruv.R",
    shell:
        """
        Rscript {params.script} {input} {params.ruv_k} {output} >& {log}
        """

rule atac_pca:
    input:
        counts = f"{atac_dir}/models/{{atac_set}}/bamscale/raw_coverages.tsv",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
    log:
        f"{log_dir}/{{atac_set}}_atac_pca.log",
    output:
        png = f"{atac_dir}/models/{{atac_set}}/pca.png",
        svg = f"{atac_dir}/models/{{atac_set}}/pca.svg",
    params:
        formula = lambda wildcards: atac_map[wildcards.atac_set]['formula'],
        script = f"{atac_script_dir}/atac_pca.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        "{params.formula}" \
        {output} > {log} 2>&1
        """

rule atac_edger_dca:
    input:
        design = lambda wildcards: dca_map[wildcards.contrast]['design'],
        fit = lambda wildcards: dca_map[wildcards.contrast]['fit'],
    log: f"{log_dir}/{{contrast}}_atac_edger_dca.log",
    output: f"{atac_dir}/contrasts/{{contrast}}/{{contrast}}.tsv",
    params:
        contrast_str = lambda wildcards: dca_map[wildcards.contrast]['contrast_str'],
        script = f"{atac_script_dir}/atac_edger_dca.R",
    shell:
        """
        Rscript {params.script} \
        {input.design} \
        {input.fit} \
        "{params.contrast_str}" \
        {output} > {log} 2>&1
        """
