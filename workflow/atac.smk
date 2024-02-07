#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###                       ATAC-seq Snakemake File                            ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

# ---   Rules For Workflow 1: All Library Processing   --- #
# -------------------------------------------------------- #


# - Fastp
#   - Removes tail 1 nucleotide
#   - Auto-detects and removes nextera adapter

rule atac_fastp:
    input:
        r1 = f"{atac_fastq_dir}/{{library}}_raw_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_raw_R2.fastq.gz",
    log:
        cmd  = f"{log_dir}/{{library}}_atac_fastp.log",
        html = f"{log_dir}/{{library}}_atac_fastp.html",
    output:
        r1 = f"{atac_fastq_dir}/{{library}}_proc_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_proc_R2.fastq.gz",
        json = f"{atac_qc_dir}/{{library}}_atac_fastp.json",
    params:
        script  = f"{atac_script_dir}/trim.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.r2} \
        {output.json} \
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


# Align trimmed reads using bowtie2
# Alignment is per cite:corces2018, fragment lengths <= 2000, very sensitive. Alignment threads and align_load limit memory usage and avoid errors.


rule atac_align_bt2:
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



# Currently this filter just excludes unlocalized contigs. Sex chromosomes and mitochondrial reads are retained at this step.


rule make_atac_keep_bed:
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


# - Snakemake

rule atac_fastqc:
    input: f"{atac_fastq_dir}/{{library}}_{{processing}}_{{read}}.fastq.gz",
    log: f"{log_dir}/{{library}}_{{processing}}_{{read}}_fastqc.log",
    output: f"{atac_qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = atac_qc_dir,
        script = f"{atac_script_dir}/fastqc_wrapper.sh",
	threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """


# - Snakemake

rule atac_idx:
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam"
    output: f"{atac_qc_dir}/{{library}}_{{build}}_{{species}}_idxstat.txt"
    shell: "samtools idxstats {input} > {output}"


# - Snakemake

#input: f"{atac_dir}/{{species}}/bams/{{library}}_{{processing}}.bam",
rule atac_samtools_stats:
    input:
        f"{atac_bam_dir}/{{library}}_{{build}}_{{processing}}.bam",
    log: f"{log_dir}/{{library}}_{{build}}_{{processing}}_samtool_stats.log",
    output:
        stat = f"{atac_qc_dir}/{{library}}_{{build}}_{{processing}}_samstats.txt",
        flagstat = f"{atac_qc_dir}/{{library}}_{{build}}_{{processing}}_flagstat.txt",
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


# :PROPERTIES:
# :ID:       96efb30b-67c7-4df9-8c85-e2bd2fc6707f
# :END:
# Broad peak calling per [[cite:&reske2020]]
# Reference Narrow peak as in cite:corces2018 and cite:hendrickson2017


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


# :PROPERTIES:
# :ID:       957607ff-67e8-48b5-b26a-112b21e564e2
# :END:


rule atac_std_peaks:
    input:
        narrowPeak = f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.narrowPeak",
        atac_genome_bed = f"{ref_dir}/{{build}}_atac_keep.bed",
    log: f"{log_dir}/{{library}}_{{build}}_{{bam_set}}_atac_std_peaks.log",
    output:
        f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.narrowPeak_std",
    params:
        blacklist = lambda wildcards: build_map[wildcards.build]['blklist'],
        script = f"{atac_script_dir}/corces_peak_filter.py",
    shell:
        """
        python {params.script} \
        --blacklist_bed {params.blacklist} \
        --chrom_size_bed {input.atac_genome_bed} \
        --narrowPeak_file {input.narrowPeak} \
        --out_narrowPeak {output} > {log} 2>&1
        """

rule atac_multiqc:
    input:
        lambda wildcards: expand(f"{atac_qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
                                 library = atac_libs_map[wildcards.atac_group]['libs'],
                                 processing = atac_libs_map[wildcards.atac_group]['fastq_processing'],
                                 read = ['R1','R2']),
        lambda wildcards: expand(f"{atac_qc_dir}/{{library}}_{{build}}_{{processing}}_{{stats}}.txt",
                                 library = atac_libs_map[wildcards.atac_group]['libs'],
                                 build = atac_libs_map[wildcards.atac_group]['build'],
                                 processing = atac_libs_map[wildcards.atac_group]['bam_processing'],
                                 stats = ['samstats','flagstat']),
    output:
        f"{atac_qc_dir}/{{atac_group}}/{{atac_group}}_atac_multiqc.html",
    params:
        outdir = f"{atac_qc_dir}/{{atac_group}}",
        out_name = f"{{atac_group}}_atac_multiqc",
    shell:
        """
        multiqc {input} \
        --force \
        --outdir {params.outdir} \
        --filename {params.out_name}.html
        """

rule atac_agg_samstat:
    input:
        lambda wildcards: expand(f"{atac_qc_dir}/{{library}}_{{build}}_{{processing}}_samstats.txt",
                                 library = atac_libs_map[wildcards.atac_group]['libs'],
                                 build = atac_libs_map[wildcards.atac_group]['build'],
                                 processing = atac_libs_map[wildcards.atac_group]['bam_processing']),
    output:
        f"{atac_qc_dir}/{{atac_group}}/agg_samstats.tsv",
    run:
        import os
        import re

        data = []

        # Loop over the input files
        for filename in input:
            # Extract the library ID from the filename
            library_id = os.path.basename(filename).split("_")[0]

            # Open the log file
            with open(filename, "r") as f:
                lines = f.readlines()

                # Find the required lines using regular expressions
                reads = duplicated = total = mapped = error = None
                for line in lines:
                    if re.match(r"SN\traw total sequences:", line):
                        reads = int(re.search(r"\d+", line).group())
                    elif re.match(r"SN\treads duplicated:", line):
                        duplicated = int(re.search(r"\d+", line).group())
                    elif re.match(r"SN\ttotal length:", line):
                        total = int(re.search(r"\d+", line).group())
                    elif re.match(r"SN\tbases mapped \(cigar\):", line):
                        mapped = int(re.search(r"\d+", line).group())
                    elif re.match(r"SN\terror rate:", line):
                        error = float(re.search(r"\d+\.\d+[eE][+-]\d+", line).group())

                # Append the data to the list
                data.append({
                    "reads": reads,
                    "duplicated": duplicated,
                    "total": total,
                    "mapped": mapped,
                    "error": error,
                    "library": library_id
                })

        # Write the data to a file
        header = "total_reads\tduplicated_reads\ttotal_bases\tmapped_bases\terror_rate\tlibrary\n"
        with open(output[0], "w") as f:
            f.write(header)
            for d in data:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(d["reads"], d["duplicated"], d["total"], d["mapped"], d["error"], d["library"]))

rule atac_agg_fastp:
    input:
        lambda wildcards: expand(f"{atac_qc_dir}/{{library}}_atac_fastp.json",
                                 library = atac_libs_map[wildcards.atac_group]['libs']),
    output:
        f"{atac_qc_dir}/{{atac_group}}/{{atac_group}}_agg_fastp.tsv",
    run:
        import json
        import csv

        # Define the output TSV file
        tsv_file = output[0]

        # Initialize a list to store the q30_rates and libraries
        q30_rates = []
        libraries = []

        # Loop over each input file
        for input_file in input:
            # Load the JSON file
            with open(input_file) as f:
                data = json.load(f)

                # Extract the q30_rate from the after_filtering section of the summary
                q30_rate = data["summary"]["after_filtering"]["q30_rate"]

                # Extract the library name from the input file name
                library = input_file.split("/")[-1].split("_")[0]

                # Append the q30_rate and library to the lists
                q30_rates.append(q30_rate)
                libraries.append(library)

        # Write the q30_rates and libraries lists to the output TSV file
        with open(tsv_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["library", "q30_rate"])
            for library, q30_rate in zip(libraries, q30_rates):
                writer.writerow([library, q30_rate])

rule atac_keep_libs:
    input:
        libs = libraries_full_rds,
        samstats = lambda wildcards: f"{atac_qc_dir}/{wildcards.atac_group}/{wildcards.atac_group}_agg_samstats.tsv",
    output:
        f"{atac_qc_dir}/{{atac_group}}/{{atac_group}}_keep_libs.tsv",
    params:
        mil_reads = lambda wildcards: atac_libs_map[wildcards.atac_group][ds_milreads],
        script = f"{atac_script_dir}/atac_keep_libs.R",
    shell:
        """
        Rscript {params.script} \
        --agg_samstats_tsv {input.samstats} \
        --libraries_full_rds {input.libs} \
        --mil_reads {params.mil_reads} \
        --out_tsv {output}
        """

rule make_dca_design:
    input: libraries_full_rds,
    log: f"{log_dir}/{{atac_set}}_make_dca_design.log",
    output: f"{atac_dir}/models/unadjusted/{{atac_set}}/design.rds",
    params:
        formula = lambda wildcards: atac_models_map[wildcards.atac_set]['formula'],
        libs = lambda wildcards: atac_models_map[wildcards.atac_set]['libs'],
        script = f"{atac_script_dir}/make_dca_design.R",
    shell:
        """
        Rscript {params.script} {input} "{params.formula}" "{params.libs}" \
        {output} \
        > {log} 2>&1
        """


# :PROPERTIES:
# :ID:       f0124001-2d9f-47a3-a55a-7004bc5db0ee
# :END:


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
        Rscript {params.script} \
        --in_peak_bed {input} \
        --out_anno_bed {output} \
        --txdb "{params.txdb}" \
        > {log} 2>&1
        """


# Filters MACS2 peaks by Corces normalized score per million and then sums total open genome from filtered peaks.

rule atac_global_open_genome:
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

rule corces_group_peak_filter:
    input:
        peaks = lambda wildcards: expand(f"{atac_dir}/peaks/{{library}}_{{build}}_{{bam_set}}_peaks.narrowPeak_std",
                                         library = atac_models_map[wildcards.atac_set]['libs'],
                                         build = atac_models_map[wildcards.atac_set]['build'],
                                         bam_set = atac_models_map[wildcards.atac_set]['bam_set']),
    log: f"{log_dir}/{{atac_set}}_corces_group_peak_filter.log",
    output:
        f"{atac_dir}/models/unadjusted/{{atac_set}}/corces_peaks_keep.bed",
    params:
        corces_count = lambda wildcards: atac_models_map[wildcards.atac_set]['corces_count'],
        overlap = lambda wildcards: atac_models_map[wildcards.atac_set]['corces_overlap'],
        script = f"{atac_script_dir}/corces_group_peak_filter.py",
    shell:
        """
        python {params.script} \
        --count {params.corces_count} \
        --out_peakset {output} \
        --overlap {params.overlap} \
        --peak_list {input.peaks} > {log} 2>&1
        """

#bed = lambda wildcards: f"{atac_dir}/{{atac_set}}/{{atac_set}}_union.bed",

rule bamscale:
    input:
        bams = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam",
                                        build = atac_models_map[wildcards.atac_set]['build'],
                                        library = atac_models_map[wildcards.atac_set]['libs']),
        bais = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam.bai",
                                        build = atac_models_map[wildcards.atac_set]['build'],
                                        library = atac_models_map[wildcards.atac_set]['libs']),
        bed= f"{atac_dir}/models/unadjusted/{{atac_set}}/corces_peaks_keep.bed",
    log: f"{log_dir}/{{atac_set}}_bamscale.log",
    params:
        out_dir = f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale",
        tmp_dir = f"/tmp/{{atac_set}}",
        script = f"{atac_script_dir}/bamscale.sh",
    output:
        f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale/FPKM_normalized_coverages.tsv",
        f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale/raw_coverages.tsv",
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

rule atac_pca:
    input:
        counts = f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale/raw_coverages.tsv",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
    log:
        f"{log_dir}/{{atac_set}}_atac_pca.log",
    output:
        rds = f"{atac_dir}/models/unadjusted/{{atac_set}}/pca.rds",
        pdf = f"{atac_dir}/models/unadjusted/{{atac_set}}/pca.pdf",
    params:
        formula = lambda wildcards: atac_models_map[wildcards.atac_set]['formula'],
        script = f"{atac_script_dir}/atac_pca.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        "{params.formula}" \
        {output} > {log} 2>&1
        """

rule atac_combat_batch_correction:
    input:
        counts = f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale/raw_coverages.tsv",
        design = f"{atac_dir}/models/unadjusted/{{atac_set}}/design.rds",
        libraries_full = f"{datamodel_dir}/lists/libraries_full.rds",
    log: f"{log_dir}/{{atac_set}}_atac_combat_batch_correction.log",
    output:
        f"{atac_dir}/models/combat/{{atac_set}}/dge_adjusted.rds",
    params:
        batch_var = lambda wildcards: atac_models_map[wildcards.atac_set]['batch_var'],
        covars = lambda wildcards: atac_models_map[wildcards.atac_set]['cohort'],
        out_dir = f"{atac_dir}/models/combat/{{atac_set}}",
        script = f"{atac_script_dir}/atac_combat_batch_correction.R",
    shell:
        """
        Rscript {params.script} \
        --batch_var "{params.batch_var}" \
        --counts_tsv {input.counts} \
        --covars "{params.covars}" \
        --design_rds {input.design} \
        --libraries_full_rds {input.libraries_full} \
        --out_dir {params.out_dir} > {log} 2>&1
        """



# - BAMscale counts adjusted with RUVseq per [[cite:&gontarz2020]]


rule atac_ruv_batch_correction:
    input:
        counts = f"{atac_dir}/models/unadjusted/{{atac_set}}/bamscale/raw_coverages.tsv",
        design = f"{atac_dir}/models/unadjusted/{{atac_set}}/design.rds",
        libraries_full = f"{datamodel_dir}/lists/libraries_full.rds",
    log: f"{log_dir}/{{atac_set}}_atac_ruv_batch_correction.log",
    output:
        f"{atac_dir}/models/ruv/{{atac_set}}/dge_ruv2.rds",
    params:
        formula = lambda wildcards: atac_models_map[wildcards.atac_set]['formula'],
        out_dir = f"{atac_dir}/models/ruv/{{atac_set}}",
        script = f"{atac_script_dir}/atac_ruv.R",
    shell:
        """
        Rscript {params.script} \
        --counts_tsv {input.counts} \
        --design_rds {input.design} \
        --formula "{params.formula}" \
        --libraries_full_rds {input.libraries_full} \
        --out_dir {params.out_dir} \
        >& {log}
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



# From RUVseq-adjusted BAMscale peak counts, differential chromatin accessibility is quantified with edgeR


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

rule atac_ruv_pca:
    input:
        counts = f"{atac_dir}/models/{{atac_set}}/ruv/ruv_{{ruv_k}}_counts.rds",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
    log:
        f"{log_dir}/{{atac_set}}_atac_ruv_{{ruv_k}}_pca.log",
    output:
        png = f"{atac_dir}/models/{{atac_set}}/ruv/ruv_{{ruv_k}}_pca.png",
        svg = f"{atac_dir}/models/{{atac_set}}/ruv/ruv_{{ruv_k}}_pca.svg",
    params:
        formula = lambda wildcards: atac_map[wildcards.atac_set]['formula'],
        script = f"{atac_script_dir}/atac_ruv_pca.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        "{params.formula}" \
        {output} > {log} 2>&1
        """
