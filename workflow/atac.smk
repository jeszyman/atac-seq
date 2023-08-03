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
    log: f"{log_dir}/{{library}}_{{build}}_{{species}}_align_bt2.log",
    params:
        prefix = f"{ref_dir}/{{build}}_bowtie2/{{build}}",
        script = atac_script_dir + "/align_bt2.sh",
        threads = 4,
    output:
        f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_raw.bam",
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
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_raw.bam",
    log: f"{log_dir}/{{library}}_{{species}}_{{build}}_atac_dedup.log",
    output: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_dedup.bam",
    params:
        script = f"{atac_script_dir}/dedup.sh",
        threads = threads,
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
        bam = f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_dedup.bam",
        keep = f"{ref_dir}/{{build}}_atac_keep.bed",
    log: f"{log_dir}/{{library}}_{{species}}_{{build}}_filter_atac_bams.log",
    output: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam",
    params:
        script = f"{atac_script_dir}/filter_atac_bams.sh",
        threads = threads
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
        duplicated_bams = expand(f"{atac_bam_dir}/{{library}}_raw.bam", library = RAW_ATAC_LIBS),
        processed_bams = expand(f"{atac_bam_dir}/{{library}}_dedup.bam", library = RAW_ATAC_LIBS),
        txdb = f"{{build}}_ensembl_txdb",
    log: f"{log_dir}/{{build}}_atacseq_qc.log",
    output: f"{qc_dir}/{{build}}_atac_qc.rdata",
    params:
        script = f"{atac_script_dir}/atacseq_qc.R",
    shell:
        """
        Rscript {params.script} \
        "{input.duplicated_bams}" \
        "{input.processed_bams}" \
        {input.txdb} \
        {output} > {log} 2>&1
        """

rule mouse_atac_multiqc:
    benchmark: f"{bench_dir}/atac_multiqc.benchmark",
    input:
        expand(f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.html",
               library = RAW_MOUSE_ATAC_LIBS,
               processing = ["raw", "proc"],
               read = ["R1","R2"]),
        expand(f"{atac_dir}/mouse/qc/{{library}}_{{processing}}_samstats.txt",
               library = RAW_MOUSE_ATAC_LIBS,
               processing = ["raw", "filt"],
               stat = ["samstats", "flagstats"]),
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

rule agg_samstat:
    input:
        lambda wildcards: expand(f"{atac_dir}/{wildcards.species}/qc/{{library}}_filt_samstats.txt", library=get_libraries(wildcards.species)),
    output:
        f"{atac_dir}/{{species}}/qc/{{species}}_atac_samstats.tsv",
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

rule atac_qc_table:
    input:
        f"{atac_dir}/{{species}}/qc/{{species}}_corces_qc.tsv",
        f"{atac_dir}/{{species}}/qc/{{species}}_frip.tsv",
        f"{atac_dir}/{{species}}/qc/{{species}}_atac_samstats.tsv",
    log: f"{log_dir}/{{species}}_atac_qc_table.log",
    output: f"{atac_dir}/{{species}}/qc/{{species}}_atac_qc.tsv",
    params:
        script = f"{atac_script_dir}/atac_qc_table.R",

rule peak_union:
    input:
        lambda wildcards: expand(f"{atac_dir}/peaks/{{library}}_{{build}}_peaks.{{peak_type}}",
                                 library = atac_map[wildcards.atac_set]['libs'],
                                 build = atac_map[wildcards.atac_set]['build'],
                                 peak_type = atac_map[wildcards.atac_set]['peak_type']),
    log: f"{log_dir}/{{atac_set}}_peak_union.bed",
    output: f"{atac_dir}/beds/{{atac_set}}_union.bed",
    params:
        script = f"{atac_script_dir}/peak_union.R",
    shell:
        """
        Rscript {params.script} "{input}" {output} >& {log}
        """

rule bamscale:
    input:
        bams = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam",
                                        build = atac_map[wildcards.atac_set]['build'],
                                        library = atac_map[wildcards.atac_set]['libs']),
        bais = lambda wildcards: expand(f"{atac_dir}/bams/{{library}}_{{build}}_filt.bam.bai",
                                        build = atac_map[wildcards.atac_set]['build'],
                                        library = atac_map[wildcards.atac_set]['libs']),
        bed = lambda wildcards: f"{atac_dir}/beds/{{atac_set}}_union.bed",
    log: f"{log_dir}/{{atac_set}}_bamscale.log",
    params:
        out_dir = f"{atac_dir}/bamscale/{{atac_set}}",
        tmp_dir = f"/tmp/{{atac_set}}",
    output:
        f"{atac_dir}/bamscale/{{atac_set}}/{{atac_set}}.FPKM_normalized_coverages.tsv",
        f"{atac_dir}/bamscale/{{atac_set}}/{{atac_set}}.raw_coverages.tsv"
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp {input.bams} {input.bais} {params.tmp_dir}

        # set the directory containing the input BAM files
        bam_dir={params.tmp_dir}

        # get a list of BAM files in the directory
        bam_files=($(ls "$bam_dir"/*.bam))

        # build the BAMscale command with the --bam flags
        bams=""
        for bam in "${{bam_files[@]}}"
        do
        bams+="--bam $bam "
        done

        BAMscale cov --bed {input.bed} \
        --prefix {wildcards.atac_set} --outdir {params.out_dir} --threads 16 \
        $bams
        rm -rf {params.tmp_dir}
        """

rule atac_ruv:
    input:
        counts = f"{atac_dir}/bamscale/{{atac_set}}/{{atac_set}}.raw_coverages.tsv",
        datmod = f"{datamodel_dir}/lists/libraries_full.rds",
    log: f"{log_dir}/{{atac_set}}_ruvk_{{ruv_k}}.log",
    output: f"{atac_dir}/ruv/{{atac_set}}_ruvk_{{ruv_k}}.rds",
    params:
        group = lambda wildcards: atac_map[wildcards.atac_set]['group'],
        ruv_k = lambda wildcards: wildcards.ruv_k,
        script = f"{atac_script_dir}/atac_ruv.R",
    shell:
        """
        Rscript {params.script} {input} "{params.group}" {params.ruv_k} {output} >& {log}
        """
