#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###                       ATAC-seq Snakemake File                            ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

rule fastp:
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
        &> {log.cmd}
        """

rule atac_index:
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
    input:
        r1 = f"{atac_fastq_dir}/{{library}}_proc_R1.fastq.gz",
        r2 = f"{atac_fastq_dir}/{{library}}_proc_R2.fastq.gz",
        index = f"{ref_dir}/{{build}}_bowtie2/{{build}}.1.bt2",
    log: f"{log_dir}/{{library}}_{{build}}_{{species}}_align_bt2.log",
    params:
        prefix = f"{ref_dir}/{{build}}_bowtie2/{{build}}",
        script = atac_script_dir + "/align_bt2.sh",
        threads = 8,
    output:
        f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_raw.bam",
    resources:
        load = 50,
        gpu = 1,
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

rule downsample_bam:
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam",
    log: f"{log_dir}/{{library}}_{{species}}_{{build}}_ds{{milreads}}_bam.log",
    output:
        ds = f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_ds{{milreads}}.bam",
        index = f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_ds{{milreads}}.bam.bai",
    params:
        script = f"{atac_script_dir}/downsample_bam.sh",
        threads = threads,
        milreads = 9,
    shell:
        """
        {params.script} \
        {input} {params.milreads} {params.threads} {output.ds} &> {log}
        """

rule macs2_broad:
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_{{proc}}.bam",
    log: f"{log_dir}/{{library}}_{{species}}_{{build}}_{{proc}}_peaks.broadPeak",
    output: f"{atac_dir}/{{species}}/peaks/{{library}}_{{build}}_{{proc}}_peaks.broadPeak",
    params:
        gsize = lambda wildcards: human_gsize if wildcards.species == "human" else mouse_gsize,
        outdir = f"{atac_dir}/{{species}}/peaks",
        script = f"{atac_script_dir}/macs2_broad.sh",
    shell:
        """
        name=$(basename -s .bam {input})
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

rule macs2_narrow:
    input: f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_{{proc}}.bam",
    log: f"{log_dir}/{{library}}_{{build}}_{{species}}_{{proc}}_macs2_narrow.log",
    output: f"{atac_dir}/{{species}}/peaks/{{library}}_{{build}}_{{proc}}_multi_peaks.narrowPeak",
    params:
        gsize = lambda wildcards: human_gsize if wildcards.species == "human" else mouse_gsize,
        outdir = f"{atac_dir}/{{species}}/peaks",
        script = f"{atac_script_dir}/macs2_narrow.sh",
    shell:
        """
        name=$(basename -s .bam {input})
        {params.script} \
        {input} \
        $name \
        {params.gsize} \
        {params.outdir} &> {log}
        """

rule peak_union:
    input:
        lambda wildcards: expand(f"{atac_dir}/{{species}}/peaks/{{library}}_{{build}}_peaks.broadPeak",
                                 species = wildcards.species,
                                 build = wildcards.build,
                                 library = joins[wildcards.join])
    log: f"{log_dir}/{{species}}_{{build}}_{{join}}_peak_union.bed",
    output: f"{atac_dir}/{{species}}/{{join}}_{{build}}_union.bed",
    params:
        in_dir = f"{atac_dir}/{{species}}",
        script = f"{atac_script_dir}/peak_union.R",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output} >& {log}
        """

rule bamscale_corces_mouse:
    input:
        bams = expand(f"{atac_dir}/mouse/bams/{{library}}_mm10_ds9.bam",
                      library = FILT_MOUSE_ATAC_LIBS),
        bais = expand(f"{atac_dir}/mouse/bams/{{library}}_mm10_ds9.bam.bai",
                      library = FILT_MOUSE_ATAC_LIBS),
        bed = f"{atac_dir}/mouse/peaks/mouse_corces_peaks_keep.bed"
    log: f"{log_dir}/mm10_mouse_bamscale.log",
    params:
        out_dir = f"{atac_dir}/mouse/dca",
        script = f"{atac_script_dir}/mouse_bamscale.sh",
        tmp_dir = f"/tmp/mouse_bamscale",
    output:
        f"{atac_dir}/mouse/dca/mouse_mm10.FPKM_normalized_coverages.tsv"
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp {input.bams} {input.bais} {params.tmp_dir}
        {params.script} {params.tmp_dir} {input.bed} {params.out_dir} &> {log}
        """

rule bamscale_corces_human:
    input:
        bams = expand(f"{atac_dir}/human/bams/{{library}}_hg38_ds9.bam",
                      library = FILT_HUMAN_ATAC_LIBS),
        bais = expand(f"{atac_dir}/human/bams/{{library}}_hg38_ds9.bam.bai",
                      library = FILT_HUMAN_ATAC_LIBS),
        bed = f"{atac_dir}/human/peaks/human_corces_peaks_keep.bed"
    log: f"{log_dir}/hg38_human_bamscale.log",
    params:
        out_dir = f"{atac_dir}/human/dca",
        script = f"{atac_script_dir}/human_bamscale.sh",
        tmp_dir = f"/tmp/human_bamscale",
    output:
        f"{atac_dir}/human/dca/human_hg38.FPKM_normalized_coverages.tsv"
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp {input.bams} {input.bais} {params.tmp_dir}
        {params.script} {params.tmp_dir} {input.bed} {params.out_dir} &> {log}
        """

rule peak_filtering:
    input:
        chrs = f"{datamodel_dir}/lists/{{species}}_peak_chrs.txt",
        libs = f"{datamodel_dir}/lists/libraries_full.rds",
        peaks = lambda wildcards: expand(f"{atac_dir}/{{species}}/peaks/{{library}}_{{build}}_{{proc}}_multi_peaks.narrowPeak",
                                         species=wildcards.species,
                                         library=get_species_params(wildcards.species)["libraries"],
                                         build=get_species_params(wildcards.species)["build"],
                                         proc=["ds9"]),
    log: f"{log_dir}/{{species}}_peak_filtering.log",
    output:
        all = f"{atac_dir}/{{species}}/peaks/{{species}}_corces_peaks_all.bed",
        clust = f"{atac_dir}/{{species}}/peaks/{{species}}_corces_peaks_clust.bed",
        max = f"{atac_dir}/{{species}}/peaks/{{species}}_corces_peaks_max.bed",
        keep = f"{atac_dir}/{{species}}/peaks/{{species}}_corces_peaks_keep.bed",
    params:
        corces_min = 5,
        lib_peaks_min = 1000,
        script = f"{cardiac_script_dir}/peak_filtering.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.all} \
        {output.clust} \
        {output.max} \
        {output.keep} \
        {params.corces_min} \
        {params.lib_peaks_min} > {log} 2>&1
        """

rule peak_annotation:
    input: f"{atac_dir}/{{species}}/dca/{{species}}_{{build}}.raw_coverages.tsv",
    log: f"{log_dir}/{{species}}_{{build}}_peak_annotation.log",
    output: f"{atac_dir}/{{species}}/dca/{{species}}_{{build}}_annotation.tsv",
    params:
        txdb = lambda wildcards: "TxDb.Mmusculus.UCSC.mm10.knownGene" if wildcards.species == 'mouse' else "TxDb.Hsapiens.UCSC.hg38.knownGene",
        bmart_dataset = lambda wildcards: "mmusculus_gene_ensembl" if wildcards.species == 'mouse' else "hsapiens_gene_ensembl",
        script = f"{cardiac_script_dir}/peak_annotation.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {params.bmart_dataset} {params.txdb} \
        {output} > {log} 2>&1
        """

rule make_ensembl_txdb:
    input: f"{ref_dir}/{{build}}.gtf.gz",
    output: f"{ref_dir}/{{build}}_ensembl_txdb",
    params: script = f"{atac_script_dir}/make_ensembl_txdb.R",
    shell:
        """
        Rscript {params.script} {input}
        cp /tmp/db {output}
        """

rule bamscale:
    input:
        bams = lambda wildcards: expand(f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam",
                                        species=wildcards.species,
                                        build=wildcards.build,
                                        library=joins[wildcards.join]),
        bais = lambda wildcards: expand(f"{atac_dir}/{{species}}/bams/{{library}}_{{build}}_filt.bam.bai",
                                        species=wildcards.species,
                                        build=wildcards.build,
                                        library=joins[wildcards.join]),
        bed = f"{atac_dir}/{{species}}/{{join}}_{{build}}_union.bed",
    log: f"{log_dir}/{{build}}_{{species}}_{{join}}_bamscale.log",
    params:
        out_dir = f"{atac_dir}/{{species}}/dca",
        tmp_dir = lambda wildcards: f"/tmp/{wildcards.join}",
    output:
        f"{atac_dir}/{{species}}/dca/{{join}}_{{species}}_{{build}}.FPKM_normalized_coverages.tsv"
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
        --prefix {wildcards.join}_{wildcards.species}_{wildcards.build} --outdir {params.out_dir} --threads 16 \
        $bams
        rm -rf {params.tmp_dir}
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
    input: lambda wildcards: f"{atac_dir}/{wildcards.species}/bams/{wildcards.library}_{'mm10' if wildcards.species == 'mouse' else 'hg38'}_{wildcards.processing}.bam",
    log: f"{log_dir}/{{library}}_{{processing}}_{{species}}_samtool_stats.log",
    output:
        stat = f"{atac_dir}/{{species}}/qc/{{library}}_{{processing}}_samstats.txt",
        flagstat = f"{atac_dir}/{{species}}/qc/{{library}}_{{processing}}_flagstat.txt",
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

checkpoint insert_size:
    input: expand(f"{atac_bam_dir}/{{library}}_dedup.bam", library = RAW_ATAC_LIBS),
    log: f"{log_dir}/insert_size.log",
    output:
        tsv = f"{qc_dir}/insert_sizes.tsv",
        plot = f"{qc_dir}/insert_sizes.pdf",
    params:
        peak_cut = atac_peak_cut,
        script = f"{atac_script_dir}/insert_size.R",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.tsv} \
        {output.plot} \
        {params.peak_cut} > {log} 2>&1
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

rule homer_bed_from_dca:
    input:
        dca = f"{atac_dir}/{{species}}/dca/{{species}}_atac_k{{rvu_k}}_{{contrast}}.tsv",
        anno = f"{atac_dir}/{{species}}/dca/{{species}}_annotation.tsv",
    log: f"{log_dir}/homer_bed_from_dca_{{species}}_k{{rvu_k}}_{{contrast}}.log",
    output:
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_up_all.bed",
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_down_all.bed",
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_up_promoter.bed",
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_down_promoter.bed",
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_up_enhancer.bed",
        f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_down_enhancer.bed",
    params:
        qval = 0.05,
        script = f"{cardiac_script_dir}/homer_bed_from_dca.R",
    shell:
        """
        Rscript {params.script} {input} {params.qval} {output} > log 2>&1
        """

rule homer_genome_enrich:
    input:
        bed = f"{atac_dir}/{{species}}/homer/bed/homer_{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}.bed",
        fasta = lambda wildcards: get_genome_fasta(wildcards.species),
    log: f"{log_dir}/homer_genome_enrich_{{species}}_{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}.log",
    output:
        dir = directory(f"{atac_dir}/{{species}}/homer/genome/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}"),
        known_tsv = f"{atac_dir}/{{species}}/homer/genome/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}/knownResults.txt",
        denovo_html = f"{atac_dir}/{{species}}/homer/genome/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}/homerResults.html",
    params:
        script = f"{cardiac_script_dir}/homer_genome_enrich.sh",
        threads = 4,
    shell:
       """
       {params.script} {input.bed} {input.fasta} {output.dir} {params.threads} &> {log}
       cp {output.known_tsv} $(dirname {output.known_tsv})_knownResults.txt
       """

rule homer_denovo_tsv:
    input: f"{atac_dir}/{{species}}/homer/genome/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}/homerResults.html",
    log: f"{log_dir}/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}_homer_denovo.log",
    output: f"{atac_dir}/{{species}}/homer/genome/{{species}}_k{{rvu_k}}_{{contrast}}_{{direction}}_{{set}}_homer_denovo.tsv",
    run:
        import os
        from bs4 import BeautifulSoup
        import csv

        # Define the file path
        file_path = input[0]

        # Expand the tilde to the user home directory
        file_path = os.path.expanduser(file_path)

        # Read the HTML file
        with open(file_path, 'r') as f:
            contents = f.read()

        # Parse the HTML
        soup = BeautifulSoup(contents, 'html.parser')

        # Find the table
        table = soup.find('table')

        # Find all rows
        rows = table.find_all('tr')

        # Prepare to write to TSV
        with open(output[0], 'w') as f:
            writer = csv.writer(f, delimiter='\t')

            for row in rows:
                # Find all columns
                cols = row.find_all('td')

                # Write columns to the TSV, excluding the SVG column (the second one, index 1)
                writer.writerow([col.text for i, col in enumerate(cols) if i != 1])
