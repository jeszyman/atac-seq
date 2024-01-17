#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###                       Integration Testing for ATAC-seq                   ###
###                 Workflow 1: All read processing and alignment            ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

data_dir = config["data_dir"]
atac_repo = "~/repos/atac"
atac_script_dir = "scripts"
atac_dir = f"{data_dir}/atac"
atac_fastq_dir = f"{atac_dir}/fastqs"
atac_bam_dir = f"{atac_dir}/bams"
threads = 4
raw_atac_libs = ['lib240']

log_dir = f"{data_dir}/logs"
ref_dir = f"{data_dir}/ref"
atac_qc_dir = f"{data_dir}/qc/atac"
libraries_full_rds = "~/cards/data-model/lists/libraries_full.rds"
datamodel_dir = "~/cards/data-model"

#########1#########2#########3#########4#########5#########6#########7#########8


rule all:
    input:
        # Run fastp on raw fastqs                              (rule atac_fastp)
        expand(f"{atac_fastq_dir}/{{library}}_raw_{{read}}.fastq.gz",
               library = raw_atac_libs,
               read = ['R1','R2']),

        # FastQC                                              (rule atac_fastqc)
        expand(f"{atac_qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
               library = raw_atac_libs,
               processing = ['raw', 'proc'],
               read = ['R1','R2']),

        # Alignment
        ##
        ## Make index                                          (rule atac_index)
        expand(f"{ref_dir}/{{build}}_bowtie2/{{build}}.1.bt2",
               build = 'hg38'),

        ## Align, deduplicate, and filter                 (rules atac_align_bt2,
        ##                                                           atac_dedup,
        ##                                                     filter_atac_bams)
        expand(f"{atac_dir}/bams/{{library}}_{{build}}_{{processing}}.bam",
               library = raw_atac_libs,
               build = 'hg38',
               processing = ['raw', 'dedup', 'filt']),

        # Alignment QC
        ##
        ## Samtools                                   (rule atac_samtools_stats)
        expand(f"{atac_qc_dir}/{{library}}_{{build}}_{{processing}}_{{stat}}.txt",
               library = raw_atac_libs,
               build = 'hg38',
               processing = ['raw', 'dedup', 'filt'],
               stat = ['samstats','flagstat']),
        ##
        ## MultiQC                                           (rule atac_multiqc)
        #f"{atac_qc_dir}/multiqc/atac_multiqc.html"

# ---   Include Statements   --- #
# ------------------------------ #

include: "atac.smk"
