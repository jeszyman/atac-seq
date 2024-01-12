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

rule all:
    input:
        # Run fastp on raw fastqs
        "/mnt/ris/szymanski/Active/atac_int_test/atac/fastqs/lib240_proc_R1.fastq.gz",
        #expand(f"{atac_fastq_dir}/{{library}}_raw_{{read}}.fastq.gz",
        #       library = raw_atac_libs,
        #       read = ['R1','R2']),
        expand(f"{atac_dir}/bams/{{library}}_{{build}}_raw.bam",
               library = raw_atac_libs,
               build = 'hg38'),

# ---   Include Statements   --- #
# ------------------------------ #

include: "atac.smk"
