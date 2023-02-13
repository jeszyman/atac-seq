# Snakemake variables
# Function
# Run command
for file in $data_dir}/atac/atac-fastq/*.fastq.gz;
do
    fastqc --outdir=$data_dir}/results/qc $file &>> "$data_dir}/log/fastqc_log.txt"
done
