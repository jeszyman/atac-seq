#########1#########2#########3#########4#########5#########6#########7#########8
#
# Function for flexbar processing
flexbar_atac() {
    base=$(basename -s _R1.fastq.gz $1)
    flexbar \
        --adapter-pair-overlap ON \
        --adapter-preset Nextera \
        --pre-trim-right 1 \
        --reads "${1}" \
        --reads2 "${2}" \
        --target "${3}/${base}_flex" \
        --threads ${4} \
        --zip-output GZ
}

# Snakemake parameters
input_r1="$1"
input_r2="$2"
params_outdir="$3"
params_threads="$4"

# Run
flexbar_atac "${input_r1}" "${input_r2}" "${params_outdir}" "${params_threads}"
