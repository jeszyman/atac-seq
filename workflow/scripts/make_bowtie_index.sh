cp /mnt/ris/jschwarz/cardiac-radiobiology/ref/keep.bed resources/keep.bed
#########1#########2#########3#########4#########5#########6#########7#########8

# Function
make_bt2_index(){
    index_dir=$(dirname $3)
    mkdir -p $index_dir
    bowtie2-build -f \
                  --threads $1 \
                  $2 \
                  $3
}

# Snakemake variables
input_fa="$1"
params_prefix="$2"
params_threads="$3"

# Run
make_bt2_index $params_threads $input_fa $params_prefix
