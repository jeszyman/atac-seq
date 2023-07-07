# Snakemake variables
input_r1="$1"
input_r2="$2"
params_prefix="$3"
params_threads="$4"
output_bam="$5"


# Function
bt2_align(){
    bowtie2 --maxins 2000 --threads $4 --very-sensitive --mm -x $3 -1 $1 -2 $2 |
        samtools view -@ 8 -f 2 -F 524 -q 40 -b -o - - |
        samtools sort -@ 8 -o $5 -
    samtools index -@ 8 $5
}


# Run
bt2_align $input_r1 $input_r2 $params_prefix $params_threads $output_bam

#########1#########2#########3#########4#########5#########6#########7#########8

# Function
bt2_align(){
    bowtie2 --maxins 2000 --threads $1 --very-sensitive -x $2 -1 $3 -2 $4 | samtools view -bS - > $5
}

# Snakemake variables
input_r1="$1"
input_r2="$2"
params_prefix="$3"
params_threads="$4"
output_bam="$5"

# Run
bt2_align "$params_threads" "$params_prefix" "$input_r1" "$input_r2" "$output_bam"
