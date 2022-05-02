#########1#########2#########3#########4#########5#########6#########7#########8

# Function

atac_bam_processing(){
    #
    # Dedup
    samtools sort -@ $1 -n -o - $2 | 
    samtools fixmate -m - - | 
    samtools sort -@ $1 -o - - | 
    samtools markdup -@ $1 -r - $3
    #
    # Filter to aligned, properly paired reads
    samtools view -@ $1 -b -f 3 -h -o $4 $3 
    #
    # Filter to autosomes and remove blacklisted regions
    samtools view -@ $1 -b -h -L $5 -o - $4 |
    samtools sort -@ $1 -n -o - - | 
    samtools fixmate -m - - |
    samtools sort -@ $1 -o $6 -
    samtools index $6
}

# Snakemake variables
input_bam="$1"
params_keep_bed="$2"
params_threads="$3"
output_dedup_bam="$4"
output_qfilt_bam="$5"
output_regfilt_bam="$6"

# Run command
atac_bam_processing "$params_threads" \
                    "$input_bam" \
                    "$output_dedup_bam" \
                    "$output_qfilt_bam" \
                    "$params_keep_bed" \
                    "$output_regfilt_bam"
samtools index "$output_regfilt_bam"
